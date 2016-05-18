#!/usr/bin/env python

import sys
import math
import bisect
import argparse
import urllib.request
import shutil
import os.path
import logging
from netCDF4 import Dataset
import numpy as np
from datetime import datetime,timedelta
import json

# http://www.weathergraphics.com/identifiers/master-location-identifier-database-20130801.csv

datasets = {
    'surface-height': {
        'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/hgt.sfc.gauss.nc'
        , 'file-name': 'hgt.sfc.gauss.nc'
        , 'desc': 'Geopotential height at the surface'
    },
    'land-sea-mask' : {
        'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/land.sfc.gauss.nc'
        , 'file-name': 'land.sfc.gauss.nc'
        , 'desc': 'land-sea mask'
    },
    'clear sky surface solar radiation 4-daily' : {
        'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/csdsf.sfc.gauss.2016.nc'
        , 'file-name': 'csdsf.sfc.gauss.2016.nc'
        , 'desc': '4-daily influx of solar radiation'
    }
}

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

Re    = 6378100.  # Radius of earth (m)
Atm_h = 8300.     # vertical extent of air layer (m)
So    = 1367.     #solar constant (watts per square meter)
Earth_emissivity = 0.6
Sigma = 5.67e-8   # Boltzmann constant
Tzero = 273.15    # absolute zero
G     = 9.807
R     = 8.3143    # universal gas constant
Ma     = 0.0289644  # molar mass of dry air
Mw    = 0.01801528 # molar mass of water
Psea  = 101325     # sea level standard pressure
Tsea  = 288.15     # sea level standard temp

def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def calculate_pressure_heights(pressure, temp):
    return -math.log (pressure) * R * temp / (Ma * G)

def coriolis_force(lat):
    min_val = 0.4e-4
    val = -2 * 7.292e-5 * math.sin(math.radians(lat))
    if val >= min_val or val <= -min_val:
        return val
    else:
        return min_val if val > 0 else -min_val

def v_wind(lat, grad_x):
    return grad_x * G / (coriolis_force(lat) * 1.11e5 )

def u_wind(lat, grad_y):
    return -grad_y * G / (coriolis_force(lat) * 1.11e5 )


def evaporation_per_minute(T_air, evap_coeff):
    """
    Very primitive calculation of evaporation of water (kilogramms/m^3 per minute)
    with temp <= 0 Celsium no evaporation
    with temp > 0 increase relative humidity is linear relative to air temp and max with t=50C
    return evaporation of water from surface (kg/ (m^3 * minute) )
    """
    t = T_air - Tzero
    if t <= 0:
        return 0
    t = 50 if t > 50 else t
    T_coeff = t / 50.0
    Pconst = 1e-8 # kg
    return evap_coeff * T_coeff * Pconst

def partial_atmospheric_water_pressure(density, T):
    """
    return pressure in Pa
    """
    return density * R * T / Mw

def equilibrium_water_pressure(P, T):
    """
    Using approximation Buck formula
    P - pressure in Pa
    T - temp
    """
    Pb = P / 100 #convert to millibars
    t = T - Tzero
    tcoeff = 17.502*t / (240.97 + t)
    Ewp =  1.0007 + 3.46e-6 * Pb * 6.1121 * pow (math.e, tcoeff) #equilibrium water pressure
    return  Ewp

def relative_humidity(density, T, P):
    """
    density - density of water vapor in atmosphere
    T - temp in K
    P - surface pressure Pa
    """
    return partial_atmospheric_water_pressure(density, T) / equilibrium_water_pressure(P, T)

def dew_point(T, RH):
    """
    Approximation formula to calculate dew point
    return temperature of dew point
    """
    return T - (100 - RH * 100) / 5


def altitude_pressure(altitude):
    return Psea * math.exp(-G * Ma * altitude / (R * Tsea) )

class SimulationGrid:
    def __init__(self, lats, lons, heights, land_sea_mask):
        self.lats = np.asarray(lats)
        self.reversed_lats = list(reversed(lats))
        self.lons = np.asarray(lons)
        self.U_wind = np.zeros(shape=(len(lats),len(lons)))
        self.V_wind = np.zeros(shape=(len(lats), len(lons)))
        self.air_temp = np.zeros(shape=(len(lats), len(lons)))
        self.air_temp.fill(14 + Tzero)
        self.mb500_height = np.zeros(shape=(len(lats), len(lons)))
        self.water_density = np.zeros(shape=(len(lats), len(lons)))
        self.humidity = np.zeros(shape=(len(lats), len(lons)))
        self.cloud_cover = np.zeros(shape=(len(lats), len(lons)))
        self.data = []
        for j,lat in enumerate(lats):
            sl = []
            for i,lon in enumerate(lons):
                lsm = land_sea_mask[j, i]
                sl.append ({
                    'height' : heights[j, i]
                    , 'evaporation-coeff' : 0.3 if lsm == 1 else 1.0
                })
            self.data.append(sl)

    def len_x(self):
        return len(self.lons)

    def len_y(self):
        return len(self.lats)

    def cell_size(self):
        """ length of cell side in meters """
        return 40075.16e3 /  self.len_x()

    def to_xy(self, lat, lon):
        x = bisect.bisect_right(self.lons, lon)
        y = bisect.bisect_right(self.reversed_lats, lat)
        y = len(self.reversed_lats) - 1 - y
        return np.clip(x-1, 0, self.len_x()) , np.clip(y-1, 0, self.len_y() - 1)


    def temperature_diff_per_minute(self, x, y, datetime, print_val=False):
        cell = self.data[y][x]
        temp = self.air_temp[y][x]
        lat = self.lats[y]
        lon = self.lons[x]
        clouds = self.cloud_cover[y][x]
        S_raw = solar_energy_influx(lat, lon, datetime)
        albedo = 0.3 + (0.6 * clouds)      #reflection coeff, increases with cloud
        S_absorbed = (1 - albedo) * S_raw
        E_emitted = 4.0 * Earth_emissivity * Sigma * math.pow(temp,4)  / 4.79
        val =  (S_absorbed ) / (4 * 8.3 * 1.2 * 1000)
        if print_val:
            eprint("T:", S_absorbed)
        return 2*val

    def get_cell_antipode(self, x, y):
        median_x = self.len_x() // 2
        #assert(median * 2 == self.len_x() , 'This will work just for even longitude size' )
        xa = median_x + x if x < median_x else x - median_x
        median_y = self.len_y() // 2
        xy_ = median_y - y  if median_y > y else y - median_y
        xy = median_y + xy_ if median_y > y else median_y - xy_
        return self.data[xa][xy-1]


    def get_source_cell_x(self, x, dir_x):
        if dir_x > 0:
            return x - 1 if x > 0 else  self.len_x() - 1 # left
        elif dir_x == 0:
            return x
        else: # dir_x < 0
            return x + 1 if x < self.len_x() - 1 else 0  #right

    def get_source_cell_y(self, y, dir_y):
        if dir_y > 0:
            return y - 1 if y > 0 else y #up
        elif dir_y == 0:
            return y
        else:
            return y + 1 if y < self.len_y() - 1 else y  #down

    def get_surface_pressure(self, x, y):
        h = self.data[y][x]['height']
        return altitude_pressure(h)

    def get_relative_humidity(self, x, y):
        Ps = self.get_surface_pressure(x, y)
        return relative_humidity(self.water_density[y][x], self.air_temp[y][x], Ps)

def update_plot(grid):
    eprint("pllot update")
    o = {}
    o['frame'] = 1
    o['temp'] = grid.air_temp.flatten().tolist()
    o['uwind'] = grid.U_wind.flatten().tolist()
    o['vwind'] = grid.V_wind.flatten().tolist()
    json.dump(o, sys.stdout)
    sys.stdout.write('\n')

def plor_print_header(lats, lons):
    o = {}
    o['title'] = "Temp and winds plot"
    o['lats'] = lats.tolist()
    o['lons'] = lons.tolist()
    json.dump(o, sys.stdout)
    sys.stdout.write('\n')

def get_or_download(e, force_reload=False):
    """ Get or download entry from datasets  """
    if(force_reload or not os.path.exists( e['file-name'] )):
        logging.info('Downloading %s from: %s', e['desc'], e['url'])
        with urllib.request.urlopen(e['url']) as response, open(e['file-name'], 'wb') as out_file:
            logging.info('Size: %s bytes', response.info()['Content-Length'])
            shutil.copyfileobj(response, out_file)
    else:
        logging.debug('Using cached version of \"' + e['file-name'] + '\"')
    return e['file-name']


def declination_of_sun_angle(day_of_year):
    """ Declination of sun angle in radians """
    N = day_of_year
    a = 0.98565 * (N + 10) + 1.914 * math.sin( math.radians( 0.98565 * (N - 2) ) )
    return -math.asin(0.39779 * math.cos( math.radians( a ) ) )

def solar_zenith_angle(lat, day_of_year, hour_angle):
    lat = math.radians(lat)
    hour_angle = math.radians(hour_angle)
    da = declination_of_sun_angle(day_of_year)
    a = math.sin(lat) * math.sin(da) + math.cos(lat) * math.cos(da) * math.cos (hour_angle)
    return math.acos(a)


def solar_energy_influx(lat, lon, dt):
    midnight = dt.replace(hour=0, minute=0, second=0, microsecond=0)
    year_begin = midnight.replace(month=1, day=1)
    day_of_year = (dt - year_begin).days
    minutes_from_midnight = (dt - midnight).seconds / 60
    """  Calculate solar energy influx at given coordinate and datetime """
    minutes_to_degrees = -minutes_from_midnight * 360 / (24 * 60) + 180
    hour_angle = lon - minutes_to_degrees
    sa = solar_zenith_angle(lat, day_of_year, hour_angle)
    cos_sa = math.cos(sa)
    return So * math.cos(sa) if cos_sa > 0 else 0


def advection_diff(u, v, T0, Txs, Tys, L, dt):
    """
    calculate parameter increase due to advection
    u   - wind speed by x
    v   - wind speed by y
    T0  - initial cell quantity
    Txs - x-source cell quantity
    Tys - y-source cell quantity
    L   - length of cell side in meters
    dt  - time in seconds
    """
    if u == 0 or v == 0:
        return 0,0,0
    dx = abs(u) * dt
    dy = abs(v) * dt
    S = L * L
    A_0 = (L - dy) * (L - dx) / S
    T_0 = A_0 * T0
    A_y = (L - dx) * dy / S
    T_y = A_y * Tys
    A_x = (L - dy) * dx / S
    T_x = A_x * Txs
    cx = dx / (dx + dy)
    cy = 1 - cx
    A_xy = 1 - (A_0 + A_y + A_x)
    T_x += Txs * cx * A_xy
#    Tsxy = Txs * cx + Tys * cy
#    T_xy = A_xy * Tsxy
    T_y += Tys * cy * A_xy
    T1 = T_0 + T_y + T_x # + T_xy
    if abs(T1 - T0) > 10 :
        S = L*L
        eprint("X:", (L - dx) * dy/ S, Txs)
        eprint("Y:", (L - dy) * dx/ S, Tys)
        eprint("0:", (L - dy) * (L - dx)/ S, T0)
        eprint("XY:", dx * dt/ S, Tsxy)
        eprint(T1 - T0, u, v, "T", T0-Tzero, Txs-Tzero, Tys-Tzero, Tsxy-Tzero, "C", A_0, A_x, A_y, A_xy)
    return T1 - T0, T_x, T_y


def integrate_step(grid, cur_time, step_minutes):
    temp_diff = np.zeros(shape=(grid.len_y(), grid.len_x()))
    u_wind_diff = np.zeros(shape=(grid.len_y(), grid.len_x()))
    v_wind_diff = np.zeros(shape=(grid.len_y(), grid.len_x()))
    water_density_diff = np.zeros(shape=(grid.len_y(), grid.len_x()))
    for x in range(grid.len_x()):
        for y in range(grid.len_y()):
            temp_diff[y][x] += grid.temperature_diff_per_minute(x,y,cur_time)  * step_minutes
            x_src = grid.get_source_cell_x(x, grid.U_wind[y][x])
            y_src = grid.get_source_cell_y(y, grid.V_wind[y][x])
            dt,dtx,dty = advection_diff(grid.U_wind[y][x],
                                        grid.V_wind[y][x],
                                        grid.air_temp[y][x],
                                        grid.air_temp[y][x_src],
                                        grid.air_temp[y_src][x],
                                        grid.cell_size(),
                                        step_minutes * 60)
            temp_diff[y][x] += dt

            evap_cf = grid.data[y][x]['evaporation-coeff']
            evap = evaporation_per_minute(grid.air_temp[y][x], evap_cf) * step_minutes
            water_density_diff[y][x] = evap

            dw,dwx,dwy = advection_diff(grid.U_wind[y][x],
                                        grid.V_wind[y][x],
                                        grid.water_density[y][x],
                                        grid.water_density[y][x_src],
                                        grid.water_density[y_src][x],
                                        grid.cell_size(),
                                        step_minutes * 60)

            water_density_diff[y][x] += dw
            water_density_diff[y][x_src] -= dwx
            water_density_diff[y_src][x] -= dwy

            grid.mb500_height[y][x] = calculate_pressure_heights(0.5, grid.air_temp[y][x])
            #cell_ant = grid.get_cell_antipode(x,y)
            #    cell_ant['air-temperature'] -= temp_diff
    total_absorbed = temp_diff.sum()
    temp_a = grid.air_temp * grid.air_temp
    temp_a = temp_a * temp_a * temp_a * temp_a
    sum_b = temp_a.sum()
    temp_a /= sum_b
    emitted_per_cell = temp_a * total_absorbed

    eprint(total_absorbed, sum_b)
    temp_diff -= emitted_per_cell

    grid.air_temp+=temp_diff
    grid.water_density += water_density_diff
    grad = np.gradient(grid.mb500_height)
    for x in range(grid.len_x()):
        for y in range(grid.len_y()):
            lat = grid.lats[y]
            v_wind_diff[y][x] = v_wind(lat,grad[1][y][x]/1.85)
            u_wind_diff[y][x] = u_wind(lat,grad[0][y][x]/1.85)
            grid.humidity[y][x] = grid.get_relative_humidity(x,y)
            Tcell = grid.air_temp[y][x]
            dew_pt_T = dew_point(Tcell, grid.humidity[y][x])
            if Tcell < dew_pt_T:
                grid.water_density[y][x] /= 2 # form clouds
                grid.humidity[y][x] = grid.get_relative_humidity(x,y)
                grid.cloud_cover[y][x] += 0.1


    grid.U_wind = u_wind_diff
    grid.V_wind = v_wind_diff
    eprint("G0",grad[0].min(),grad[0].max())
    eprint("G1",grad[1].min(),grad[1].max())
    eprint("U",grid.U_wind.min(),grid.U_wind.max())
    eprint("V",grid.V_wind.min(),grid.V_wind.max())
    eprint("T",grid.air_temp.min() - Tzero ,grid.air_temp.max() - Tzero)
    eprint("W",grid.water_density.min(),grid.water_density.max())
    eprint("RH",grid.humidity.min(),grid.humidity.max())
    eprint("Cl",grid.cloud_cover.min(),grid.cloud_cover.max())


def integrate(grid, step_minutes, num_days, start_date, plot_data=False):
    delta = timedelta(minutes=step_minutes)
    end_date = start_date + timedelta(days=num_days)
    cur_time = start_date
    x,y = grid.to_xy(52.53,13.35)


    j = 0
    if plot_data:
        plor_print_header(grid.lats, grid.lons)
    while cur_time < end_date:
        if plot_data and j % 2 == 1:
            update_plot(grid )
        j+=1
        cur_time = cur_time + delta
        #x,y = grid.to_xy(0,0)
        integrate_step(grid, cur_time, step_minutes)
        #cell = grid.data[x][y]
        eprint ("{}:Temperature in Berlin: {}, absorbed solar EL {}, RH: {}, {}, {}".format(
            cur_time,
            grid.air_temp[y][x] - Tzero,
            grid.temperature_diff_per_minute(x,y,cur_time, True),
            grid.get_relative_humidity(x, y),
            equilibrium_water_pressure(Psea, grid.air_temp[y][x]  ),
            partial_atmospheric_water_pressure(grid.water_density[y][x], grid.air_temp[x][y])
        ))
        i = 0
        eprint("Ttoal:", np.average(grid.air_temp) - Tzero )

def run(args):
    logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s',
                        level=logging.DEBUG if args.debug else logging.INFO )
    logging.info('Starting program..')
    for key in datasets:
        get_or_download(datasets[key], args.force_reload)
    with Dataset(datasets['surface-height']['file-name'],'r') as hgt_file,\
         Dataset(datasets['land-sea-mask']['file-name'],'r') as land_file:
        lat_list = list(hgt_file.variables['lat'])
        lon_list = list(hgt_file.variables['lon'])
        hgt = hgt_file.variables['hgt'][0,:,:]
        land = land_file.variables['land'][0,:,:]
        grid = SimulationGrid(lat_list, lon_list, hgt, land)
        integrate(grid, 30, 30, datetime(2012,4,1), args.plot)

def parse_args(argv):
    parser = argparse.ArgumentParser('Import initial data for dewpoint weather simulator')
    parser.add_argument("-d","--debug", help="print debug statements", action='store_true')
    parser.add_argument("-p","--plot", help="output plot data to stdout", action='store_true')
    parser.add_argument("--force-reload",
                        help="force reload, ignore cached files",
                        action='store_true')
    return parser.parse_args(argv)

def print_gradeint():
    ar = np.array([[0,0,0],[1,1,1],[0,0,0]])
    gr = np.gradient(ar)
    eprint(ar)
    eprint(gr[0])
    eprint(gr[1])

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = parse_args(argv[1:])
    if (not args):
        return 2
    run(args)
#    print_gradeint()

if __name__ == "__main__":
    sys.exit(main())
