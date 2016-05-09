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
from mpl_toolkits.basemap import Basemap
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime,timedelta

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

Re    = 6378100.  # Radius of earth (m)
Atm_h = 8300.     # vertical extent of air layer (m)
So    = 1367.     #solar constant (watts per square meter)
Earth_emissivity = 0.6
Sigma = 5.67e-8   # Boltzmann constant
Tzero = 273.15    # absolute zero
G     = 9.807
R     = 8.3143    # universal gas constant
M     = 0.0289644 # molar mass of dry air

def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def calculate_pressure_heights(pressure, temp):
    return -math.log (pressure) * R * temp / (M * G)

def coriolis_force(lat):
    val = 2 * 7.292e-5 * math.sin(math.radians(lat))
    return val

def v_wind(lat, grad_x):
    return grad_x * G / (coriolis_force(lat) * 1.11e5 * 1.85)

def u_wind(lat, grad_y):
    return -grad_y * G / (coriolis_force(lat) * 1.11e5* 1.85)

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
        self.data = []
        for j,lon in enumerate(lons):
            sl = []
            for i,lat in enumerate(lats):
                lsm = land_sea_mask[i, j]
                sl.append ({
                    'height' : heights[i, j]
                    , 'humidity' : 0
                    , 'evaporation-coeff' : 0.5 if lsm == 1 else 1.0
                    , 'cloud-cover' : 0.        # 0..1
                })
            self.data.append(sl)

    def len_x(self):
        return len(self.lons)

    def len_y(self):
        return len(self.lats)

    def to_xy(self, lat, lon):
        x = bisect.bisect_right(self.lons, lon)
        y = bisect.bisect_right(self.reversed_lats, lat)
        y = len(self.reversed_lats) - 1 - y
        return np.clip(x-1, 0, self.len_x()) , np.clip(y-1, 0, self.len_y() - 1)


    def temperature_diff_per_minute(self, x, y, datetime):
        cell = self.data[x][y]
        temp = self.air_temp[y][x]
        lat = self.lats[y]
        lon = self.lons[x]
        S_raw = solar_energy_influx(lat, lon, datetime)
        albedo = 0.3 + (0.6 * cell['cloud-cover'])      #reflection coeff, increases with cloud
        S_absorbed = (1 - albedo) * S_raw
        E_emitted = 4.0 * Earth_emissivity * Sigma * math.pow(temp,4)  / 4.45
        val =  (S_absorbed - E_emitted) / (4 * 8.3 * 1.2 * 1000)
        return  2*val

    def get_cell_antipode(self, x, y):
        median_x = self.len_x() // 2
        #assert(median * 2 == self.len_x() , 'This will work just for even longitude size' )
        xa = median_x + x if x < median_x else x - median_x
        median_y = self.len_y() // 2
        xy_ = median_y - y  if median_y > y else y - median_y
        xy = median_y + xy_ if median_y > y else median_y - xy_
        return self.data[xa][xy-1]


def plot_heights(bmaps, grid):
    Z = grid.air_temp
    U = grid.U_wind
    V = grid.V_wind

    plt.contourf(bmaps['global_x'],bmaps['global_y'], Z,\
                 extend='both',antialiasing=False)
    plt.hold(True)
    plt.quiver(bmaps['global_x'],bmaps['global_y'],U, V)
    bmaps['global'].drawcoastlines()
#    plt.colorbar()
    plt.show()


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
    minutes_to_degrees = minutes_from_midnight * 360 / (24 * 60) + 180.0
    hour_angle = minutes_to_degrees - lon
    sa = solar_zenith_angle(lat, day_of_year, hour_angle)
    cos_sa = math.cos(sa)
    return So * math.cos(sa) if cos_sa >0 else 0


def create_basemaps(lats, lons):
    """ Setup global basemaps for eventual plotting """
    print("Creating basemaps for plotting")

    long, latg = np.meshgrid(lons,lats)

    # Set up a global map
    bmap_globe = Basemap(projection='cea',llcrnrlat=-70, urcrnrlat=70,\
                         llcrnrlon=0,urcrnrlon=360,lat_ts=20,resolution='c')
    xg,yg = bmap_globe(long,latg)

    return {'global' : bmap_globe,
            'global_x' : xg,
            'global_y' : yg,
            }

def integrate_step(grid, cur_time, step_minutes):
    for x in range(grid.len_x()):
        for y in range(grid.len_y()):
            cell = grid.data[x][y]
            temp_diff = grid.temperature_diff_per_minute(x,y,cur_time)  * step_minutes
#            if temp_diff > 0:
            grid.mb500_height[y][x] = calculate_pressure_heights(0.5, grid.air_temp[y][x])
            grid.air_temp[y][x] += temp_diff
            #cell_ant = grid.get_cell_antipode(x,y)
            #    cell_ant['air-temperature'] -= temp_diff

    grad = np.gradient(grid.mb500_height)
    print("0",grad[0].min(),grad[0].max())
    print("1",grad[1].min(),grad[1].max())
    for x in range(grid.len_x()):
        for y in range(grid.len_y()):
            cell = grid.data[x][y]
            lat = grid.lats[y]
            grid.V_wind[y][x] =  v_wind(lat,grad[1][y][x])
            grid.U_wind[y][x] =  u_wind(lat,grad[0][y][x])
    print("U",grid.U_wind.min(),grid.U_wind.max())
    print("V",grid.V_wind.min(),grid.V_wind.max())
    print("T",grid.air_temp.min(),grid.air_temp.max())



def integrate(bm, grid, step_minutes, num_days, start_date):
    delta = timedelta(minutes=step_minutes)
    end_date = start_date + timedelta(days=num_days)
    cur_time = start_date
    x,y = grid.to_xy(52.53,13.35)

    plt.figure() #figsize=(10,10))
    j = 0
    while cur_time < end_date:
        if j % 10 == 0:
            plot_heights(bm, grid )
        j+=1
        cur_time = cur_time + delta
        #x,y = grid.to_xy(0,0)
        integrate_step(grid, cur_time, step_minutes)
        #cell = grid.data[x][y]
        print ("{}:Temperature in Berlin: {}".format(cur_time,grid.air_temp[x][y] - Tzero) )
        i = 0
        print("Ttoal:", np.average(grid.air_temp) - Tzero )

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
        bm =  create_basemaps(lat_list, lon_list)
        hgt = hgt_file.variables['hgt'][0,:,:]
        land = land_file.variables['land'][0,:,:]
        grid = SimulationGrid(lat_list, lon_list, hgt, land)
        integrate(bm, grid, 10, 10, datetime(2012,1,1))
#        lon_list.append(360.)
#        bm =  create_basemaps(lat_list, lon_list)
#        hgt_tran = np.transpose(hgt)
        # Duplicate the last longitude column for periodicity
#        plot_heights(bm, grid )
        #time_ar = [solar_energy_influx(52.53, 13.38, x, 30*6) for x in range(0, 24 * 360) ]
        #print(sum(time_ar) / (24*60))


def parse_args(argv):
    parser = argparse.ArgumentParser('Import initial data for dewpoint weather simulator')
    parser.add_argument("-d","--debug", help="print debug statements", action='store_true')
    parser.add_argument("--force-reload",
                        help="force reload, ignore cached files",
                        action='store_true')
    return parser.parse_args(argv)

def print_gradeint():
    ar = np.array([[1,2,3],[4,5,6],[7,8,9]])
    gr = np.gradient(ar)
    print(ar)
    print(gr[0])
    print(gr[1])

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = parse_args(argv[1:])
    if (not args):
        return 2
    run(args)
    #print_gradeint()

if __name__ == "__main__":
    sys.exit(main())
