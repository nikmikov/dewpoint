#!/usr/bin/env python

import sys
import math
import bisect
import argparse
import urllib.request
import shutil
import os.path
from netCDF4 import Dataset
import numpy as np
import json
import csv

datasets = {
    'stations': {
        'url': 'http://www.weathergraphics.com/identifiers/master-location-identifier-database-20130801.csv'
        , 'file-name': 'stations-database-input.csv'
        , 'desc': 'Database of meteostations with coordinates'
    },
    'surface-height': {
        'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/hgt.sfc.gauss.nc'
        , 'file-name': 'hgt.sfc.gauss.nc'
        , 'desc': 'Geopotential height at the surface'
    },
    'land-sea-mask' : {
        'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/land.sfc.gauss.nc'
        , 'file-name': 'land.sfc.gauss.nc'
        , 'desc': 'land-sea mask'
    }
}

def get_or_download(e, force_reload=False):
    """ Get or download entry from datasets  """
    if(force_reload or not os.path.exists( e['file-name'] )):
        print('Downloading {} from: {}'.format( e['desc'], e['url']))
        with urllib.request.urlopen(e['url']) as response, open(e['file-name'], 'wb') as out_file:
            print('Size: {} bytes'.format( response.info()['Content-Length']))
            shutil.copyfileobj(response, out_file)
    else:
        print('Using cached version of \"' + e['file-name'] + '\"')
    return e['file-name']


def run(args):
    print('Starting program..')
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
    parser.add_argument("--force-reload", default=False,
                        help="force reload, ignore cached files",
                        action='store_true')
    parser.add_argument("output_dir", help="output directory", nargs=1)
    return parser.parse_args(argv)

def import_stations(data_dir, force_reload):
    stations_input = get_or_download(datasets['stations'], force_reload)
    d = os.path.abspath(data_dir)
    if not os.path.exists(d):
        os.makedirs(d)
    output_name = os.path.join(d, 'stations.csv')
    with open(stations_input, newline='', encoding='windows-1252') as csvinput,\
         open(output_name, 'w', newline='') as csvoutput:
        reader = csv.reader(csvinput, delimiter=',',skipinitialspace=True )
        fieldnames = ['iata', 'country', 'city', 'station_name', 'lat', 'lon']
        writer = csv.DictWriter(csvoutput, delimiter=';', fieldnames=fieldnames)
        uniq_set = set("")
        for row in reader:
            o = {
                'iata' : row[13]
                , 'country': row[2]
                , 'city' : row[5]
                , 'station_name' : row[6]
                , 'lat'  : row[29]
                , 'lon'  : row[30]
            }
            if o['iata'] not in uniq_set:
                uniq_set.add(o['iata'])
                writer.writerow(o)
        print("Station list imported!")

def import_grid(data_dir, force_reload):
    land_sea_mask = get_or_download(datasets['land-sea-mask'], force_reload)
    surface_height = get_or_download(datasets['surface-height'], force_reload)
    d = os.path.abspath(data_dir)
    if not os.path.exists(d):
        os.makedirs(d)
    output_name = os.path.join(d, 'grid.json')
    j = {}
    with Dataset(datasets['surface-height']['file-name'],'r') as hgt_file,\
         Dataset(datasets['land-sea-mask']['file-name'],'r') as land_file,\
         open(output_name, 'w') as jsonoutput:
        j = {}
        j['lats'] = np.asarray(hgt_file.variables['lat']).tolist()
        j['lons'] = np.asarray(hgt_file.variables['lon']).tolist()
        heights = np.asarray(hgt_file.variables['hgt'][0,:,:])
        j['height'] = np.flipud(heights).transpose().flatten().tolist()
        landseamask =  np.asarray(land_file.variables['land'][0,:,:])
        j['land_sea_mask'] = np.flipud(landseamask).transpose().flatten().tolist()
        json.dump(j, jsonoutput)
        print("Grid imported!")

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = parse_args(argv[1:])
    if (not args):
        return 2


    import_stations(args.output_dir[0], args.force_reload)
    import_grid(args.output_dir[0], args.force_reload)


if __name__ == "__main__":
    sys.exit(main())
