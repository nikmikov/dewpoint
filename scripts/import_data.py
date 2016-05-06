#!/usr/bin/env python

import sys
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

datasets = {
    'geopotentail-height': {
        'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/hgt.sfc.gauss.nc'
        , 'file-name': 'hgt.sfc.gauss.nc'
        , 'desc': 'Geopotential height data'
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

def get_or_download(e, force_reload=False):
    if(force_reload or not os.path.exists( e['file-name'] )):
        logging.info('Downloading %s from: %s', e['desc'], e['url'])
        with urllib.request.urlopen(e['url']) as response, open(e['file-name'], 'wb') as out_file:
            logging.info('Size: %s bytes', response.info()['Content-Length'])
            shutil.copyfileobj(response, out_file)
    else:
        logging.debug('Using cached version of \"' + e['file-name'] + '\"')
    return e['file-name']


def plot(bmaps, Z):
    plt.figure(figsize=(16,16))
    plt.contourf(bmaps['global_x'],bmaps['global_y'], Z , np.linspace(0, 1500,200), extend='both',antialiasing=False)
    bmaps['global'].drawcoastlines()
    plt.colorbar()
    plt.show()
    print('DONE')


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

def run(args):
    logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s',
                        level=logging.DEBUG if args.debug else logging.INFO )
    logging.info('Starting program..')
    for key in datasets:
        get_or_download(datasets[key], args.force_reload)
    with Dataset(datasets['geopotentail-height']['file-name'],'r') as infile:
        lat_list = list(infile.variables['lat'][1:-1])
        lon_list = list(infile.variables['lon'][:])
#        lon_list.append(360.)
        bm =  create_basemaps(lat_list, lon_list)
        hgt = infile.variables['hgt'][0,1:-1,:]
        hgt_last = hgt[0,:]
#        hgt_tran = np.transpose(hgt)
        new_hgt = np.vstack((hgt,hgt_last))
        # Duplicate the last longitude column for periodicity
        plot(bm, hgt )


def parse_args(argv):
    parser = argparse.ArgumentParser('Import initial data for dewpoint weather simulator')
    parser.add_argument("-d","--debug", help="print debug statements", action='store_true')
    parser.add_argument("--force-reload",
                        help="force reload, ignore cached files",
                        action='store_true')
    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = parse_args(argv[1:])
    if (not args):
        return 2
    run(args)

if __name__ == "__main__":
    sys.exit(main())
