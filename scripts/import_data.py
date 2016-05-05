#!/usr/bin/env python

import sys
import argparse
import urllib.request
import shutil
import os.path
import logging

datasets = {
    'geopotentail-heigt': {
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


def run(args):
    logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s',
                        level=logging.DEBUG if args.debug else logging.INFO )
    logging.info('Starting program..')
    for key in datasets:
        get_or_download(datasets[key], args.force_reload)

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
