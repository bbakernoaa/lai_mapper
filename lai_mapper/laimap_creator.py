#!/usr/bin/env python

import argparse
import gc
#from utility import warp_like
import os
from subprocess import run

import laiprocessor
#import make_map
import nasa_utils as nu
import pandas as pd
import xarray as xr


def normit(dset, vmin, vmax):
    return (dset - vmin) / (vmax - vmin)


def write_ncf(dset, output_name):
    comp = dict(zlib=True, complevel=5)
    encoding = {dset.name: comp}
    if 'time' not in dset.dims:
        final = dset.expand_dims('time')
    else:
        final = dset
    final.attrs['github_url'] = 'https://github.com/bbakernoaa/lai_mapper'
    final.attrs['contact'] = 'barry.baker@noaa.gov'
    final.attrs['date_created'] = pd.to_datetime('today').strftime('%Y-%m-%d')
    final.attrs['transform'] = (0.004, 0.0, -179.9999999971294, 0.0, -0.004,
                                69.99791666038062)
    final.attrs['crs'] = '+init=epsg:4326'
    final.attrs['res'] = (0.004, 0.004)
    final.to_netcdf(output_name, encoding=encoding)


def main_laiprocessor(dates, user, passwd, sat):
    #VIIRS_VNP15A2H.001_2019329.nc
    if sat == 'VIIRS':
        product = 'VNP15A2H.001'
    else:
        product = 'MCD15A2H.006'


#    print(dates)
#    print(sat,product)
    avail_dates = nu.get_modis_available_days(dates,sat=sat, product=product)
    print(avail_dates)
    for day in avail_dates:#dates:  #avail_dates:
        year = day.strftime('%Y')
        doy = day.strftime('%j').zfill(3)
        print(
            '#########################################################################'
        )
        print('CREATING Baker-Schepanski MAP FOR: ', sat, day)
        print(
            '#########################################################################'
        )
        try:
            yyyymmdd = day.strftime('%Y%m%d')
            testfile = '{}_{}_{}.nc'.format(sat,product,yyyymmdd)
            if os.path.isfile(testfile):
                print('file exists moving on')
            else:
                ssm_path = laiprocessor.driver(day.strftime('%Y'),
                                      day.strftime('%j'),
                                      user,
                                      passwd,
                                      sat=sat)
        except:
            print('proceeding....')

            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='download and gridding MODIS data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-s',
        '--start_doy',
        type=str,
        default=None,
        help='start day of processing. *Note: leave blank for Real-time format YYYY-mm-dd')
    parser.add_argument(
        '-e',
        '--end_doy',
        type=str,
        default=None,
        help='end day of processing. *Note: leave blank for Real-time: format YYYY-mm-dd')
    parser.add_argument('-u',
                        '--user',
                        default='barrybaker',
                        help='NASA USERNAME')
    parser.add_argument('-p',
                        '--password',
                        default='PASSWORD',
                        help='NASA PASSWORD')
    parser.add_argument(
        '--satellite',
        default='VIIRS',
        help='Satellite product: MOLA (AQUA), MOTA (TERRA), VIIRS')
    args = parser.parse_args()

    start_doy = args.start_doy
    end_doy = args.end_doy
    start = pd.to_datetime(start_doy)
    end = pd.to_datetime(end_doy)
#    start = pd.to_datetime('{}{}'.format(str(args.year),
#                                         str(start_doy).zfill(3)),
#                           format='%Y%j')
#    end = pd.to_datetime('{}{}'.format(str(args.year),
#                                       str(end_doy).zfill(3)),
#                         format='%Y%j')

    dates = pd.date_range(start=start, end=end)
    main_laiprocessor(dates=dates,
                     user=args.user,
                     passwd=args.password,
                     sat=args.satellite)
