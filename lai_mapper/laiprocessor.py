import gc
import os
import time
from glob import glob
import rioxarray
from rioxarray.merge import merge_arrays
import numpy as np
import rasterio as rio
import xarray as xr
import pandas as pd
from joblib import Parallel, delayed
from nasa_utils import (SessionWithHeaderRedirection, create_netrc,
                        get_modis_files_to_download, get_nasa_data,
                        get_nsidc_files_to_download)
from rasterio.transform import from_origin
from rasterio.warp import Resampling
from utility import merge_tile_data, warp_to_wgs84, write_array_tif

path = '.'

def driver(year, doy, user, passwd, sat='MOLA'):
    
    region=None
    
    if sat is 'MOLA' or sat is 'MOTA':
        s = 'MOTA'
        p1 = 'MOD15A2H.006'
        ext='.hdf'
    else:
        s = 'VIIRS'
        p1 = 'VNP15A2H.001'
        ext='.h5'
        
    #ESTABLISH NASA PASSWORD AND USER
    session = SessionWithHeaderRedirection(user,passwd)
    #make sure user has .netrc file in home directory
    create_netrc(user,passwd)


    target_crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    mcd43a1_path = os.path.join(os.getcwd(), '{}_{}_{}{}.tif'.format(s, p1, year, str(doy).zfill(3)))
    if not os.path.exists(mcd43a1_path.replace('.tif','.nc')):
        output_path = './'
        print('###################################################')
        print('Downloading {} {} Tiles for: {}'.format(s,p1,doy))
        print("###################################################")
        for i in range(4):
            fns, basenames_a1 = get_modis_files_to_download(year, doy, None, output_path, ext,sat=s,product=p1)
            if len(fns) > 0:
                Parallel(n_jobs=5,verbose=10)(delayed(download)(user,passwd,fn) for fn in fns)
            else:
                print('{} {} is all downloaded ......'.format(s,p1))
                break
            
            
        #check for tiff tiles already created
        print('################################################')
        print('          Writing {} to TIFs'.format(p1))
        print('################################################')
        tiffnames1 = [fname.replace(ext,'.tif') for fname in basenames_a1]
#        print(tiffnames1,basenames_a1)
        if s == 'VIIRS':
            print(basenames_a1)
            a1_reprojected = Parallel(n_jobs=24,verbose=10)(delayed(_convert_hdf_tif_vnp15A2H)(hdf,tif) for hdf,tif in zip(basenames_a1,tiffnames1))
        else:
            a1_reprojected = Parallel(n_jobs=24,verbose=10)(delayed(_convert_hdf_tif_mod15A2H)(hdf,tif) for hdf,tif in zip(basenames_a1,tiffnames1))
        gc.collect()
        #        print('################################################')
        #        print('          WARPING {} TIF Files'.format(p1))
        #        print('################################################')
        #        tiffwarp1 = [fname.replace('.tif','_warped.tif') for fname in tiffnames1]
        #        Parallel(n_jobs=12,verbose=10)(delayed(warp_to_wgs84)(fname) for fname in tiffnames1)
        
        #        print('################################################')
        #        print('        Merging {} Warped TIFS to File'.format(p1))
        #        print('################################################')
        #        pp = p1.split('.')[0] + '.A'
        #        tiffm1 = [fname.replace('.tif','_warped.tif') for fname in tiffnames1]
        #        tiffm1 = glob('{}{}{}.*_warped.tif'.format(pp,year, str(doy).zfill(3)))
        #        merge_tile_data(tiffm1, mcd43a1_path)
        #        gc.collect()

        # merge all to lat lon grid
        lai_tiles = [read_hdf_vnp15ah2(fname) for fname in tiffnames1] 
        lai_merged = merge_arrays(lai_tiles)
#        lai_merged.to_netcdf('lai_merged.nc')
        #lai_merged.name = 'LAI'
        lai = lai_merged.rio.reproject(target_crs,shape=(5000,10000),resampling=Resampling.bilinear)
        print('################################################')
        print('        Cleaning Up Files....')
        print('################################################')
        #files = glob('{}.A{}{}.*.tif'.format(year, str(doy).zfill(3)))
        Parallel(n_jobs=12,verbose=10)(delayed(os.remove)(fname) for fname in basenames_a1)
        Parallel(n_jobs=12,verbose=10)(delayed(os.remove)(fname) for fname in tiffnames1)
        #        files =	glob('{}.A{}{}.*.tif'.format(p2.plyear, str(doy).zfill(3)))
#        Parallel(n_jobs=12,verbose=10)(delayed(os.remove)(fname) for fname in tiffwarp1)
        
        # Print Cleaning up HDF files
        print('################################################')
        print('          Cleaning Up RAW HDF Files')
        print('################################################')
        #        Parallel(n_jobs=10,verbose=10)(delayed(os.remove)(fname) for fname in basenames_a1)
        #        Parallel(n_jobs=10,verbose=10)(delayed(os.remove)(fname) for fname in basenames_a2)
        
        # write to netcdf
        #lai = xr.open_rasterio(mcd43a1_path)
        lai.attrs['long-name'] = 'Leaf Area Index'
        #lai.attrs['_FillValue'] = 0.
        #lai.attrs['valid_range'] = [0,100]
        #lai.attrs['scale_factor'] = 1
        print('###############################################')
        print( '                      writting file')
        print('###############################################')
        lai = lai.where(lai > 0,0)
        das = dict(lai=lai)
        ds = xr.Dataset(das).rename(dict(x='longitude',y='latitude',band='level'))
        ds['longitude'].attrs['long_name'] = 'longitude'
        ds['longitude'].attrs['standard_name'] = 'longitude'
        ds['longitude'].attrs['units'] = 'degrees_east'
        ds['longitude'].attrs['axis'] = 'X'
        ds['latitude'].attrs['long_name'] = 'latitude'
        ds['latitude'].attrs['standard_name'] = 'latitude'
        ds['latitude'].attrs['units'] = 'degrees_north'
        ds['latitude'].attrs['axis'] = 'Y'

        d = pd.to_datetime('{}{}'.format(year,doy),format='%Y%j')
        hours = (d - pd.to_datetime('{}'.format(year), format='%Y'))/ np.timedelta64(1,'h')
        ds['time'] = hours
        ds.time.attrs['standard_name'] = 'time'
        ds.time.attrs['long_time'] = 'Time'
        ds.time.attrs['units'] = pd.to_datetime('{}-01-01'.format(year),format='%Y-%m-%d').strftime('hours since %Y-%m-%d %H:%M:%S')
        ds.time.attrs['calendar'] = 'standard'
        ds = ds.expand_dims('time').set_coords('time')
        yyyymmdd = d.strftime('%Y%m%d')
        newname = '{}_{}_{}.nc'.format(s,p1,yyyymmdd)
        write_ncf(ds,newname)
        return mcd43a1_path.replace('.tif','.nc')
        
                  
def write_ncf(dset,output_name):
    print('Writing:', output_name)
    comp = dict(zlib=True,complevel=5)
    encoding = {}
    for i in dset.data_vars.keys():
        encoding[i] = comp
    dset.attrs['github_url'] = 'https://github.com/bbakernoaa/laimapper'
    dset.attrs['contact'] = 'barry.baker@noaa.gov'
    dset.attrs['date_created'] = pd.to_datetime('today').strftime('%Y-%m-%d')
    dset.to_netcdf(output_name,encoding=encoding)

def download(user,passwd,fn):
    print("==========downloading:{}==========".format(fn))
    for attempt in range(1):  # keep retrying
        try:
            a = get_nasa_data('barrybaker', 'Tallguy1', fn)
            if a is None:
                print("attempt:{}".format(attempt))
                break
            else:
                if attempt == 10:
                    print("too many attempts! for {}".format(fn))
                    break
        except:
            time.sleep(90)
            print("attempt:{}".format(attempt+1))


def _convert_hdf_tiff(basename):
    # from utility import
    hdf_file_name = basename.split(os.sep)[-1]
    tif_file_name = hdf_file_name[:-3] + "tif"
    # print(fn)
    convert_hdf_tif_modis(hdf_file_name, tif_file_name)

def _modis_hv_extent(fname):
    tile = fname.split(".")[2]
    h = int(tile[1:3])
    v = int(tile[4:])
    ulx = -20015109.354 + h * 1111950.5196666666
    uly = -10007554.677 + (18 - v) * 1111950.5196666666
    res = 463.312716527917

def calc_omega_ns(mcd43a3,mcd43a1):
    wn = mcd43a3 / (mcd43a1 + 1e-3)
    wn = wn.where(xr.ufuncs.isfinite(wn))
    a = 0
    b = 1
    wnmax = 35.
    wnmin = 0.
    wns = (a-b)*(wn - wnmax) / (wnmin - wnmax) + b
    return wns

def _convert_hdf_tiff(basename):
    # from utility import
    hdf_file_name = basename.split(os.sep)[-1]
    tif_file_name = hdf_file_name[:-3] + "tif"
    # print(fn)
    convert_hdf_tif_modis(hdf_file_name, tif_file_name)

def _modis_hv_extent(fname):
    tile = fname.split(".")[2]
    h = int(tile[1:3])
    v = int(tile[4:])
    ulx = -20015109.354 + h * 1111950.5196666666
    uly = -10007554.677 + (18 - v) * 1111950.5196666666
    res = 463.312716527917

def _convert_hdf_tif_vnp15A2H(hdf_file_name,tif_file_name):
    crs = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
    tile = hdf_file_name.split(".")[2]
    h = int(tile[1:3])
    v = int(tile[4:])
    ulx = -20015109.354 + h * 1111950.5196666666
    uly = -10007554.677 + (18 - v) * 1111950.5196666666
    res = 463.312716527917
    transform = from_origin(ulx - res / 2, uly - res / 2, res, res)
    s = rio.open(hdf_file_name)
    if len(s.subdatasets) < 1:
        iso = s.read().squeeze()
    else:
        #field = 'HDF5:{}://HDFEOS/GRIDS/VNP_Grid_VNP15A2H/Data_Fields/Lai'.format(hdf_file_name)
        field = 'HDF5:{}://HDFEOS/GRIDS/VIIRS_Grid_LAIFPAR/Data_Fields/Lai'.format(hdf_file_name)
        band = [i for i in s.subdatasets if field in i][0]
#        print(rio.open(band).read().shape)
        iso = rio.open(band).read().squeeze()
        #iso = xr.open_rasterio(band)[0,:,:]
    #iso.data = (iso.where((iso != iso.nodatavals) & (iso < 249)) * 0.1 + iso.offsets).data
    #iso = iso.expand_dims('band')
    #return iso
    write_array_tif(iso, crs, transform, tif_file_name)

def read_hdf_vnp15ah2(fname):
    crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    #s = rio.open(fname)
    #band = [i for i in s.subdatasets if "/Lai" in i][0]
    iso = xr.open_rasterio(fname)
    iso.data = (iso.where(( iso < 249)) * 0.1).data
    iso.name = 'lai'
    return iso
