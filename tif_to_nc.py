#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 2023

@author: Todd Clark
"""

import xarray as xr
import rioxarray
import xesmf as xe
import numpy as np
from osgeo import gdal
import rasterio
from rasterio.vrt import WarpedVRT
import os

def rescale_resolution(input_file, output_file, scaling_factor):
    '''

    Parameters
    ----------
    input_file : string
        name of the file you want to rescale
    output_file : string
        name of file and location you want to rescale the file too
    scaling_factor : float
        example : 0.5 will make tiffs length and width half as man pixels rounded down

    Returns
    -------
    None.

    '''
    #open input file
    ds = gdal.Open(input_file)

    # Calculate the new dimensions
    new_width = int(ds.RasterXSize * scaling_factor)
    new_height = int(ds.RasterYSize * scaling_factor)

    # set options
    options = gdal.WarpOptions(
        format="GTiff",
        outputType= gdal.GDT_Float32,
        width=new_width,
        height=new_height,
        resampleAlg=gdal.GRA_Average,  # You can change the resampling method as needed
    )
    
    #transfor,
    gdal.Warp(output_file, ds, options=options)

    ds = None  # Close the input file
        
def create_mask(input_file, output_file, threshold):
    '''
    creates mask with ones representing values over a threshold and zeros for values under and equal to the threshold
    Parameters
    ----------
    input_file : string
        location of file you want to make mask of
    output_file : string
        location that you want to save the mask to
    threshold : float
        minimum value to be set at one

    Returns
    -------
    None.

    '''
    ds = gdal.Open(input_file)

    # Read the data
    data = ds.ReadAsArray()

    # Create a binary mask based on the threshold
    mask = (data > threshold).astype(np.uint8)

    # Save the mask as a GeoTIFF file
    mask_driver = gdal.GetDriverByName('GTiff')
    mask_ds = mask_driver.Create(output_file, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte)
    mask_ds.GetRasterBand(1).WriteArray(mask)

    # Set the geo referencing information for the mask
    mask_ds.SetProjection(ds.GetProjection())
    mask_ds.SetGeoTransform(ds.GetGeoTransform())

    ds = None  # Close the input file
    mask_ds = None  # Close the mask file

def transform_5070_to_4326(tif_path_5070, tif_path_4326):
    '''
    takes tiff input in 5070 and transforms it to 4326

    Parameters
    ----------
    tif_path_5070 : string
        path to input 5070 tif file.
    tif_path_4326 : string
        path to ouput 4326 tif file

    Returns
    -------
    None.

    '''
    
    
    input_ds = gdal.Open(tif_path_5070)
    driver = gdal.GetDriverByName("GTiff")
    output_ds = gdal.Warp(tif_path_4326, input_ds, dstSRS="EPSG:4326")

def convert_to_geochem_nc(tif_file_path, nc_output_path, template_ds):
    '''
    

    Parameters
    ----------
    tif_file_path : TYPE
        DESCRIPTION.
    nc_output_path : TYPE
        DESCRIPTION.
    template_ds : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    #open file as xarray object using rioxarray
    input_ds = rioxarray.open_rasterio(tif_file_path)
    
    #convert to data set
    input_ds = input_ds.to_dataset(name = "input")
    input_ds.rio.write_grid_mapping(inplace = True) # I don't know what this line does
    
    #rename axes so data will work with regridder
    input_ds = input_ds.rename({'x':'lon','y':'lat'})
    
    #create regrid object
    regridder = xe.Regridder(input_ds, template_ds, "conservative")
    
    #regrid xarray object
    template_ds = regridder(input_ds)
    template_ds.to_netcdf(nc_output_path)
    
def tif_to_nc(fulres_4326_path, lowres_4326_path, f_plya_nc_path, tif_5070_path, template_path,scaling_factor = 0.01):
    
    #define template
    template_ds=xr.open_dataset(template_path)
    
    #run methods
    transform_5070_to_4326 (tif_5070_path, fulres_4326_path)
    rescale_resolution(fulres_4326_path, lowres_4326_path, scaling_factor)
    convert_to_geochem_nc(lowres_4326_path, f_plya_nc_path, template_ds)

# if __name__ == "__main__":
#     main()
