"""
Created on Sun Oct 26 2023

@author: Todd Clark
"""

import rioxarray
import xesmf as xe
import numpy as np
import numpy.matlib as npm
from osgeo import gdal
import rasterio

def replace_none_with_0(input_file_path, output_file_path):
    xr_array = rioxarray.open_rasterio(input_file_path)
    
    #replace possible null values
    xr_array = xr_array.where(xr_array != xr_array._FillValue, 0) #remove array null value
    xr_array = xr_array.where(xr_array.notnull(),0) #remove Nan, None, or NaT
    
    xr_array.rio.to_raster(output_file_path)

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
        outputType = gdal.GDT_Float32, #outputType= gdal.GDT_Int64, this needs to be changed based on what is going to be done with the data for salt I am using a float because the data may have been fractionalized when being warped.
        width=new_width,
        height=new_height,
        resampleAlg=gdal.GRA_Sum,  # You can change the resampling method as needed
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
    takes tiff input in 5070 and transforms it to 4326 (I don't know what error is created here)

    Parameters
    ----------
    tif_path_5070 : string
        path to input 5070 tif file.
    tif_path_4326 : string
        path to output 4326 tif file
    Returns
    -------
    None.
    '''
    input_ds = gdal.Open(tif_path_5070)
    
        # set options
    options = gdal.WarpOptions(
        format="GTiff",
        outputType = gdal.GDT_Float32,
        srcSRS = "EPSG:5070", 
        dstSRS ="EPSG:4326",
        resampleAlg = gdal.GRA_Sum
    )
    
    gdal.Warp(tif_path_4326, input_ds, options = options) # goes from srcSRS to dstSRS
    input_ds = None #close the data set

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
    # reggrider weights
    
    #rename axes so data will work with regridder
    input_ds = input_ds.rename({'x':'lon','y':'lat'})
    
    #create regrid object
    regridder = xe.Regridder(input_ds, template_ds, "conservative")
    ############################################################################################################################################
    # I AM NOT SURE IF THE xe.Regridder NORMALIZES FLUX VALUE TO CELL AREA
    # SEE REFERENCE https://journals.ametsoc.org/view/journals/mwre/127/9/1520-0493_1999_127_2204_fasocr_2.0.co_2.xml
    # I also don't know why I'm getting /uufs/chpc.utah.edu/common/home/u1285966/mambaforge-pypy3/envs/xesmf_env/lib/python3.11/site-packages/xesmf/backend.py:56: UserWarning: Latitude is outside of [-90, 90]
    # warnings.warn('Latitude is outside of [-90, 90]') because the data is in epsg 
    ############################################################################################################################################
    
    #regrid xarray object
    regrided_ds = regridder(input_ds) #when I plot this it looks like data disappears around the ocean
    #regrided_ds = regrided_ds.where(regrided_ds > 0, 0) # this replaces all negative values with zeros. (I am not sure where the negatives come form)
    
    # multiply by area to get what I should be total salt (could be error in )
    
    #regrided_ds = regrided_ds * template_ds["AREA"] #havent tested this # this doesn't work. Never use this
    
    regrided_ds.to_netcdf(nc_output_path)
    
    
def build_area_array(file_path):
    with rasterio.open(file_path) as testif:
        rast = testif.read(1)
        gt = testif.transform
        pix_width = gt[0]
        ulX = gt[2]
        ulY = gt[5]
        rows = testif.height
        cols = testif.width
        lrX = ulX + gt[0] * cols
        lrY = ulY + gt[4] * rows

    lats = np.linspace(ulY,lrY,rows+1)

    a = 6378137
    b = 6356752.3142

    # Degrees to radians
    lats = lats * np.pi/180

    # Intermediate vars
    e = np.sqrt(1-(b/a)**2)
    sinlats = np.sin(lats)
    zm = 1 - e * sinlats
    zp = 1 + e * sinlats

    # Distance between meridians
    #        q = np.diff(longs)/360
    q = pix_width/360

    # Compute areas for each latitude in square km
    areas_to_equator = np.pi * b**2 * ((2*np.arctanh(e*sinlats) / (2*e) + sinlats / (zp*zm))) / 10**6
    areas_between_lats = np.diff(areas_to_equator)
    areas_cells = np.abs(areas_between_lats) * q

    areagrid = np.transpose(npm.repmat(areas_cells,cols,1))
    return areagrid
    