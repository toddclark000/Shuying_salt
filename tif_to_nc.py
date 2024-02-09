"""
Created on Wed Jan 17 2024

@author: Todd Clark
"""
import xarray as xr
from helper_functions import *
    
def tif_to_nc(fulrez_4326_path, lowrez_4326_path, output_nc_path, tif_5070_path, tif_5070_path_none_with_0, fulrez_4326_path_none_with_0, template_path,scaling_factor = 0.01):
    """
    Converts a TIFF file to a NetCDF file using a given template.

    Args:
        fulres_4326_path (str): The path to where the intermediate full-resolution TIFF will be stored in EPSG:4326 projection.
        lowres_4326_path (str): The path to where the intermediate low-resolution TIFF will be stored file in EPSG:4326 projection.
        output_nc_path (str): The path to save the output NetCDF file.
        tif_5070_path (str): The path to the input TIFF file in EPSG:5070 projection is stored.
        template_path (str): The path to the template NetCDF file.
        scaling_factor (float, optional): The scaling factor to apply during the resolution rescaling. Defaults to 0.01

    Returns:
        None

    Examples:
        tif_to_nc('fulres.tif', 'lowres.tif', 'result.nc', 'tif_5070.tif', 'template.nc', scaling_factor=0.1)
    """

    #define template
    template_ds=xr.open_dataset(template_path)
    
    #run methods
    replace_none_with_0(tif_5070_path, tif_5070_path_none_with_0)
    transform_5070_to_4326(tif_5070_path_none_with_0, fulrez_4326_path)
    replace_none_with_0(fulrez_4326_path, fulrez_4326_path_none_with_0)
    rescale_resolution(fulrez_4326_path_none_with_0, lowrez_4326_path, scaling_factor)
    #convert_to_geochem_nc(lowrez_4326_path, output_nc_path, template_ds)

def tif_to_nc_without_replacing_zeros(fulres_4326_path, lowres_4326_path, output_nc_path, tif_5070_path, template_path,scaling_factor = 0.01):
    """
    Converts a TIFF file to a NetCDF file using a given template.

    Args:
        fulres_4326_path (str): The path to where the intermediate full-resolution TIFF will be stored in EPSG:4326 projection.
        lowres_4326_path (str): The path to where the intermediate low-resolution TIFF will be stored file in EPSG:4326 projection.
        output_nc_path (str): The path to save the output NetCDF file.
        tif_5070_path (str): The path to the input TIFF file in EPSG:5070 projection is stored.
        template_path (str): The path to the template NetCDF file.
        scaling_factor (float, optional): The scaling factor to apply during the resolution rescaling. Defaults to 0.01

    Returns:
        None

    Examples:
        tif_to_nc('fulres.tif', 'lowres.tif', 'result.nc', 'tif_5070.tif', 'template.nc', scaling_factor=0.1)
    """

    #define template
    template_ds=xr.open_dataset(template_path)
    
    #run methods
    #replace_none_with_0(tif_5070_path, tif_5070_path_none_with_0)
    transform_5070_to_4326(tif_5070_path, fulres_4326_path)
    rescale_resolution(fulres_4326_path, lowres_4326_path, scaling_factor)
    convert_to_geochem_nc(lowres_4326_path, output_nc_path, template_ds)
