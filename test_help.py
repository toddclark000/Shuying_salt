#from helper_functions import *
from osgeo import gdal
gdal.DontUseExceptions() #if you say UseExceptions() instead it makes it so you have to use a lower scaling_factor
import numpy as np
import pandas as pd
import rioxarray
import matplotlib.pyplot as plt

#from within the project
from tif_to_nc import tif_to_nc
from tif_to_nc import tif_to_nc_without_replacing_zeros

#TODO: I need to fix the sig-figs on range print outs

# create slice files
def slice_raster (input_file_path, output_folder, plot = False):
    
    #open array snd convert none values to zero
    xr_array = rioxarray.open_rasterio(input_file_path)
    xr_array = xr_array.where(xr_array != xr_array._FillValue,0)
    
    slices = np.arange(0, int(3.5e6)+1, int(0.5e6))
    
    slice_paths = []
    for index in range(len(slices)-1):
        
        slice_output_path = f"{output_folder}raster_fulrez_5070_{slices[index]:.0e}_to_{slices[index+1]:.0e}.tif"
        slice_paths.append(slice_output_path)
        
        xr_slice = xr_array.where(xr_array.y>slices[index]).where(xr_array.y<slices[index+1])
        xr_slice = xr_slice.where(xr_array.notnull(),0)
        xr_slice.rio.to_raster(slice_output_path)
        
        print(f"finished_{slices[index]:.0e}_to_{slices[index+1]:.0e}")
        if plot:
            plt.figure()
            xr_slice.plot()
            plt.title(f"{slices[0]:.0e}_to_{slices[index+1]:.0e}")
    
    return {"slice_indices": slices, "slice_paths": slice_paths}
    
# convert array of slice files
def convert_slices (
    slice_path_array, output_folder,slice_indices,
    scaling_factor = 0.024,
    template_path = r"template_data/dust_emissions_05.20210906.nc"):
    
    resolutions = ["fulrez", "lowres", "cf_nc"]
    epsg_s = [5070, 4326]
    
    slices = slice_indices
    
    output_nc_paths = []
    all_paths = []
    range_key = []
    
    for index, slice_path in enumerate(slice_path_array):
        
        slice_range = f'{slices[index]:.1e}_to_{slices[index+1]:.1e}'
        
        fulrez_4326_path = f"{output_folder}raster_{resolutions[0]}_{epsg_s[1]}_{slice_range}.tif"
        fulrez_4326_path_none_with_0 = f"{output_folder}raster_{resolutions[0]}_{epsg_s[1]}_{slice_range}_none_with_0.tif"
        lowrez_4326_path = f"{output_folder}raster_{resolutions[1]}_{epsg_s[1]}_{slice_range}.tif"
        output_nc_path = f"{output_folder}raster_{resolutions[2]}_{epsg_s[1]}_{slice_range}.nc"
        input_tif_5070_path = slice_path
        tif_5070_path_none_with_0 = f"{output_folder}raster_{resolutions[0]}_{epsg_s[0]}_{slice_range}_none_with_0.tif"
        
        tif_to_nc(
            fulrez_4326_path,
            lowrez_4326_path,
            output_nc_path,
            input_tif_5070_path,
            tif_5070_path_none_with_0,
            fulrez_4326_path_none_with_0,
            template_path,
            scaling_factor=scaling_factor)
        
        output_nc_paths.append(output_nc_path)
        all_paths.append([fulrez_4326_path, fulrez_4326_path_none_with_0, lowrez_4326_path, tif_5070_path_none_with_0, output_nc_path])
        range_key.append(slice_range)
        
        print(f"finished converting {slice_range}")
        
    return range_key, output_nc_paths, all_paths

def run_stats(keys, paths_to_rasters_2d):
    for index, key in enumerate(keys):
        paths_to_rasters_1d = paths_to_rasters_2d[index]
        slice_data = []
        index_names = []
        
        for path in paths_to_rasters_1d:
            xr_array = rioxarray.open_rasterio(path)
            xr_array = xr_array.where(xr_array != xr_array._FillValue,0) #remove none values so stats can calculate ok
            
            row_stats = {
                "Total_in_slice":int(xr_array.sum()),
                "Lowest_lon":int(xr_array.x.min()),
                "Highest_lon":int(xr_array.x.max()),
                "Lowest_lat":int(xr_array.y.min()),
                "Highest_lat":int(xr_array.y.max())
            }
            
            name = path[20:len(path)-4]
            slice_data.append(row_stats)
            index_names.append(name)
            print(f'    finished taking stats for {name}')
   
        slice_df = pd.DataFrame(slice_data, index = index_names)
        slice_df.to_csv(f"slice_stats/{key}.csv")
        print(f"finished taking stats for {key}")
    

# get stats