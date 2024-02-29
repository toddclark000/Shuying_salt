from tif_to_nc import tif_to_nc

fulres_4326_path = r"intermediate_data/fulres_4326.tif"
lowres_4326_path = r"intermediate_data/lowres_4326.tif"
output_nc_path = r"output_data/output_nc.nc"
input_tif_5070_path = r"input_data/1992_2015/1992.tif"
template_path = r"template_data/dust_emissions_05.20210906.nc"
scaling_factor = 0.024 #may need to be switched to 0.025

tif_to_nc(fulres_4326_path, lowres_4326_path, output_nc_path, input_tif_5070_path, template_path, scaling_factor=scaling_factor)