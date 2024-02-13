"""Creates a HEMCO emissions file for use in GEOS-Chem with a single 
point-source of H2O2."""
#imports
import xarray as xr
import numpy as np
import os
from datetime import datetime, timedelta
from geopy.distance import geodesic
import pytz

#position values for each added source as a list of tuples of `(lat, lon)`
target_latlons = [(50.5, -110.5), (31.4, -87.7), (35.3, -119.1)]

#time zones for each tower (recognised by pytz package). Must be in the same 
#order as the lat/lons
tower_tzs = ["America/Edmonton", "America/Chicago", "America/Los_Angeles"]

#start and end times for the H2O2 emission (in hours of the day local time i.e. 18 = 6pm)
ton = 8
toff = 18

#path to save the output to
out_path = "Diurnal_H2O2_Emission.nc"

#Emission flux of H2O2 from each tower in kg/s (this will be converted to kg/m2/s by 
#calculating the area of each specified grid box and dividing the emission rate 
#by the area)
emiss_rate = 0.612

#horizontal resolution for the final emissions file
target_resolution = 0.1

#NaN fill value for the emissions file
nan_val = 1E20
###############################################################################
#make time index
times = np.arange(0.0,24.0,1.0)
#reference date for the data. This shouldn't matter too much as GEOS-Chem will 
#reuse the diurnal data for each day of the model if the right settings are 
#speciefied in HEMCO_Config.rc
ref_datetime = "2016-01-01 00:00:00"

time_array = xr.DataArray(data= times, name = "time", dims = "time",
                          attrs = {"long_name":"time",
                                   "units" : f"hours since {ref_datetime}",
                                   "calendar" : "standard",
                                   "axis" : "T"})

#make lat and lon indices
lats = np.arange(-90+(target_resolution/2), 90+(target_resolution/2), 
                 target_resolution, dtype = float)
lat_array = xr.DataArray(data= lats, name = "lat", dims = "lat",
                          attrs = {"long_name":"Latitude",
                                   "units" : "degrees_north",
                                   "axis" : "Y"})
lons = np.arange(-180+(target_resolution/2), 180+(target_resolution/2), 
                 target_resolution, dtype = float)
lon_array = xr.DataArray(data= lons, name = "lon", dims = "lon",
                          attrs = {"long_name":"Longitude",
                                   "units" : "degrees_east",
                                   "axis" : "X"})

coords_dict = {"time" : time_array,
               "lat" : lat_array,
               "lon" : lon_array}

#make the emissions data
empty_data = np.empty((len(times), len(lats), len(lons)))
empty_data.fill(float(0.0))

empty_array = xr.DataArray(data= empty_data, name = "H2O2", dims = ["time",
                                                                    "lat", 
                                                                    "lon"],
                          attrs = {"long_name":"H2O2",
                                   "units" : "kg/m2/s",
                                   "_FillValue" : nan_val,
                                   "missing_value" : nan_val})

                                                                     
#Set up the COARDS-compliant netCDF dataset
emiss_ds = xr.Dataset(data_vars = {"H2O2" : empty_array}, 
                      coords=coords_dict,
                      attrs = {"Title" : "Emissions data simulating H2O2 towers for the removal of CH4.",
                               "Contact" : "Alfred Mayhew (a.mayhew@utah.edu)",
                               "References" : "",
                               "Conventions" : "COARDS",
                               "Filename" : out_path.split(os.sep)[-1],
                               "History" : datetime.strftime(datetime.now(), 
                                                             "%a %b %d %H:%M:%S %Y"),
                               "Format" : "NetCDF-4",
                               "Grid" : f"GEOS_{target_resolution}x{target_resolution}",
                               "Delta_Lon" : target_resolution,
                               "Delta_Lat" : target_resolution,
                               "SpatialCoverage" : "global"})

#add emissions in specified locations
def get_emiss_list(time_zone):
    "Determines when to turn emissions on or off in local time"
    tz_obj = pytz.timezone(time_zone)
    ref_dt = datetime.strptime(ref_datetime, "%Y-%m-%d %H:%M:%S")  
    ref_dt = tz_obj.localize(ref_dt)
    
    local_ton = ref_dt + timedelta(hours = ton)
    local_toff = ref_dt + timedelta(hours = toff)
    
    utc_ton = local_ton.astimezone(pytz.UTC).time().hour
    utc_toff = local_toff.astimezone(pytz.UTC).time().hour
    
    #make a list of emission rates, accounting for running over midnight
    if utc_toff > utc_ton:
        emiss_func = lambda x : emiss_rate if ((x >= utc_ton) and (x <= utc_toff)) else 0
    elif utc_toff < utc_ton:
        emiss_func = lambda x : emiss_rate if ((x >= utc_ton) or (x <= utc_toff)) else 0
    return [emiss_func(x) for x in times]
    
    
for i, (tlat, tlon) in enumerate(target_latlons):
    diurnal_vals = get_emiss_list(tower_tzs[i])
    sel_latlon = emiss_ds.sel(lat = tlat, lon = tlon, method = "nearest")
    upper_latlon = emiss_ds.sel(lat = tlat + target_resolution, 
                                lon = tlon + target_resolution, 
                                method = "nearest")
    sel_lat = sel_latlon.lat.values.tolist()
    sel_lon = sel_latlon.lon.values.tolist()
    up_lon = upper_latlon.lon.values.tolist()
    
    #calculate the emission flux by dividing the rate by the area of the target
    #grid cell (note that I'm ignoring curvature in the area calculation 
    #because the selected grid cell should be very small)
    edge_dist = geodesic((sel_lat, sel_lon), (sel_lat, up_lon))
    cell_area = edge_dist.m**2
    
    emiss_ds["H2O2"].loc[dict(lon=sel_lon, 
                              lat=sel_lat)] = [x/cell_area for x in diurnal_vals]

#fill nan values with the nan value
emiss_ds = emiss_ds.fillna(nan_val)

#save the dataset to a netcdf
emiss_ds.to_netcdf(out_path, format = "NETCDF4")

print(f"""The following line will need to be added to the HEMCO_Config file (after the file has been copied to the GC run directory)):
`0 Tower_H2O2 $CFDIR/{out_path.split(os.sep)[-1]} H2O2 2016/1/1/0-23 C xyL=HEIGHT_OF_TOWERm kg/m2/s H2O2 - 2 50`""")
