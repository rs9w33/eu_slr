# %%
import xarray as xr
import os
import geopandas as gpd
import regionmask
import matplotlib.pyplot as plt

# %% read in shapefile and all the netcdf files 
shp = gpd.read_file(r'ices_ecoregions.gpkg')
shp = shp.to_crs(4326)

#exploded=shp.explode(index_parts=False)
#shp['mean_sealevel'] =''

#%%
files = [os.path.join('.', 'reanalysis', i) for i in os.listdir(r'.\reanalysis') if '2010' in i]
dx = xr.open_mfdataset(r'.\reanalysis\*.nc')

# %% mask the  netcdf files 
shp_mask = regionmask.mask_geopandas(
                shp, 
                dx.station_x_coordinate, 
                dx.station_y_coordinate,
                wrap_lon=180)
ds_mask = dx.where(~shp_mask.isnull(), drop=True)

# Calculate the mean along the time dimension pver the masked dataset
yearly_mean = ds_mask['waterlevel'].groupby('time.year').mean(dim='time')

# %% plotting to make sure we have the correct data
y = yearly_mean['station_y_coordinate'].values
x = yearly_mean['station_x_coordinate'].values
# Create a scatter plot
plt.scatter(x, y, color='red', marker='o')

# %%
for index, row in shp.iterrows():
    # Get the polygon geometry
    polygon = gpd.GeoSeries(row.geometry, crs='EPSG:4326')

    regmask = regionmask.mask_geopandas(
        polygon,
        dx.station_y_coordinate,
        dx.station_x_coordinate,
        wrap_lon=180
    )
    basin = dx.where(~regmask.isnull(),drop=True)
    basinAvg = basin.waterlevel.groupby('time.year').mean(dim='time').mean(dim='stations')
    
    clipped_data = yearly_mean.sel(latitude=slice(polygon.bounds.miny, polygon.bounds.maxy), 
                          longitude=slice(polygon.bounds.minx, polygon.bounds.maxx))

    # Calculate average for a specific variable (e.g., 'temperature')
    avg = clipped_data['waterlevel'].mean(dim='time')

    shp['mean_waterlevel'][index]= avg
# %%
