# %%
import xarray as xr
import os, regionmask
import geopandas as gpd
import pandas as pd

# read in shapefile and all the netcdf files 
shp = gpd.read_file(r'ices_ecoregions.gpkg')
shp = shp.to_crs(4326)

files = [os.path.join('.', 'reanalysis', i) for i in os.listdir(r'.\reanalysis') if '2010' in i]
dx = xr.open_mfdataset(r'.\reanalysis\*.nc')
# %%
# iterate over each basin
for index, row in shp.iterrows():
    if index == 1: 
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
        
        df = basinAvg.to_dataframe().transpose().reset_index().drop('index', axis=1)
        gdf = gpd.GeoDataFrame(row).transpose()

        basindf = gpd.GeoDataFrame(
            pd.concat([gdf.reset_index(), df.reset_index()], axis=1))
        basindf.crs = 'EPSG:4326'

        name = row['Ecoregion'].replace(' ', '_')

        if not os.path.exists(os.path.join('.', 'shapefile')):
            os.makedirs(os.path.join('.', 'shapefile'))

        #cannot have ints as column headers so convert to str before saving
        basindf.columns = [str(col) for col in basindf.columns]
        basindf.drop(columns=['index'], inplace=True)

        #add unit
        columns=basindf.columns.tolist()
        columns.insert(6, 'Unit')
        basindf_n = basindf.reindex(columns=columns)
        basindf_n['Unit'] = 'sea_surface_height_in_meters'

        basindf_n.to_file(os.path.join('.','shapefile', name + '.shp'))

        if not os.path.exists(os.path.join('.', 'csv')):
            os.makedirs(os.path.join('.', 'csv'))

        dft = df.transpose()
        dft.rename(columns={dft.columns[0]: 'sea_surface_height_in_meters'}).to_csv(os.path.join('.', 'csv', name + '.csv'))

        

        
    
# %%
