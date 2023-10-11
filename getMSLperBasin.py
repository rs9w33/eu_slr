# %%
import xarray as xr
import os, regionmask
import geopandas as gpd
import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
# %%
# read in shapefile and all the netcdf files 
shp = gpd.read_file(r'ices_ecoregions.gpkg')
shp = shp.to_crs(4326)

 #%%
dx = xr.open_mfdataset(r'.\reanalysis\*.nc')
mdx = xr.open_mfdataset(r'.\slr\*.nc')
# %%
# iterate over each basin
for index, row in shp.iterrows():
    # Get the polygon geometry
    polygon = gpd.GeoSeries(row.geometry, crs='EPSG:4326')

    #processing GTSM
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

    #calculate std
    #basinStd = basin.waterlevel.groupby('time.year').std(dim='time').std(dim='stations')
    #dfstd = basinStd.to_dataframe().transpose().reset_index().drop('index', axis=1)
    #dfs = dfstd.transpose()
    #dfs.rename(columns={dft.columns[0]: 'standard_deviation_ssh'}, inplace=True)
    #result = pd.concat([dft, dfs], axis=0)

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

    basindf_n.to_file(os.path.join('./shapefile/', name + '_gtsm.shp'))

    if not os.path.exists(os.path.join('.', 'csv')):
        os.makedirs(os.path.join('.', 'csv'))

    #------saving GTSM
    dft = df.transpose()
    dft.rename(columns={dft.columns[0]: 'sea_surface_height_in_meters'}, inplace=True)
    dft.to_csv(os.path.join('./csv/', name + '_gtsm.csv'))

    #------Calculate regression
    print('calculation regression ssh ', name)
    dft['years'] = dft.index
    # Create the independent variable (X)
    X = dft[['years']]
    # Create an instance of the LinearRegression model
    model = LinearRegression()
    # Fit the model
    model.fit(X, dft['sea_surface_height_in_meters'])
    # Generate predictions
    predictions = model.predict(X)

    if not os.path.exists(os.path.join('.', 'png')):
        os.makedirs(os.path.join('.', 'png'))

    plt.figure()
    plt.plot(dft['years'], dft['sea_surface_height_in_meters'], label='Data')
    plt.plot(dft['years'], predictions, label='Regression Line', color='red')
    plt.xlabel('Years')
    plt.ylabel('Sea Surface Height in Meters')
    plt.title(name.replace('_', ' ')+ ': Sea Surface Height in Meters')
    plt.legend()
    plt.savefig(os.path.join('./png/', name + '_gtsm.png'))  # Save the figure

    #processing dewi's MSL data
    mregmask = regionmask.mask_geopandas(
        polygon,
        mdx.station_y_coordinate,
        mdx.station_x_coordinate,
        wrap_lon=180
    )
    mbasin = mdx.where(~mregmask.isnull(),drop=True)
    mbasinAvg = mbasin.mean_sea_level.groupby('time.year').mean(dim='time').mean(dim='stations')
    #mbasinStd = mbasin.waterlevel.groupby('time.year').std(dim='time').std(dim='stations')
    
    mdf = mbasinAvg.to_dataframe().transpose().reset_index().drop('index', axis=1)
    #mdfstd = mbasinAvg.to_dataframe().transpose().reset_index().drop('index', axis=1)
    mgdf = gpd.GeoDataFrame(row).transpose()

    mbasindf = gpd.GeoDataFrame(
        pd.concat([mgdf.reset_index(), mdf.reset_index()], axis=1))
    mbasindf.crs = 'EPSG:4326'

    mname = row['Ecoregion'].replace(' ', '_')

    if not os.path.exists(os.path.join('.', 'shapefile')):
        os.makedirs(os.path.join('.', 'shapefile'))

    #cannot have ints as column headers so convert to str before saving
    mbasindf.columns = [str(col) for col in mbasindf.columns]
    mbasindf.drop(columns=['index'], inplace=True)

    #add unit
    mcolumns=mbasindf.columns.tolist()
    mcolumns.insert(6, 'Unit')
    mbasindf_n = mbasindf.reindex(columns=mcolumns)
    mbasindf_n['Unit'] = 'mean_sea_level_in_meters'

    mbasindf_n.to_file(os.path.join('./shapefile/', mname + '_sea_level.shp'))

    if not os.path.exists(os.path.join('.', 'csv')):
        os.makedirs(os.path.join('.', 'csv'))

    #------saving MSL
    mdft = mdf.transpose()
    mdft.rename(columns={mdft.columns[0]: 'mean_sea_level_in_meters'}, inplace=True)
    mdft.to_csv(os.path.join('./csv/', mname + '_sea_level.csv'))

    #------Calculate regression
    print('calculation regression dewis MSL data', name)
    mdft['years'] = mdft.index
    # Create the independent variable (X)
    X = mdft[['years']]
    # Create an instance of the LinearRegression model
    mmodel = LinearRegression()
    # Fit the model
    mmodel.fit(X, mdft['mean_sea_level_in_meters'])
    # Generate predictions
    mpredictions = mmodel.predict(X)

    if not os.path.exists(os.path.join('.', 'png')):
        os.makedirs(os.path.join('.', 'png'))

    plt.figure()
    plt.plot(mdft['years'], mdft['mean_sea_level_in_meters'], label='Data')
    plt.plot(mdft['years'], mpredictions, label='Regression Line', color='red')
    plt.xlabel('Years')
    plt.ylabel('Mean Sea Level in Meters')
    plt.title(mname.replace('_', ' ')+ ': Mean Sea Level in Meters')
    plt.legend()
    plt.savefig(os.path.join('./png/', mname + '_sea_level.png'))
# %%
