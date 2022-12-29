from occurrences import Occurrences, Raster
import pandas as pd
import geopandas as gpd
import rioxarray
from shapely import Point
import numpy as np
import xarray
import matplotlib.pyplot as plt

raster = rioxarray.open_rasterio('../data/wc2.1_30s_bio_1.tif')

# print(raster[1, 2])

data = pd.read_csv('../data/occurrence.txt', sep="\t", low_memory=False)
geometry = [Point(xy) for xy in zip(data.decimalLongitude, data.decimalLatitude)]
# geometry = [x.buffer(0.2) for x in geometry]
data = gpd.GeoDataFrame(data, crs="EPSG:4326", geometry=geometry)

data.within(raster)

# # print(data.loc[0:1, :])
# occs = Occurrences(data.loc[:5000, :])
# thinned = occs.thin(resolution = 10, thin_proportion=0.95)
# # # print(thinned)
# # x = occs.within(raster)
# # y = x[~x.geometry.is_empty]
# # print(y)
# # print(plt.show(thinned.plot))
# # print(type(thinned))
# thinned.plot()

# rast[x, y]

# for x in rast
# z = [x for x in y]
# print(y[0,0])
# print(thinned)

# print(range(len(rast.coords['longitude'])))
# print(y.rio.total_bounds)

# print(y[1:3,])

# print(ds_1)
#
# print(occs.within(y))

# print(y.coords)

# print(raster.coords)

# print(y)
# print(y.shape)

# y.intersection(x.loc[0:3, :])
# occs.any_within(y)
# print(y)
# r = Raster(raster)
# print(r.contains(data))

# print(y.shape)
# print(y.crs)