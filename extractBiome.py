#!/usr/bin/env python

import os
import sys
import numpy as np
import h5py as h5
import rasterio
from shapely.geometry import mapping
from rasterio.mask import mask
from shapely import geometry
from shapely.geometry import Polygon
import shapely.ops as ops
import progressbar
import pyproj 
from functools import partial
import shapely.speedups
import shapely.vectorized

print('Reading GeoTIFF GLOBCOVER')
src = rasterio.open('/home/cfranken/GLOBCOVER_L4_200901_200912_V2.3.tif')
print('Done')

f = h5.File('./TROPO_SIF_2018-03-08_ungridded.nc')
lat_bnd = f['lat_bnds'][:]
lon_bnd = f['lon_bnds'][:]
n = len(lat_bnd[0,:])
print('Read Lat/Lon bound')
geoms = [Polygon(np.vstack((lon_bnd[:,j],lat_bnd[:,j])).T.tolist()) for j in range(n)]
print('Geoms done')
# Create Biome array
biomes = np.zeros((n,24))
index_biome = np.asarray([11,14,20,30,40,50,60,70,90,100,110,120,130,140,150,160,170,180,200,210,220,230])

bar = progressbar.ProgressBar(max_value=n)
for j in range(n):
    #pp = zip(lon_bnd[:,j],lat_bnd[:,j])
    #poly = Polygon(pp)
    #geoms = [Polygon(np.vstack((lon_bnd[:,j],lat_bnd[:,j])).T.tolist())]
    #mask = shapely.vectorized.contains(geoms, x, y)
    #geom_area = ops.transform(
    #partial(pyproj.transform,pyproj.Proj(init='EPSG:4326'),pyproj.Proj(proj='aea',lat1=geoms.bounds[1],lat2=geoms.bounds[3])),geoms)

    #print(geom_area.area/1000.**2) 
    #print(geoms[0])
    try:
    	out_image, out_transform = mask(src, [geoms[j]], crop=True)
    	out_array = out_image[out_image>0]
    	n_total = len(out_array)
    	biomes[j,-1]=n_total
    	for i in range(23):
    		biomes[j,i] = np.count_nonzero(out_array == index_biome[i])/n_total*100.
    		#print(i,biomes[j,i])

    except:
        None
    if j%100==0:
    	bar.update(j)
f['biomes']=biomes
f.close()