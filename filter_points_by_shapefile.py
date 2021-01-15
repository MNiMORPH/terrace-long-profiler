# script to filter the csv file to remove any points which are not within the shapefile boundary.
# Also masks the terrace relief raster to the shapefile.

import numpy as np
import pandas as pd
from glob import glob
import sys
from shapely.geometry import shape, Point
import fiona
import os
import rasterio
import rasterio.mask

#data_dir = '/home/fiona/MississippiTerraces/'
data_dir = 'D:\\MississippiTerraces\\'

subdirs = next(os.walk(data_dir))[1]

for dir in subdirs:
    print(dir)
    if 'UMV_DEM5m_10' in dir:

        # get the filename
        s = dir.split('_')
        print(s)
        fname = 'UMV_DEM5m_'+s[2]
        print(fname)

        # get the shapefile as a polygon
        c = fiona.open(data_dir+dir+'\\'+fname+'_terraces.shp')
        shp_geom = []
        new_ids = []
        i = 0
        for pol in c:
            shp_geom.append(shape(pol['geometry']))
            new_ids.append(i)
            i+=1

        # mask the relief raster to this geometry
        with rasterio.open(data_dir+dir+'\\'+fname+'_final_terrace_relief_final.bil') as src:
            out_image, out_transform = rasterio.mask.mask(src, shp_geom, crop=True)
            out_meta = src.meta

        # write to an output
        out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

        with rasterio.open(data_dir+dir+'\\'+fname+'_terrace_relief_masked.tif', "w", **out_meta) as dest:
            dest.write(out_image)

        #read in the points csv file
        pts_csv = data_dir+dir+'\\'+fname+'_final_terrace_info.csv'
        pts = []
        df = pd.read_csv(pts_csv)

        X = df.X.values
        Y = df.Y.values
        #nodes = df.node.values
        keep_index = []
        id_index = []
        # find which points are within the shapefile
        for i in range(len(X)):
            sys.stdout.write("point %d of %d       \r" %(i, len(X)))
            sys.stdout.flush()
            point = Point(X[i], Y[i])
            for j, this_shp in enumerate(shp_geom):
                if (point.within(this_shp)):
                    keep_index.append(i)
                    id_index.append(new_ids[j])

        output_df = df.loc[keep_index]
        # write the new IDs to a new column
        output_df['new_ID'] = id_index
        output_df.to_csv(data_dir+dir+'\\'+fname+'_final_terrace_info_filtered.csv', index=False)
