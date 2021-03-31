#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 10:16:57 2021

@author: vegveg
"""

import os
import numpy as np
import pandas as pd
import rasterstats as rs
import geopandas as gpd
import rasterio as rio 
#import cartopy 
#import cartopy.crs as ccrs
import matplotlib.pyplot as plt

os.chdir("/home/vegveg/providence_mapping/code")

# =============================================================================
# functions
# =============================================================================
def zonalmetric(shp, rast, classids, stats):
    r = rio.open(rast)
    s = gpd.read_file(shp)
    rs_agg = []
    for stat in stats:
        for classid in classids:
            rs_stat = rs.zonal_stats(s, r.read()[classid,:,:],
                                     affine = r.transform,
                                     stats = stat)
            rs_stat = [np.asarray([x[attribute] for x in rs_stat]) for attribute in [stat]][0]
            rs_stat[rs_stat == np.array(None)] = np.nan
            rs_agg.append(rs_stat.astype(float))
        return rs_agg
    
# =============================================================================
# run zonal stats
# =============================================================================
shp = "../data/Los Angeles Neighborhood Map/geo_export_3b62483e-90a4-4acf-9359-a830f3086143_utm_clip_clean.shp"
shp = "../data/acs_2016_5yr_bg_06_california_polygons_joined_clip.geojson"
rast = "../data/lariac_10m_fractions.tif"
classids = [1, 2]
stats = ['mean']

# run zonal stats
rs_agg = zonalmetric(shp, rast, classids, stats)
rs_agg = np.array(rs_agg).T

# merge output w/ shapefile
s = gpd.read_file(shp)
sr = pd.concat([s,
                pd.DataFrame(rs_agg, columns = ['tree', 'grass'])],
               axis = 1)

# create total veg column
sr['veg'] = sr['tree'] + sr['grass']

# plot
p = sr.plot('veg')
plt.colorbar(p)


# =============================================================================
# plot
# =============================================================================
