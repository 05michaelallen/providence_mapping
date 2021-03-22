# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio as rio 
import rasterio.mask
import matplotlib.pyplot as plt

os.chdir("/Users/mallen/Documents/q11/providence_mapping/code")

# =============================================================================
# functions
# =============================================================================
def join(shp, join_layer):
    joined = pd.concat([shp, join_layer], axis = 1)
    return joined

# =============================================================================
# main 
# =============================================================================
# get all files
shp = gpd.read_file("../data/acs_2016_5yr_bg_06_california_polygons.geojson")
join_layer1 = pd.read_csv("../data/5yr_bg_b19301e1_pcincome.csv").iloc[:,0] 
join_layer2 = pd.read_csv("../data/5yr_bg_b01001e1_population.csv").iloc[:,0]
clipper = gpd.read_file("../data/la_appeears_request.shp").to_crs(epsg = 32611)

# join the loaded layers
joined = join(shp, join_layer1)
joined = join(joined, join_layer2)

# reproject and clip
joined_reproj = joined.to_crs(epsg = 32611)
joined_reproj_clip = gpd.clip(joined_reproj, clipper)

# output
joined_reproj_clip.to_file("../data/acs_2016_5yr_bg_06_california_polygons_joined_clip.geojson", driver = "GeoJSON")
