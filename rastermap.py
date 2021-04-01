#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 10:32:57 2021

@author: vegveg
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio as rio 
import cartopy 
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

os.chdir("/home/vegveg/providence_mapping/code")

# =============================================================================
# functions
# =============================================================================
def scale_bar(ax, zone, length, location = (0.85, 0.05), linewidth = 2, fontsize = 12):
    """
    ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar.
    """

    #Get the extent of the plotted area in coordinates in metres
    tmc = ccrs.UTM(zone, southern_hemisphere = False)
    x0, x1, y0, y1 = ax.get_extent(tmc)
    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length:
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)        
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length)

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform = tmc, color = '0', linewidth = linewidth, solid_capstyle = "butt")
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform = tmc, color = '0', fontsize = fontsize,
            horizontalalignment = 'center', verticalalignment = 'bottom')

# =============================================================================
# import and process
# =============================================================================
shp = "../data/Los Angeles Neighborhood Map/geo_export_3b62483e-90a4-4acf-9359-a830f3086143_utm_wgs84.shp"
rast = "../data/lariac_20m_fractions_noclip.tif"
rast_link = rio.open(rast)
rast_read = rast_link.read()

# =============================================================================
# plots
# =============================================================================
labels = False
roffset = 25000
coffset = 20000
projection = ccrs.UTM(11)
bbox = [rast_link.bounds[0], rast_link.bounds[2],
        rast_link.bounds[1], rast_link.bounds[3]]
bbox_zoom = [rast_link.bounds[0] + roffset - 16000, rast_link.bounds[2] - roffset - 16000,
             rast_link.bounds[1] + coffset, rast_link.bounds[3] - coffset]

fig, ax0 = plt.subplots(1, 1, figsize = [6.5, 6.5],
                        #constrained_layout = True,
                        subplot_kw = dict(projection = projection))

ax0.set_extent(bbox, projection)

ax0.add_feature(cartopy.feature.OCEAN, zorder = 0)
ax0.add_feature(cartopy.feature.LAND, zorder = 0)

# add plot
p0 = ax0.imshow(rast_read[0] * 100, vmin = 0, vmax = 50, cmap = plt.cm.Greens,
                extent = bbox, origin = 'upper', zorder = 0)
ax0.set_xlim(bbox_zoom[0], bbox_zoom[1])
ax0.set_ylim(bbox_zoom[2], bbox_zoom[3])
ax0.tick_params(top = True, right = True, zorder = 00)
ax0.grid(zorder = 0)

# labels
ax0.set_xticks(np.linspace(np.round(bbox_zoom[0], -4), np.round(bbox_zoom[1]-1000, -4), 4))
ax0.set_yticks(np.linspace(np.round(bbox_zoom[2], -4), np.round(bbox_zoom[3], -4), 4))
ax0.set_xlabel("Easting, m")
ax0.set_ylabel("Northing, m")
cb = fig.colorbar(p0, ax = ax0, shrink = 0.6,
                  ticks = [0, 10, 20, 30, 40, 50])
cb.ax.tick_params(labelsize = 12)
cb.set_label(label = "Tree Cover, %", size = 12, rotation = 90)

scale_bar(ax0, 11, 3, location = [0.05, 0.02])

#plt.tight_layout()
plt.savefig("../plots/tree_rast_zoom.png", dpi = 800)