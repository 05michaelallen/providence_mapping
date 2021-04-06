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
import cartopy 
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import cmocean

os.chdir("/home/vegveg/providence_mapping/code")

# =============================================================================
# functions
# =============================================================================
def zonalmetric(shp, rast, classids, stats):
    """
    Take zonal stats of raster based on polygons in a shapefile. Note: this is meant for classifications.

    Parameters
    ----------
    shp : str
        relative path to shapefile
    rast : str
        relative path to raster
    classids : list of ints
        class ids to query
    stats : list of str
        stats to perform, red into rasterstats

    Returns
    -------
    rs_agg : list
        list of numpy arrays with size classids x stats

    """
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
# run zonal stats
# =============================================================================
shp = "../data/Los Angeles Neighborhood Map/geo_export_3b62483e-90a4-4acf-9359-a830f3086143_utm_wgs84.shp"
#shp = "../data/acs_2016_5yr_bg_06_california_polygons_joined_clip.geojson"
shp = "../data/Median_Household_Income_(2016)/Median_Household_Income_(2016)_wgs84_buffer.shp"
rast = "../data/lariac_10m_fractions_noclip.tif"
rast_link = rio.open(rast)
classids = [1, 2]
stats = ['mean']

# run zonal stats
rs_agg = zonalmetric(shp, rast, classids, stats)
rs_agg = np.array(rs_agg).T * 100

# merge output w/ shapefile
s = gpd.read_file(shp)
sr = pd.concat([s,
                pd.DataFrame(rs_agg, columns = ['tree', 'grass'])],
               axis = 1)

# create total veg column
sr['veg'] = sr['tree'] + sr['grass']

# adjust income
sr['MHI2016'] = sr['MHI2016'] / 1000

# =============================================================================
# plot
# =============================================================================
col = 'MHI2016'
col = 'veg'
labels = True
roffset = 25000
coffset = 20000
projection = ccrs.UTM(11)
bbox = [rast_link.bounds[0] + roffset - 16000, rast_link.bounds[2] - roffset - 16000,
        rast_link.bounds[1] + coffset, rast_link.bounds[3] - coffset]

fig, ax0 = plt.subplots(1, 1, figsize = [12, 9],
                        #constrained_layout = True,
                        subplot_kw = dict(projection = projection))

ax0.set_extent(bbox, projection)

ax0.add_feature(cartopy.feature.OCEAN, zorder = 0)
ax0.add_feature(cartopy.feature.LAND, zorder = 0)
#ax0.add_feature(cartopy.feature.COASTLINE, zorder = 10, linewidth = 4)

# add plot
if col == 'veg':
    s0 = sr.plot(column = col, vmin = 0, vmax = 40, ax = ax0, cmap = plt.cm.Greens,
                 edgecolor = '0.5', linewidth = 0.2)
    label = "Green Cover, %"
else:
    s0 = sr.plot(column = col, vmin = 0, vmax = 150, ax = ax0, cmap = cmocean.cm.deep_r, #plt.cm.cividis,
                 edgecolor = '0.5', linewidth = 0.2, missing_kwds = dict(color = 'lightgrey'))
    label = "2016 Median Household Income (*1000)"

patch_col = ax0.collections[0]
cb = fig.colorbar(patch_col, ax = ax0, shrink = 0.6)
cb.ax.tick_params(labelsize = 12)
cb.set_label(label = label, size = 12, rotation = 90)

scale_bar(ax0, 11, 3, location = [0.05, 0.02])

ax0.set_xticks(np.linspace(np.round(bbox[0]+5000, -4), np.round(bbox[1]-1000, -4), 4))
ax0.set_yticks(np.linspace(np.round(bbox[2], -4), np.round(bbox[3], -4), 4))
ax0.set_xlabel("Easting, m")
ax0.set_ylabel("Northing, m")
ax0.tick_params(top = True, right = True, zorder = 0)
ax0.grid(zorder = 0)

# add labels
if labels:
    s_lab = gpd.read_file("../data/Los Angeles Neighborhood Map/geo_export_3b62483e-90a4-4acf-9359-a830f3086143_utm_wgs84.shp")
    for i in range(len(s_lab)):
            sri = s_lab.iloc[i,:]
            if float(sri['sqmi']) > 3:
                x = sri['geometry'].centroid.x
                y = sri['geometry'].centroid.y
                if (x > bbox[0]) and (x < bbox[1]) and (y > bbox[2]) and (y < bbox[3]) and (sri['name'] != "Boyle Heights") and (sri['name'] != "Commerce") and (sri['name'] != "Alhambra") and (sri['name'] != "Encino") and (sri['name'] != "Hollywood Hills West") and (sri['name'] != "Sherman Oaks") and (sri['name'] != "South Pasadena"):
                    ax0.text(x, y, sri['name'], va = 'center', ha = 'center', fontsize = 10, color = '0',
                             path_effects = [PathEffects.withStroke(linewidth = 1.5, foreground = "0.9", alpha = 1)])

plt.tight_layout()
plt.savefig("../plots/" + col + "_tract_zoom.png", dpi = 800)