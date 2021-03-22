#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 10:13:42 2021

@author: mallen
"""

import os
import numpy as np
import rasterio as rio 

os.chdir("/home/vegveg/providence_mapping/code")
target_px_size = 10
source_px_size = 2
n_classes = 8

# =============================================================================
# 
# =============================================================================
# import lc data
lc = rio.open("../data/la_cover_urban_sub_utm_majority2m_clip_mask.tif")
lcr = lc.read()[0]

# reclassify zeros as nan
lcr = lcr.astype(np.float32)
lcr[lcr == 0] = np.nan

# create numpy array to fill
fracs = np.full([n_classes,
                 int(lcr.shape[0]/(target_px_size/source_px_size)), 
                 int(lcr.shape[1]/(target_px_size/source_px_size))], 
                -999,
                dtype = np.float32)

# grab and modify metadata
meta = lc.meta.copy()
meta.update({'dtype': 'float32',
             'count': n_classes,
             'nodata': -999})


# =============================================================================
# aggregation
# =============================================================================
total = lcr.copy()
total[total > 0] = 1
total[total < 1] = 0

# main loop to calculate fractions
for c in range(n_classes):
    c = c + 1 # classes are base 1
    print(c)
    for i in range(fracs.shape[1]):
        bi = int(i * (target_px_size/source_px_size))
        ti = int(bi + (target_px_size/source_px_size))
        for j in range(fracs.shape[2]):
            bj = int(j * (target_px_size/source_px_size))
            tj = int(bj + (target_px_size/source_px_size))
            fracs[c-1, i, j] = np.nansum(lcr[bi:ti,bj:tj] == c) / np.nansum(total[bi:ti,bj:tj])

with rio.open("../data/lariac_" + str(target_px_size) + "m_fractions.tif", 'w', **meta) as dst:
    dst.write(fracs)
