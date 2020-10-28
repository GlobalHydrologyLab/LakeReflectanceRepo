#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:07:02 2018

@author: simontopp
"""
#%%
import time
import ee
import os
#import pandas as pd
#import feather
ee.Initialize()

#Source necessary functions.
exec(open('GEE_pull_functions.py').read())

#water = ee.Image("JRC/GSW1_1/GlobalSurfaceWater").select('occurrence').gt(80)

## Bring in EE Assets
# Deepest point for CONUS Hydrolakes from Xiao Yang
# Code available https://zenodo.org/record/4136755#.X5d54pNKgUE
dp = (ee.FeatureCollection('users/eeProject/HydroLAKES_newct_NA_e5e5da51107417d472452564cef3eb4f')
  .filterMetadata('type','equals','dp')
  .filterMetadata('distance', "greater_than", 60))

#CONUS Boundary
us = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")\
    .filterMetadata('country_na', 'equals', 'United States')

## WRS Tiles in descending (daytime) mode for CONUS
wrs = ee.FeatureCollection('users/sntopp/wrs2_asc_desc')\
    .filterBounds(us)\
    .filterMetadata('MODE', 'equals', 'D')
    
pr = wrs.aggregate_array('PR').getInfo()


#%%
dummyBands = ee.Image(-99).rename('Null_CS')\
    .addBands(ee.Image(-99).rename('Null_TIR2'))

def addDummy(i):
    return i.addBands(dummyBands)

l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
l7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')\
    .map(addDummy)
l5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')\
    .map(addDummy)

#Standardize band names between the various collections and aggregate 
#them into one image collection
bn8 = ['B1','B2','B3', 'B4', 'B5','B6','B7', 'B10','B11','pixel_qa']
bn57 = ['Null_CS', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6','Null_TIR2', 'pixel_qa']
bns = ['Aerosol','Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', 'TIR1','TIR2','pixel_qa']
  
ls5 = l5.select(bn57, bns)
ls7 = l7.select(bn57, bns)
ls8 = l8.select(bn8, bns)

ls = ee.ImageCollection(ls5.merge(ls7).merge(ls8))\
    .filter(ee.Filter.lt('CLOUD_COVER', 75))\
    .filterBounds(us)  


## Set up a counter and a list to keep track of what's been done already
counter = 0
done = []    

#%%

## In case something goofed, you should be able to just 
## re-run this chunk with the following line filtering out 
## what's already run. 
pr = [i for i in pr if i not in done]

for tiles in pr:
    tile = wrs.filterMetadata('PR', 'equals', tiles)
    # For some reason we need to cast this to a list and back to a
    # feature collection
    lakes = dp.filterBounds(tile.geometry())\
        .map(dpBuff)
        
    #lakes = ee.FeatureCollection(lakes.toList(10000))
    stack = ls.filterBounds(tile.geometry().centroid())
    out = stack.map(RefPull).flatten().filterMetadata('cScore_clouds','less_than',.5)
    dataOut = ee.batch.Export.table.toDrive(collection = out,\
                                            description = str(tiles),\
                                            folder = 'LakeReflRepo',\
                                            fileFormat = 'csv',\
                                            selectors = ['Aerosol','Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', 'TIR1','TIR2','pixel_qa', 'hillShadow', 'sd_NirSD','cScore_clouds','pCount_dswe3','pCount_dswe1','system:index','Hylak_id','distance'])
    #Check how many existing tasks are running and take a break if it's >15  
    maximum_no_of_tasks(25, 120)
    #Send next task.
    dataOut.start()
    counter = counter + 1
    done.append(tiles)
    print('done_' + str(counter) + '_' + str(tiles)))
        
#%%

