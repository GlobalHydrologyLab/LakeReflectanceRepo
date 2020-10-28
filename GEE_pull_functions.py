#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Google Earth Engine Reflectance Pull Functions
Created on Mon Apr  9 14:24:13 2018
@author: simontopp
"""

# Add filler panchromatic band to landsat 5 images.

###These are functions for unpacking the bit quality assessment band for TOA  
def Unpack(bitBand, startingBit, bitWidth):
  #unpacking bit bands
  #see: https://groups.google.com/forum/#!starred/google-earth-engine-developers/iSV4LwzIW7A
  return (ee.Image(bitBand)\
  .rightShift(startingBit)\
  .bitwiseAnd(ee.Number(2).pow(ee.Number(bitWidth)).subtract(ee.Number(1)).int()))
  
def UnpackAll(bitBand, bitInfo):
  unpackedImage = ee.Image.cat([Unpack(bitBand, bitInfo[key][0], bitInfo[key][1]).rename([key]) for key in bitInfo])
  return unpackedImage

def clipImage(image):
  return image.clip(lake.geometry())

## These functions all go into calculating the USGS Dynamic water

# def AddFmask(image):

#     bitInfo = {
#     'Cloud': [5, 1],
#     'CloudShadow': [3, 1], 
#     'SnowIce': [4, 1],
#     'Water': [2, 1]
#     }
      
#     temp = UnpackAll(image.select(['pixel_qa']), bitInfo)
      
#     fmask = (temp.select(['Water']).rename(['fmask'])
#     .where(temp.select(['SnowIce']), ee.Image(3))
#     .where(temp.select(['CloudShadow']), ee.Image(2))
#     .where(temp.select(['Cloud']), ee.Image(4))
#     .mask(temp.select(['Cloud']).gte(0))) 
#       #mask the fmask so that it has the same footprint as the quality (BQA) band
#     return(image.addBands(fmask))

def AddFmask(image):
    qa = image.select('pixel_qa')
    water = qa.bitwiseAnd(1 << 2)
    cloud = qa.bitwiseAnd(1 << 5)
    snow = qa.bitwiseAnd(1 << 4)
    cloudshadow = qa.bitwiseAnd(1 << 3)
    
    fmask = (water.gt(0).rename(['fmask'])
    .where(snow.gt(0), ee.Image(3))
    .where(cloudshadow.gt(0), ee.Image(2))
    .where(cloud.gt(0), ee.Image(4))
    .updateMask(qa.gte(0)))
    #mask the fmask so that it has the same footprint as the quality (BQA) band
    return image.addBands(fmask)
    

def Mndwi(image):
  return image.normalizedDifference(['Green', 'Swir1']).rename('mndwi')
  
def Mbsrv(image):
  return image.select(['Green']).add(image.select(['Red'])).rename('mbsrv')
  
def Mbsrn(image):
  return image.select(['Nir']).add(image.select(['Swir1'])).rename('mbsrn')

def Ndvi(image):
  return image.normalizedDifference(['Nir', 'Red']).rename('ndvi')

def Awesh(image):
  return (image.addBands(Mbsrn(image))
  .expression('Blue + 2.5 * Green + (-1.5) * mbsrn + (-0.25) * Swir2', {
    'Blue': image.select(['Blue']),
    'Green': image.select(['Green']),
    'mbsrn': Mbsrn(image).select(['mbsrn']),
    'Swir2': image.select(['Swir2'])
    }))


## The DSWE Function itself    
def Dswe(i):
  mndwi = Mndwi(i)
  mbsrv = Mbsrv(i)
  mbsrn = Mbsrn(i)
  awesh = Awesh(i)
  swir1 = i.select(['Swir1'])
  nir = i.select(['Nir'])
  ndvi = Ndvi(i)
  blue = i.select(['Blue'])
  swir2 = i.select(['Swir2'])

  t1 = mndwi.gt(0.124)
  t2 = mbsrv.gt(mbsrn)
  t3 = awesh.gt(0)
  t4 = (mndwi.gt(-0.44)
  .And(swir1.lt(900))
  .And(nir.lt(1500))
  .And(ndvi.lt(0.7)))
  t5 = (mndwi.gt(-0.5)
  .And(blue.lt(1000))
  .And(swir1.lt(3000))
  .And(swir2.lt(1000))
  .And(nir.lt(2500)))
  
  t = t1.add(t2.multiply(10)).add(t3.multiply(100)).add(t4.multiply(1000)).add(t5.multiply(10000))
  
  noWater = (t.eq(0)
  .Or(t.eq(1))
  .Or(t.eq(10))
  .Or(t.eq(100))
  .Or(t.eq(1000)))
  hWater = (t.eq(1111)
  .Or(t.eq(10111))
  .Or(t.eq(11011))
  .Or(t.eq(11101))
  .Or(t.eq(11110))
  .Or(t.eq(11111)))
  mWater = (t.eq(111)
  .Or(t.eq(1011))
  .Or(t.eq(1101))
  .Or(t.eq(1110))
  .Or(t.eq(10011))
  .Or(t.eq(10101))
  .Or(t.eq(10110))
  .Or(t.eq(11001))
  .Or(t.eq(11010))
  .Or(t.eq(11100)))
  pWetland = t.eq(11000)
  lWater = (t.eq(11)
  .Or(t.eq(101))
  .Or(t.eq(110))
  .Or(t.eq(1001))
  .Or(t.eq(1010))
  .Or(t.eq(1100))
  .Or(t.eq(10000))
  .Or(t.eq(10001))
  .Or(t.eq(10010))
  .Or(t.eq(10100)))
  
  iDswe = (noWater.multiply(0)
  .add(hWater.multiply(1))
  .add(mWater.multiply(2))
  .add(pWetland.multiply(3))
  .add(lWater.multiply(4)))
  
  return iDswe.rename('dswe')


## Calculuate hillshades to correct DWSE
def CalcHillShades(image, geo):
  MergedDEM = ee.Image("users/eeProject/MERIT").clip(geo.buffer(300))
  
  hillShade = (ee.Terrain.hillshade(MergedDEM, ee.Number(image.get('SOLAR_AZIMUTH_ANGLE')),
  image.get('SOLAR_ZENITH_ANGLE')).rename(['hillShade']))
  return hillShade
  
  ## Calculuate hillshadow to correct DWSE
def CalcHillShadows(image, geo):
  MergedDEM = ee.Image("users/eeProject/MERIT").clip(geo.buffer(3000))
  hillShadow = (ee.Terrain.hillShadow(MergedDEM, ee.Number(image.get('SOLAR_AZIMUTH_ANGLE')),
  ee.Number(90).subtract(image.get('SOLAR_ZENITH_ANGLE')), 30).rename(['hillShadow']))
  return hillShadow


## Buffer the lake sites
def dpBuff(i):
    return i.buffer(120)


## Remove geometries
def removeGeo(i):
    return i.setGeometry(None)
## Create water mask and extract lake medians

## Set up the reflectance pull
def RefPull(image):
    f = AddFmask(image).select('fmask')
    clouds = f.gte(2).rename('clouds')
    #cScore = clouds.reduceRegion(ee.Reducer.mean(),lake.geometry(), 30).get('fmask')
    d = Dswe(image).select('dswe')
    #h = CalcHillShades(image, tile.geometry()).select('hillShade')
    hs = CalcHillShadows(image, tile.geometry()).select('hillShadow')
    dswe3 = d.eq(3).rename('dswe3').selfMask().updateMask(clouds.eq(0)) 
    dummy = (image.select(['Blue'],['dswe1'])
        .updateMask(clouds.eq(0)).updateMask(d.eq(1)))
    
    #imageSD = image.select(['Blue','Green','Red','Nir'],['BlueSD','GreenSD','RedSD','NirSD'] )
   
    pixOut = (image.addBands(hs)
              .addBands(image.select(['Nir'],['NirSD']))
              .updateMask(d.eq(1))
              .updateMask(clouds.eq(0))
              .addBands(dswe3)
              .addBands(dummy)
              .addBands(clouds))
    
    combinedReducer = (ee.Reducer.median().unweighted()
    .forEachBand(pixOut.select(['Aerosol','Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', 'TIR1','TIR2','pixel_qa', 'hillShadow']))
    .combine(ee.Reducer.stdDev().unweighted().forEachBand(pixOut.select(['NirSD'])), 'sd_', False)
    .combine(ee.Reducer.count().unweighted().forEachBand(pixOut.select(['dswe3', 'dswe1'])), 'pCount_', False)
    .combine(ee.Reducer.mean().unweighted().forEachBand(pixOut.select(['clouds'])), 'cScore_', False))
         
    
    # Collect median reflectance and occurance values
    # Make a cloud score, and get the water pixel count
    lsout = pixOut.reduceRegions(lakes, combinedReducer, 30)
    
    out = lsout.map(removeGeo)
    
    return out

##Function for limiting the max number of tasks sent to
#earth engine at one time to avoid time out errors

def maximum_no_of_tasks(MaxNActive, waitingPeriod):
  ##maintain a maximum number of active tasks
  time.sleep(10)
  ## initialize submitting jobs
  ts = list(ee.batch.Task.list())

  NActive = 0
  for task in ts:
       if ('RUNNING' in str(task) or 'READY' in str(task)):
           NActive += 1
  ## wait if the number of current active tasks reach the maximum number
  ## defined in MaxNActive
  while (NActive >= MaxNActive):
      time.sleep(waitingPeriod) # if reach or over maximum no. of active tasks, wait for 2min and check again
      ts = list(ee.batch.Task.list())
      NActive = 0
      for task in ts:
        if ('RUNNING' in str(task) or 'READY' in str(task)):
          NActive += 1
  return()
    
