"""
/***************************************************************************
 CCDC at Point Plugin
                                 A QGIS plugin
 Continuous Change Detection Plugin
                              -------------------
        copyright            : (C) 2023-2023 by Daniel Moraes
        email                : moraesd90@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

"""
import numpy as np


def getImageCollection(name, geometry, date_start, date_end, cloud_filter, bands):
    import ee
    if name=='Sentinel-2':
        name = 'COPERNICUS/S2_SR_HARMONIZED'
    #get image collection
    img_col = ee.ImageCollection(name).filterBounds(geometry).filterDate(date_start,date_end)

    #add indices
    if 'ndvi' in list(map(lambda x: x.lower(), bands)):
        img_col = img_col.map(addNDVI)
    if 'nbr' in list(map(lambda x: x.lower(), bands)):
        img_col = img_col.map(addNBR)
    if 'evi' in list(map(lambda x: x.lower(), bands)):
        img_col = img_col.map(addEVI)

    #apply cloud filter
    if cloud_filter=='Sen2Cor':
        img_col_filtered = img_col.map(filterS2_level2A)
    elif cloud_filter=='s2cloudless':
        #get cloud probability collection
        s2_cloudprob = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY').filterBounds(geometry).filterDate(date_start,date_end)
        img_col_filtered = filterS2cloudless(img_col,s2_cloudprob)
    elif cloud_filter=='No Mask':
        img_col_filtered = img_col

    return img_col_filtered

def addNDVI(image):
    ndvi = image.normalizedDifference(['B8','B4']).multiply(10000).int16()
    return image.addBands(ndvi.rename('ndvi'))

def addNBR(image):
    nbr = image.normalizedDifference(['B8', 'B12']).multiply(10000).int16()
    return image.addBands(nbr.rename('nbr'))


def addEVI(image):
    evi = image.expression(
      '2.5 * ((NIR-RED) / (NIR + 6 * RED - 7.5* BLUE +1))', {
        'NIR':image.select('B8').divide(10000),
        'RED':image.select('B4').divide(10000),
        'BLUE':image.select('B2').divide(10000)    
          }).multiply(10000).int16()                                                        
    return image.addBands(evi.rename('evi'))



#SCL cloud/shadow filter
def filterS2_level2A(image):
    import ee
    SCL = image.select('SCL')
    mask01 = ee.Image(0).where((SCL.lt(8)).And(SCL.gt(3)),1)   #Put a 1 on good pixels
    #(SCL.gt(3),1)
    return image.updateMask(mask01)


#s2cloudless filter
def filterS2cloudless(S2SRCol, S2CloudCol):
    import ee
    CLOUD_FILTER = 60;
    CLD_PRB_THRESH = 50;
    NIR_DRK_THRESH = 0.2;
    CLD_PRJ_DIST = 1;
    BUFFER = 50;

    #filter images based on cloudy percentage (metadata)
    S2SRCol = S2SRCol.filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER))

    #join S2SR with S2CloudCol
    joined = ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(
      primary = S2SRCol,
      secondary = S2CloudCol,
      condition = ee.Filter.equals(
          leftField = 'system:index',
          rightField = 'system:index'
      )))
  
    def add_cloud_bands(img):
        #Get s2cloudless image, subset the probability band.
        cld_prb = ee.Image(img.get('s2cloudless')).select('probability')
        #Condition s2cloudless by the probability threshold value
        is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')
        #Add the cloud probability layer and cloud mask as image bands
        return img.addBands(ee.Image([cld_prb, is_cloud]))
  
    def add_shadow_bands(img):
        #Identify water pixels from the SCL band
        not_water = img.select('SCL').neq(6)
        #Identify dark NIR pixels that are not water (potential cloud shadow pixels)
        SR_BAND_SCALE = 1e4
        dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')
        #Determine the direction to project cloud shadow from clouds (assumes UTM projection)
        shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')))
        #Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input
        cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
            .reproject(crs=img.select(0).projection(), scale=100)
            .select('distance')
            .mask()
            .rename('cloud_transform'))
        #Identify the intersection of dark pixels with cloud shadow projection
        shadows = cld_proj.multiply(dark_pixels).rename('shadows')
        #Add dark pixels, cloud projection, and identified shadows as image bands
        return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

    def add_cld_shdw_mask(img):
        #Add cloud component bands.
        img_cloud = add_cloud_bands(img)
        #Add cloud shadow component bands.
        img_cloud_shadow = add_shadow_bands(img_cloud)
        #Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
        is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)
        #Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input
        #20 m scale is for speed, and assumes clouds don't require 10 m precision
        is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
            .reproject(crs=img.select([0]).projection(), scale= 20)
            .rename('cloudmask'))
        #Add the final cloud-shadow mask to the image
        return img_cloud_shadow.addBands(is_cld_shdw)

    def apply_cld_shdw_mask(img):
        #Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
        not_cld_shdw = img.select('cloudmask').Not();
        #Subset reflectance bands and update their masks, return the result.
        return img.updateMask(not_cld_shdw)

    s2_sr = joined.map(add_cld_shdw_mask).map(apply_cld_shdw_mask)

    return s2_sr



