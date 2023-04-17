"""
/***************************************************************************
 ChangeDetectionDialog
                                 A QGIS plugin
 Retrives CCDC's information and plot chart for a given point.
 
                             -------------------
        begin                : 2023-04-04
        git sha              : $Format:%H$
        copyright            : (C) 2023 by Daniel Moraes
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


#run ccdc
def runCCDC(geometry, img_collection, breakpointbands, tmask, numObs, chi, minYears, Lambda):
    import ee
    #run algorithm
    ccdc_result = ee.Algorithms.TemporalSegmentation.Ccdc(img_collection.select(breakpointbands), breakpointbands, tmask, numObs, chi, minYears, 2, Lambda)
    #apply reduce region to get data at point
    ccdc_result_info = ccdc_result.reduceRegion(ee.Reducer.toList(),
                                                geometry,
                                                scale=10,).getInfo() #change if using landsat
                                                
                                                
                                                
    return ccdc_result_info



