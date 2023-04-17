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



from ccdc_retriever.core import gee_data_sentinel as img_aq_s2
from ccdc_retriever.core import gee_data_landsat as img_aq_ls
from ccdc_retriever.core import plot_data
from ccdc_retriever.core import process_ccdc


def executeProcessing(name, coords, date_start, date_end, cloud_filter, bands, band_to_plot,
                      breakpointbands, tmask, numObs, chi, minYears, Lambda):
    import ee

    #creates earth engine point geometry based on coordinates
    geometry = ee.Geometry.Point(coords)

    #get image collection according to what was entered by the user
    if name=='Sentinel-2':
        img_col_filtered = img_aq_s2.getImageCollection(name, geometry, date_start, date_end, cloud_filter, bands)
    elif name=='Landsat col. 1':
        img_col_filtered = img_aq_ls.getImageCollection(img_aq_ls.get_full_collection(geometry,[date_start,date_end],None,1),geometry)
    elif name=='Landsat col. 2':
        img_col_filtered = img_aq_ls.getImageCollection(img_aq_ls.get_full_collection(geometry,[date_start,date_end],None,2),geometry)


    #get time series values at point location    
    timeseries = plot_data.getTSvalues(geometry,img_col_filtered,band_to_plot)
    #get first date of TS
    first_date = int(timeseries[1][3])
    
    #create artificial dates to plot the regression fit (plug-in dates into equations)
    artificial_dates = plot_data.createArtificialDates(first_date, date_end)

    #get ccdc result
    ccdc_result_info = process_ccdc.runCCDC(geometry, img_col_filtered, breakpointbands, tmask, numObs, chi, minYears, Lambda)

    #plot ccdc fits and time series observations
    plot_data.plotCCDCFitting(geometry, img_col_filtered, band_to_plot, ccdc_result_info, artificial_dates)
    
    
