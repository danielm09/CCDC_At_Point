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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates



#get time series values - [id,lon,lat,time,band value]
def getTSvalues(geometry, img_collection, band_to_plot):
    import ee
    timeseries = ee.List(img_collection.select(band_to_plot).getRegion(geometry,scale=10)).getInfo() #change scale for landsat
    timeseries = np.array(timeseries)

    return timeseries



#create artificial dates for plotting the regression (plug values into regression equation)
def createArtificialDates(first_date, date_end):
    import ee
    #create sequence of dates from first date to date_end, spaced by 10 days
    interval = 10 #days
    
    #convert date end to milliseconds
    date_end_millis = ee.Date(date_end).millis().getInfo()
    
    #calculate the number of intervals
    num_intervals = int((date_end_millis-first_date)/(interval*24*60*60*1000))
    
    #determine the artificial dates
    artificial_dates = [first_date+x*interval*24*60*60*1000 for x in range(num_intervals)]

    #adjust end of series
    if artificial_dates[-1]<date_end_millis:
        artificial_dates.append(date_end_millis)
    elif artificial_dates[-1]>date_end_millis:
        a.pop(-1)
    artificial_dates = np.array(artificial_dates)

    return artificial_dates


def plotCCDCFitting(geometry, img_collection, band_to_plot, ccdc_result_info, artificial_dates):
    
    #get the number of segments (fits)
    nsegments = len(ccdc_result_info['tBreak'][0])
    
    #list to store predicted values of each segment. values of segment 0 are stored in the first element of the list and so on
    predicted_values = []
    fig, ax = plt.subplots(figsize=(10, 6))
    for seg in range(nsegments):
        artificial_dates_seg = artificial_dates[(artificial_dates<ccdc_result_info['tEnd'][0][seg])&(artificial_dates>ccdc_result_info['tStart'][0][seg])]
        coefs = ccdc_result_info['{}_coefs'.format(band_to_plot)][0][seg]
        #plug values into the regression equation
        pred = [coefs[0]+coefs[1]*t+
                coefs[2]*np.cos(t*1*2*np.pi/(365.25*24*60*60*1000))+
                coefs[3]*np.cos(t*1*2*np.pi/(365.25*24*60*60*1000))+
                coefs[4]*np.cos(t*2*2*np.pi/(365.25*24*60*60*1000))+
                coefs[5]*np.cos(t*2*2*np.pi/(365.25*24*60*60*1000))+
                coefs[6]*np.cos(t*3*2*np.pi/(365.25*24*60*60*1000))+
                coefs[7]*np.cos(t*3*2*np.pi/(365.25*24*60*60*1000))
                for t in artificial_dates_seg]
           
        predicted_values.append(pred)
        ax.plot(pd.to_datetime(artificial_dates_seg,unit='ms'),pred,label='Fit {}'.format(seg+1))

    #plot actual observations
    timeseries = getTSvalues(geometry, img_collection, band_to_plot)
    times = np.stack(timeseries,axis=1)[:][-2][1:].astype('int64')
    values = np.stack(timeseries,axis=1)[:][-1][1:].astype('float')

    ax.scatter(pd.to_datetime(times,unit='ms'),values,color='red',s=1,label='Observations')

    #adjust chart elements
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.xticks(rotation=90)
    plt.ylabel(band_to_plot.upper())
    plt.xlabel('Date')
    plt.grid(alpha=0.5)
    ax.legend()

    plt.tight_layout()
    plt.show()







    
