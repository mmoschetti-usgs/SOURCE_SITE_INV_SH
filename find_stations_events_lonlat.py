#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 4 12:05

@author: mmoschetti

Script to identify stations and events with distances from lon/lat coordinates

"""

import sys
import os
import shutil
import pandas as pd
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.geodetics import (calc_vincenty_inverse, degrees2kilometers, gps2dist_azimuth, kilometer2degrees, locations2degrees)
from libcomcat.utils import read_phases
import glob
import time
import numpy as np
from obspy.clients.fdsn.header import FDSNNoDataException, FDSNException
np.warnings.filterwarnings('ignore')  # To ignore numpy deprecation warning from obspy tapering

# formatted output for writing
def format2(value):
    return "%02.0f" % value

#
cwd = os.getcwd()

# Import event and station list
#eventFile = str(cwd+'/master_event_list_2009_2017.csv')
#stationFile = str(cwd+'/master_station_list.csv')
#parameterFile = str(cwd+'/params_lonlat_dist.csv')
#stationFile = str('master_station_list.csv')
eventFile = str('event_list_with_data.csv')
stationFile = str('stations_with_data.csv')
parameterFile = str('params_lonlat_dist.csv')

df_event = pd.read_csv(eventFile)
df_station = pd.read_csv(stationFile)
df_setup = pd.read_csv(parameterFile)

print('Reading event file:', eventFile)
print('Reading station file:', stationFile)
print('Reading parameter file:', parameterFile)

# get coordinates and maximum distances from parameter file
lonCenter = df_setup.lon[0]
latCenter = df_setup.lat[0]
distSta = df_setup.distStation[0]
distEv = df_setup.distEvent[0]
minMag = df_setup.minMag[0]
print('Search for stations and events within ',distSta,' and ',distEv,' km from coordinate (',lonCenter,',',latCenter,')')
print('Enforcing minimum magnitude of ',minMag)

#
#eventFileOut=('event_list_dist'+str(distEv)+'km.csv')
#stationFileOut=('station_list_dist'+str(distEv)+'km.csv')
eventFileOut=('event_list_subset.csv')
stationFileOut=('station_list_subset.csv')
print('Output files: ',eventFileOut,stationFileOut)

# Set up IRIS client
#fdsn_client = Client('https://service.iris.edu')
fdsn_client = Client('IRIS')

# do a bulk call for the need stations
networks = df_station.network.tolist()
stations = df_station.station.tolist()
locations = df_station.location.tolist()
channels = df_station.channel.tolist()

# Some location codes might need fixing in the station file
for idx, loc in enumerate(locations):
    if loc == '__':
        locations[idx] = '--'

# preallocate pick arrays
ev_time_file = []
p_pick_time = []
s_pick_time = []
window = []
channel_list = []

# selecting events within distance threshold
cntEv=0
f = open(eventFileOut,'w')
outStr=('year,month,day,hour,minute,second,magnitude,latitude,longitude,id\n')
f.write(outStr)
for ev_num, event_id in enumerate(df_event.id):
#    print(event_id + ' (' + str(ev_num + 1) + '/' + str(len(df_event.id)) + ')')
    lonEv=df_event.longitude[ev_num]
    latEv=df_event.latitude[ev_num]
    magEv=df_event.magnitude[ev_num]
#    calc_vincenty_inverse(lat1, lon1, lat2, lon2, a=6378137.0, f=0.0033528106647474805)[source]
    dist_m,az,baz=calc_vincenty_inverse(latEv, lonEv, latCenter, lonCenter, a=6378137.0, f=0.0033528106647474805)
    dist_km=dist_m/1000
    if (dist_km < distEv and magEv>=minMag):
        print('Event within distance threshold: ',dist_km,'(',lonEv,latEv,')')
#        outStr=('year,month,day,hour,minute,second,magnitude,latitude,longitude,id\n')
        outStr=(str(df_event.year[ev_num])+','+str(format2(df_event.month[ev_num]))+','+str(format2(df_event.day[ev_num]))+','+str(format2(df_event.hour[ev_num]))+','+str(format2(df_event.minute[ev_num]))+','+str(format2(df_event.second[ev_num]))+','+str(df_event.magnitude[ev_num])+','+str(df_event.latitude[ev_num])+','+str(df_event.longitude[ev_num])+','+str(df_event.id[ev_num])+'\n')
        f.write(outStr)
        cntEv +=1
# list number of events selected and close file
print('Number of events within distance threshold of ',distEv,'km =',cntEv)
f.close()

# reference start and end times from input to parameter file
t_start = UTCDateTime('2009'+'-'+'01'+'-'+'01'+'T00:00:00.000')
t_end= UTCDateTime('2017'+'-'+'12'+'-'+'31'+'T23:59:59.000')
cntSta=0
print('Station sorting...')
f = open(stationFileOut,'w')
outStr=('network,station,location,channel\n')
f.write(outStr)
for numSta, sta in enumerate(df_station.station):
    netW=df_station.network[numSta]
    staLoc=str(df_station.location[numSta])
    if staLoc == '__':
        staLoc='--'
    staChan=df_station.channel[numSta]
#    print('   ',numSta,sta,netW,staLoc,staChan)
    try:
        inv = fdsn_client.get_stations(network=netW, station=sta,
                                           location=staLoc, channel=staChan,
                                           starttime=t_start, endtime=t_end,
                                           level='channel')
    except FDSNNoDataException:
        print('Unable to get station inventory from IRIS for', sta)
#        print('Continuing to next station...')
        continue

    net = inv[0]
    staSt = net[0]
    station_lat = staSt.latitude
    station_lon = staSt.longitude
    station_ele = staSt.elevation
    dist_m,az,baz=calc_vincenty_inverse(station_lat, station_lon, latCenter, lonCenter, a=6378137.0, f=0.0033528106647474805)
    dist_km=dist_m/1000
    if (dist_km < distSta ):
        print('Station',sta,'within distance threshold: ',dist_km,'(',station_lon,station_lat,')')
        outStr=(netW+','+sta+','+staLoc+','+staChan+'\n')
        f.write(outStr)
        cntSta +=1
# list number of stations and close file
print('Number of stations within distance threshold of ',distSta,'km =',cntSta)
f.close()

