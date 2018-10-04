#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:22:13 2018

@author: jrekoske

Prepares 'fortran_stations_subset.txt' using a list of stations, events, and a
data directory. Requires that the files 'station_list_subset.csv' and
'master_event_list_2009_2017.csv' are in the current directory. The data
directory must be called 'data' and must have folders for each event.

"""

import os
import pandas as pd
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from glob import glob

cwd = os.getcwd()

df_station = pd.read_csv('station_list_subset.csv')
df_event = pd.read_csv('master_event_list_2009_2017.csv')
fdsn_client = Client('IRIS')

networks = df_station.network.tolist()
stations = df_station.station.tolist()
locations = df_station.location.tolist()
channels = df_station.channel.tolist()

num_events = df_event.shape[0]
t_start = UTCDateTime(str(df_event.year[0])+'-'+str(df_event.month[0])+'-'+str(df_event.day[0])+'T'+str(df_event.hour[0])+':'+str(df_event.minute[0])+':'+str(df_event.second[0])+'.000')
t_end = UTCDateTime(str(df_event.year[num_events - 1])+'-'+str(df_event.month[num_events - 1])+'-'+str(df_event.day[num_events - 1])+'T'+str(df_event.hour[num_events - 1])+':'+str(df_event.minute[num_events - 1])+':'+str(df_event.second[num_events - 1])+'.000')

# Create the file for stations to later append
f = open('fortran_stations_subset.txt', 'w')
f.close()

for idx, stn in enumerate(stations):

    # Check if we have any sac files for this station in any of the event folders
    have_sac_files = False
    event_folders = glob('data/*/')

    for event_fold in event_folders:
        matching_sac_files = glob(event_fold + '/*' + stn + '*.sac')
        if matching_sac_files:
            have_sac_files = True
            break

    if have_sac_files:

        # Get station coordinates and azimuth information from IRIS
        netst = networks[idx]
        stachan = channels[idx]
        staloc = str(locations[idx])

        # Fix the location if its length is one
        if len(staloc) is 1:
            staloc = '0' + staloc

        print('inv-get_stations: station/network: ', stn, netst)
        print('staloc/start-time/end-time: ', staloc, t_start, t_end)

        try:
            inv = fdsn_client.get_stations(network=netst, station=stn,
                                           location=staloc, channel=stachan,
                                           starttime=t_start, endtime=t_end,
                                           level='channel')
        except FDSNNoDataException:
            print('Unable to get station inventory from IRIS for', stn)

            # Waveform data is useless without having station information
            # Loop through all directories and delete those SAC files as well
            files_in_data = os.listdir(cwd + '/data')
            for fil in files_in_data:
                if os.path.isdir(cwd + '/data/' + fil):
                    files_to_delete = glob(cwd + '/data/' + fil + '/*' + stn + '*.sac')
                    for del_fil in files_to_delete:
                        os.remove(del_fil)

            print('Continuing to next station.')
            continue

        net = inv[0]
        sta = net[0]
        station_lat = sta.latitude
        station_lon = sta.longitude
        station_ele = sta.elevation

        # East chan think i need to add 90 to dip, this has 0 as horizontal
        chan = sta[0]
        station_az_e = chan.azimuth
        station_dip_e = chan.dip + 90

        # north channel, add 90 to dip
        chan = sta[1]
        station_az_n = chan.azimuth
        station_dip_n = chan.dip + 90

        # think i need to add 90 to dip, this has 0 as horizontal, Z comp add 90 to dip
        chan = sta[2]
        station_az_z = chan.azimuth
        station_dip_z = chan.dip + 90

        while (len(stn) < 5):
            stn += '_'

        # format and append as list
        stabulk = 'STNID ' + stn + '\n'
        stabulk2 = 'COORD ' + str(station_lat) + ' ' + str(station_lon) + ' 0\n'
        stabulk3 = 'CHANALT 1 ' + str(station_dip_z) + '\n'
        stabulk4 = 'CHANAZI 1 ' + str(station_az_z) + '\n'
        stabulk5 = 'CHANALT 2 ' + str(station_dip_e) + '\n'
        stabulk6 = 'CHANAZI 2 ' + str(station_az_e) + '\n'
        stabulk7 = 'CHANALT 3 ' + str(station_dip_n) + '\n'
        stabulk8 = 'CHANAZI 3 ' + str(station_az_n) + '\n'

        with open('fortran_stations_subset.txt', 'a') as myfile:
            myfile.write(stabulk)
            myfile.write(stabulk2)
            myfile.write(stabulk3)
            myfile.write(stabulk4)
            myfile.write(stabulk5)
            myfile.write(stabulk6)
            myfile.write(stabulk7)
            myfile.write(stabulk8)

        myfile.close()
