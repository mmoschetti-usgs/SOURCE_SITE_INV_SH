#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 09:37:24 2018
Updated Thu Oct 4 12:05

@author: jrekoske

Faster version of setup source site files.

"""

import sys
import os
import shutil
import pandas as pd
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from libcomcat.utils import read_phases
import glob
import time
import numpy as np
from obspy.clients.fdsn.header import FDSNNoDataException, FDSNException
np.warnings.filterwarnings('ignore')  # To ignore numpy deprecation warning from obspy tapering


def channel_split(chan):
    channel_split = chan.split('.')
    network = channel_split[0]
    station = channel_split[1]
    channel = channel_split[2]
    location = channel_split[3]

    if len(location) == 1:
        location = '0' + location

    return (network, station, channel, location)


cwd = os.getcwd()

# directory clean-up
if os.path.isdir(cwd + '/data/'):
    shutil.rmtree(cwd + '/data/')
if os.path.exists('fortran_stations.txt'):
    os.remove('fortran_stations.txt')

# Import event and station list
eventFile = str(cwd+'/master_event_list_2009_2017.csv')
stationFile = str(cwd+'/master_station_list.csv')
parameterFile = str(cwd+'/setup_parameters.csv')

df_event = pd.read_csv(eventFile)
df_station = pd.read_csv(stationFile)
df_setup = pd.read_csv(parameterFile)

print('Reading event file:', eventFile)
print('Reading station file:', stationFile)
print('Reading parameter file:', parameterFile)

# Define filter parameters from setup file
fmin1 = df_setup.fmin1[0]
fmin2 = df_setup.fmin2[0]
fmax1 = df_setup.fmax1[0]
fmax2 = df_setup.fmax2[0]
pre_filt = (fmin1, fmin2, fmax1, fmax2)

# Other set up parameters
waveformlength = df_setup.waveformlength[0]
magfilt = df_setup.magfilt[0]
magfiltlim = df_setup.magfiltlowerlimit[0]
pick_type = df_setup.PorS[0]
minWindow = df_setup.windowlength[0]
outputRate = df_setup.outputSampleRate[0]

# Set up IRIS client
fdsn_client = Client('https://service.iris.edu')
print('Using minimum window duration of', minWindow, 's')
print('Decimating to output rate of', outputRate, 'Hz')

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

for ev_num, event_id in enumerate(df_event.id):

    print(event_id + ' (' + str(ev_num + 1) + '/' + str(len(df_event.id)) + ')')
    t1 = UTCDateTime(str(df_event.year[ev_num])+'-'+str(df_event.month[ev_num])+'-'+str(df_event.day[ev_num])+'T'+str(df_event.hour[ev_num])+':'+str(df_event.minute[ev_num])+':'+str(df_event.second[ev_num])+'.000')
    t2 = t1 + waveformlength

    # Load picks for this event
    pick_id = df_event.id[ev_num]
    phase_file = cwd + '/Phase_Data/' + pick_id + '_phases.csv'

    if not os.path.exists(phase_file):
        print('The phase file', phase_file, 'does not exist. Continuing to next event.')
        print('------------------------------------------------------------------')
        continue

    try:
        events, picks = read_phases(phase_file)
    except pd.errors.EmptyDataError:
        print('Phase file did not have any phase data. Continuing on to the next event.')
        print('------------------------------------------------------------------')
        continue

    bulk = []
    for chan_idx1, chan1 in enumerate(picks.Channel):

        network1, station1, channel1, location1 = channel_split(chan1)

        # Check if we have an S-pick
        if 'S' in picks.Phase[chan_idx1]:

            # Loop through the master station list, and see if there's a match
            for chan_idx2, chan2 in enumerate(channels):
                network2 = networks[chan_idx2]
                station2 = stations[chan_idx2]
                location2 = locations[chan_idx2]

                if network1 == network2 and station1 == station2 and location1 == location2 and channel1[0:1] == chan2[0:1]:
                    bulkrow = (network2, station2, location2, chan2, t1, t2)
                    if bulkrow not in bulk:
                        bulk.append(bulkrow)

    if (len(bulk) == 0):
        print('No S picks are available. Continuing on to the next event.')
        continue

    print(bulk)

    # Get all the waveforms for this event from IRIS
    success = False
    num_tries = 0
    while (success is False and num_tries < 5):
        try:
            start = time.time()
            st = fdsn_client.get_waveforms_bulk(bulk, attach_response=True)
            end = time.time()
            print('Downloaded', len(st), 'waveforms in', end - start, 'seconds.')
            success = True
        except FDSNNoDataException:
            print('No data was available.')
            break
        except FDSNException:
            print('FDSN Error trying to request data. Trying again...')
            num_tries += 1
        except KeyboardInterrupt:
            sys.exit('Keyboard interrupt.')
        except:
            print('Error trying to request data. Trying again...')
            num_tries += 1

    if not success:
        print('Failed to get waveform data. Continuing on to the next event.')
        print('------------------------------------------------------------------')
        continue

    for tr in st:

        # Detrend and taper
        tr.detrend('simple').taper(0.05, type='hann')

        # Try to remove response
        # If if response information is not available, then just continue to
        # the next trace wihtout saving trace data
        try:
            tr.remove_response(output='VEL', pre_filt=pre_filt, zero_mean=True)
        except:
            print('Unable to remove instrument response for:', tr.get_id())
            continue

        # Decimate sampling rate to outputRate
        if int(tr.stats.sampling_rate) > outputRate:
            decimateVal = int(tr.stats.sampling_rate / outputRate)
            print('Decimating to output rate of', outputRate, 'Hz, from ', tr.stats.sampling_rate, 'Hz (factor of', str(decimateVal) + ')')
            tr.decimate(decimateVal)
        else:
            print('Output sample rate below record sample rate:', outputRate, tr.stats.sampling_rate, ')')
            print('Not writing event', event_id, 'station', tr.stats.station, 'to file.')
            continue

        # Get the P and S pick times for this trace
        found_p_pick = False

        for chan_idx, chan in enumerate(picks.Channel):

            pick_network, pick_station, pick_location, pick_channel = channel_split(chan)

            if tr.stats.network == pick_network and tr.stats.station == pick_station:

                if 'P' in picks.Phase[chan_idx]:
                    pick_time_p = UTCDateTime(picks['Arrival Time'][chan_idx])
                    print('The P arrival for', tr.stats.station, 'is', pick_time_p)
                    found_p_pick = True

                if 'S' in picks.Phase[chan_idx]:
                    pick_time_s = UTCDateTime(picks['Arrival Time'][chan_idx])
                    print('The S arrival for', tr.stats.station, 'is', pick_time_s)

        tr_start = tr.stats.starttime
        pick_file_diff_s = pick_time_s - tr_start

        if found_p_pick:
            pick_file_diff_p = pick_time_p - tr_start
            pick_file_diff_ps = pick_time_s - pick_time_p
        else:
            pick_file_diff_p = ''
            pick_file_diff_ps = minWindow

        # Filter by event magnitude
        if (magfilt == 'yes'):
            mag = df_event.magnitude[ev_num]
            fmin = 10**(mag*-1/2.3+1)-(0.496-magfiltlim)
        else:
            fmin = fmin2
        if (fmin < fmin2):
            fmin = fmin2
        else:
            fmin = fmin

        tr.detrend('simple').taper(0.05, type='hann').filter('highpass', freq=fmin, corners=2, zerophase=True)

        # Format the time and station to have consistent lengths
        year = str(t1.year)
        month = str(t1.month)
        day = str(t1.day)
        hour = str(t1.hour)
        minute = str(t1.minute)
        second = str(t1.second)

        if (len(month) == 1):
            month = '0' + month

        if (len(day) == 1):
            day = '0' + day

        if (len(hour) == 1):
            hour = '0' + hour

        if (len(minute) == 1):
            minute = '0' + minute

        if (len(second) == 1):
            second = '0' + second

        stnn = tr.stats.station
        while (len(stnn) < 5):
            stnn += '_'

        ev_string = year + month + day + hour + minute + second

        # Create a directory for each event if it doesn't already exist
        if not os.path.exists(cwd + '/data/'):
            os.makedirs(cwd + '/data/')
        if not os.path.exists(cwd + '/data/' + ev_string):
            os.makedirs(cwd+'/data/' + ev_string)

        # Write out the sac file to the main data directory and the event folder
        out_file = cwd + '/data/' + ev_string + '/' + ev_string + '_' + stnn + '_' + tr.stats.channel + '.sac'
        tr.write(out_file, format='SAC')
        out_file = cwd + '/data/' + ev_string + '_' + stnn + '_' + tr.stats.channel + '.sac'
        tr.write(out_file, format='SAC')

        # Save the picks and event and station ids to use below
        ccc = tr.stats.channel
        if (ccc[2] == 'Z'):
            p_pick_time.append(str(pick_file_diff_p))
            s_pick_time.append(str(pick_file_diff_s))
            window.append(str(pick_file_diff_ps))
            ev_time_file.append(ev_string)
            channel_list.append(stnn + '.' + tr.stats.location + '.' + tr.stats.channel[0:2])

    print('------------------------------------------------------------------')


# Create the input files to the inversion code (fortran) in the correct format
# first files to create is the event list
# gather all the output files so that we see if any were skipped due to lack of picks or waveforms
print('Creating Fortran stations file.')
t_start = UTCDateTime(str(df_event.year[0])+'-'+str(df_event.month[0])+'-'+str(df_event.day[0])+'T'+str(df_event.hour[0])+':'+str(df_event.minute[0])+':'+str(df_event.second[0])+'.000')
t_end = UTCDateTime(str(df_event.year[ev_num])+'-'+str(df_event.month[ev_num])+'-'+str(df_event.day[ev_num])+'T'+str(df_event.hour[ev_num])+':'+str(df_event.minute[ev_num])+':'+str(df_event.second[ev_num])+'.000')

# Create the file for stations to later append
f = open('fortran_stations.txt', 'w')
f.close()
sacdir = cwd + '/data/'
sacfiles = os.listdir(sacdir)
sacfiles = [f for f in sacfiles if '.sac' in f]

# Loop over stations and append text file in the correct format for fortran
for xx in range(len(stations)):
    sta_xx = stations[xx]
    stafiles = [f for f in sacfiles if sta_xx in f]
    # make sure there are sac files for that station
    if (len(stafiles) > 0):
        # Get station coordinates and azimuth information from IRIS
        netst = networks[xx]
        stast = stations[xx]
        stachan = channels[xx]
        staloc = str(locations[xx])
        if len(staloc) == 1:
            staloc = '0' + staloc

        print('inv-get_stations: station/network: ', stast, netst)
        print('staloc/start-time/end-time: ', staloc, t_start, t_end)

        # Sometimes station inventory is not available

        try:
            inv = fdsn_client.get_stations(network=netst, station=stast,
                                           location=staloc, channel=stachan,
                                           starttime=t_start, endtime=t_end,
                                           level='channel')
        except FDSNNoDataException:
            print('Unable to get station inventory from IRIS for', stast)

            # Waveform data is useless without having station information
            files_to_delete = glob.glob(cwd + '/data/' + '*' + stast + '*.sac')
            for fil in files_to_delete:
                os.remove(fil)

            # Loop through all directories and delete those SAC files as well
            files_in_data = os.listdir(cwd + '/data')
            for fil in files_in_data:
                if os.path.isdir(cwd + '/data/' + fil):
                    files_to_delete = glob.glob(cwd + '/data/' + fil + '/*' + stast + '*.sac')
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
        station_dip_e = chan.dip+90

        # north channel, add 90 to dip
        chan = sta[1]
        station_az_n = chan.azimuth
        station_dip_n = chan.dip+90

        # think i need to add 90 to dip, this has 0 as horizontal, Z comp add 90 to dip
        chan = sta[2]
        station_az_z = chan.azimuth
        station_dip_z = chan.dip+90

        while (len(stast) < 5):
            stast += '_'

        # format and append as list
        stabulk = 'STNID '+stast+'\n'
        stabulk2 = 'COORD '+str(station_lat)+' '+str(station_lon)+' 0\n'
        stabulk3 = 'CHANALT 1 '+str(station_dip_z)+'\n'
        stabulk4 = 'CHANAZI 1 '+str(station_az_z)+'\n'
        stabulk5 = 'CHANALT 2 '+str(station_dip_e)+'\n'
        stabulk6 = 'CHANAZI 2 '+str(station_az_e)+'\n'
        stabulk7 = 'CHANALT 3 '+str(station_dip_n)+'\n'
        stabulk8 = 'CHANAZI 3 '+str(station_az_n)+'\n'

        with open('fortran_stations.txt', 'a') as myfile:
            myfile.write(stabulk)
            myfile.write(stabulk2)
            myfile.write(stabulk3)
            myfile.write(stabulk4)
            myfile.write(stabulk5)
            myfile.write(stabulk6)
            myfile.write(stabulk7)
            myfile.write(stabulk8)

        myfile.close()

# Create the files of events, with phase arrival and window duration
print('------------------------------------------------------------------')
print('Creating DAT files.')

dataDir = cwd + '/data'
allDataDir = os.listdir(dataDir)
cnt = 0
for file1 in os.listdir(dataDir):
    if os.path.isdir(dataDir + '/' + file1):
        evDir = cwd+'/data/' + file1
        evDirNm = file1
        outputName = cwd + '/data/' + evDirNm + 'NAMESa.DAT'
        print('NAMESa File: ', outputName)

        with open(outputName, 'w') as myfile:

            print(evDirNm, outputName)
            sacfiles = os.listdir(evDir)
            sacfilesEv = os.listdir(evDir)
            sacfilesEv = [f for f in sacfilesEv if '.sac' in f]
            sacfiles_z = [f for f in sacfilesEv if 'Z.sac' in f]

            for xx in range(len(sacfiles_z)):

                file_name = sacfiles_z[xx]
                sacfiles_all = [f for f in sacfilesEv if file_name[0:23] in f]
                file_name_e = [f for f in sacfiles_all if 'E.sac' in f or '1.sac' in f]
                file_name_n = [f for f in sacfiles_all if 'N.sac' in f or '2.sac' in f]
                time_name = file_name[0:14]
                stast = file_name[15:20]

                ind1 = [f for f, s in enumerate(channel_list) if stast in s]
                ind2 = [f for f, s in enumerate(ev_time_file) if time_name in s]
                ind_all = set(ind1) & set(ind2)
                ind_all = list(map(int, ind_all))
                ind_all = [int(i) for i in ind_all]

                # match the picks and use the input from the setup parameters file to decide if you want P or S
                ppick = [p_pick_time[i] for i in ind_all]
                spick = [s_pick_time[i] for i in ind_all]
                win = [window[i] for i in ind_all]

                if (pick_type == 'P'):
                    pick_use = ppick
                else:
                    pick_use = spick

                file_name_e = ''.join(file_name_e)
                file_name_n = ''.join(file_name_n)
                pick_use = ''.join(pick_use)
                win = win[0]

                if len(win) >= 1:
                    print('filename, win: ', file_name, win)
                    winFloat = float(win)
                    print('winFloat: ', winFloat)
                    if winFloat < minWindow:
                        win = str(minWindow)
                    pick_use = pick_use[0:5]
                    line_z = file_name+' '+pick_use+' '+win+'\n'
                    line_e = file_name_e+'\n'
                    line_n = file_name_n+'\n'
                    myfile.write(line_z)
                    myfile.write(line_e)
                    myfile.write(line_n)

        print('Completed event directory: ', evDir)
        myfile.close()

# Clean up data directory
print('Cleaning up data directory.')
sacdir = cwd + '/data/'
sacfiles = os.listdir(sacdir)
sacfiles = [f for f in sacfiles if '.sac' in f]
for filef in sacfiles:
    filefFull = cwd + '/data/' + filef
    os.remove(filefFull)
