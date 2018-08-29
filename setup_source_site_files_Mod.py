#!/usr/bin/env python

# User needs to run "conda activate comcat" before running script

## Import needed packages and functions
import obspy as op
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from os import listdir
from os.path import isfile, join
import fnmatch
import libcomcat
#from libcomcat.utils import phase_reader
from libcomcat.utils import read_phases
from datetime import date
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# Get the working directory so everything is a subset of that
cwd = os.getcwd()
# Import event and station list
df_event=pd.read_csv(cwd+'/event_list.csv')
df_station=pd.read_csv(cwd+'/station_list.csv')
df_setup=pd.read_csv(cwd+'/setup_parameters.csv')

# Define filter parameters from setup file
fmin1=df_setup.fmin1[0]
fmin2=df_setup.fmin2[0]
fmax1=df_setup.fmax1[0]
fmax2=df_setup.fmax2[0]
pre_filt = (fmin1, fmin2, fmax1, fmax2)

# example script exiting
#sys.exit()

# Other set up parameters
waveformlength=df_setup.waveformlength[0]
magfilt=df_setup.magfilt[0]
magfiltlim=df_setup.magfiltlowerlimit[0]
pick_type=df_setup.PorS[0]
minWindow=df_setup.windowlength[0]
outputRate=df_setup.outputSampleRate[0]
# REMOVE
# Call for sac files from IRIS
fdsn_client = Client('IRIS')
# some hard-coded parameters
# minimum window duration for spectra
#minWindow=10
print("Using minimum window duration of ",minWindow," s")
#print("Modify minWindow variable, if needed...")
print("Decimating to output rate of ",outputRate," Hz")

# do a bulk call for the need stations
networks = df_station.network.tolist()
stations = df_station.station.tolist()
locations = df_station.location.tolist()
channels = df_station.channel.tolist()

# preallocate pick arrays
ev_time_file=[]
p_pick_time=[]
s_pick_time=[]
window=[]
stat_time=[]

# Loop over events
for ev_num in range(len(df_event.id)):
    event_id = df_event.id[ev_num]
    t1 = UTCDateTime(str(df_event.year[ev_num])+"-"+str(df_event.month[ev_num])+"-"+str(df_event.day[ev_num])+"T"+str(df_event.hour[ev_num])+":"+str(df_event.minute[ev_num])+":"+str(df_event.second[ev_num])+".000")
    t2 = t1+waveformlength

    # creat bulk variable to call all the stations for this event
    bulk=[]
    for x in range(len(networks)):
        bulkrow=(networks[x],stations[x],locations[x],channels[x],t1,t2)
        bulk.append(bulkrow)

    # Load picks for this event
    picks_folder = cwd+'/Phase_Data/'
    pick_id = df_event.id[ev_num]
    pick_id=pick_id[0:11]
    pick_id_file = picks_folder + pick_id + "_phases.csv"
    #events,picks=phase_reader(pick_id_file)
    events,picks=read_phases(pick_id_file)
    p_indices = [f for f, s in enumerate(picks.Phase) if 'P' in s]
    s_indices = [f for f, s in enumerate(picks.Phase) if 'S' in s]

    # Get all the waveform for this event from IRIS
    st = fdsn_client.get_waveforms_bulk(bulk,attach_response=True)
    # Loop over stations
    for st_num in range(len(st)):
        # Need to save the picks in the right order
        stast=st[st_num].stats.station
        indices = [f for f, s in enumerate(picks.Channel) if stast in s]
        use_ind=set(p_indices) & set(indices)
        use_ind_s=set(s_indices) & set(indices)
        # Skip stations that do not have picks
        if (not use_ind or not use_ind_s):
            print("event "+event_id+" station "+stast+" has no p or s picks")
        else:
            # Remove response from waveform, output is velocity, detrend and taper
            #print("working on station "+str(st_num)+stast)
#            st[st_num].remove_response(output='VEL', pre_filt=pre_filt, zero_mean=True).detrend('simple').taper(0.05,type = 'hann')
            tr=st[st_num]
            tr.detrend('simple').taper(0.05,type = 'hann').remove_response(output='VEL', pre_filt=pre_filt, zero_mean=True)

# decimate sample rate to outputRate
            if int(tr.stats.sampling_rate) > outputRate :
                decimateVal=int(tr.stats.sampling_rate/outputRate)
                print("Decimating to output rate of ",outputRate," Hz, from ",tr.stats.sampling_rate," Hz (factor of ",decimateVal,")")
            else:
                print("Output sample rate below record sample rate (",outputSampleRate,",",tr.stats.sampling_rate,")")
                sys.exit()
            tr.decimate(decimateVal)

            # Get picks in ordinal time and then make them relative to the start time
            pick_time=picks['Arrival Time'][use_ind]
            lpt=str(list(pick_time))
            pick_ord=date.toordinal(date(int(lpt[2:6]),int(lpt[7:9]),int(lpt[10:12])))+int(lpt[13:15])/24+int(lpt[16:18])/1440+int(lpt[19:21])/86400

            pick_time_s=picks['Arrival Time'][use_ind_s]
            lpts=str(list(pick_time_s))
            pick_ord_s=date.toordinal(date(int(lpts[2:6]),int(lpts[7:9]),int(lpt[10:12])))+int(lpts[13:15])/24+int(lpts[16:18])/1440+int(lpts[19:21])/86400

            evtime=str(tr.stats.starttime)
            evtime_date=date.toordinal(date(int(evtime[0:4]),int(evtime[5:7]),int(evtime[8:10])))+int(evtime[11:13])/24+int(evtime[14:16])/1440+float(evtime[17:26])/86400

            # convert picks to relative time to the event start in seonds
            pick_file_diff=(pick_ord-evtime_date)*86400
            pick_file_diff_s=(pick_ord_s-evtime_date)*86400
            pick_file_diff_ps=(pick_ord_s-pick_ord)*86400

            # filter by event magnitude
            if (magfilt=='yes'):
                mag=df_event.magnitude[ev_num]
                fmin=10**(mag*-1/2.3+1)-(0.496-magfiltlim)
            else:
                fmin=fmin2
            if (fmin<fmin2):
                fmin=fmin2
            else:
                fmin=fmin

            # Copy original and then apply detrend and highpass filter
            tr_filt= tr.copy()
            tr_filt.detrend('simple').taper(0.05,type = 'hann').filter('highpass', freq=fmin, corners=2, zerophase=True)

            # need to get all output names to be the same length
            stm=str(df_event.month[ev_num])
            if (len(stm)==1):
                stm='0'+stm
            else:
                stm=stm
            stday=str(df_event.day[ev_num])
            if (len(stday)==1):
                stday='0'+stday
            else:
                stday=stday

            sthour=str(df_event.hour[ev_num])
            if (len(sthour)==1):
                sthour='0'+sthour
            else:
                sthour=sthour
            stminute=str(df_event.minute[ev_num])
            if (len(stminute)==1):
                stminute='0'+stminute
            else:
                stminute=stminute
            stsecond=str(df_event.second[ev_num])
            if (len(stsecond)==1):
                stsecond='0'+stsecond
            else:
                stsecond=stsecond

            stnn=tr.stats.station
            if (len(stnn)==4):
                stnn=stnn+"_"
            else:
                stnn=stnn
            # Create a directory for each event if it doesn't already exist
            if not os.path.exists(cwd+"/data/"):
                os.makedirs(cwd+"/data/")
            if not os.path.exists(cwd+"/data/"+str(df_event.year[ev_num])+stm+stday+sthour+stminute+stsecond):
                os.makedirs(cwd+"/data/"+str(df_event.year[ev_num])+stm+stday+sthour+stminute+stsecond)
            # Write out the sac file to the main data directory and the event folder
            out_file=cwd+"/data/"+str(df_event.year[ev_num])+stm+stday+sthour+stminute+stsecond+"/"+str(df_event.year[ev_num])+stm+stday+sthour+stminute+stsecond+"_"+stnn+"_"+tr.stats.channel+".sac"
            tr_filt.write(out_file, format='SAC')
            out_file=cwd+"/data/"+str(df_event.year[ev_num])+stm+stday+sthour+stminute+stsecond+"_"+stnn+"_"+tr.stats.channel+".sac"
            tr_filt.write(out_file, format='SAC')
            # Save the picks and event and station ids to use below
            ccc=tr.stats.channel
            if (ccc[2]=='Z'):
                p_pick_time.append(str(pick_file_diff))
                s_pick_time.append(str(pick_file_diff_s))
                window.append(str(pick_file_diff_ps))
                ev_time_file.append(str(df_event.year[ev_num])+stm+stday+sthour+stminute+stsecond)
                stat_time.append(stnn)
            else:
                print("outputing sac files")
    print("event "+str(ev_num)+" is done")


## Create the input files to the inversion code (fortran) in the correct format
# first files to create is the event list
# gather all the output files so that we see if any were skipped due to lack of picks or waveforms
t_start = UTCDateTime(str(df_event.year[0])+"-"+str(df_event.month[0])+"-"+str(df_event.day[0])+"T"+str(df_event.hour[0])+":"+str(df_event.minute[0])+":"+str(df_event.second[0])+".000")
t_end = UTCDateTime(str(df_event.year[ev_num])+"-"+str(df_event.month[ev_num])+"-"+str(df_event.day[ev_num])+"T"+str(df_event.hour[ev_num])+":"+str(df_event.minute[ev_num])+":"+str(df_event.second[ev_num])+".000")

# Create the file for stations to later append
f = open('fortran_stations.txt','w')
f.close()
sacdir=cwd+'/data/'
sacfiles=listdir(sacdir)
sacfiles=[f for f  in sacfiles if '.sac' in f]
# Loop over stations and append text file in the correct format for fortran
for xx in range(len(stations)):
    sta_xx=stations[xx]
    stafiles=[f for f  in sacfiles if sta_xx in f]
    # make sure there are sac files for that station
    if (len(stafiles)>0):
        # Get station coordinates and azimuth information from IRIS
        netst=networks[xx]
        stast=stations[xx]
        stachan=channels[xx]
        staloc=locations[xx]
        inv=fdsn_client.get_stations(network=netst,station=stast,location=staloc,channel=stachan,starttime=t_start,endtime=t_end,level='channel')
        net=inv[0]
        sta=net[0]
        station_lat=sta.latitude
        station_lon=sta.longitude
        station_ele=sta.elevation
        # East chan think i need to add 90 to dip, this has 0 as horizontal
        chan=sta[0]
        station_az_e=chan.azimuth
        station_dip_e=chan.dip+90
        #north channel, add 90 to dip
        chan=sta[1]
        station_az_n=chan.azimuth
        station_dip_n=chan.dip+90
        # think i need to add 90 to dip, this has 0 as horizontal, Z comp add 90 to dip
        chan=sta[2]
        station_az_z=chan.azimuth
        station_dip_z=chan.dip+90
        if (len(stast)==4):
            stast=stast+"_"
        else:
            stast=stast

        # format and append as list
        stabulk="STNID "+stast+"\n"
        stabulk2="COORD "+str(station_lat)+" "+str(station_lon)+" 0\n"
        stabulk3="CHANALT 1 "+str(station_dip_z)+"\n"
        stabulk4="CHANAZI 1 "+str(station_az_z)+"\n"
        stabulk5="CHANALT 2 "+str(station_dip_e)+"\n"
        stabulk6="CHANAZI 2 "+str(station_az_e)+"\n"
        stabulk7="CHANALT 3 "+str(station_dip_n)+"\n"
        stabulk8="CHANAZI 3 "+str(station_az_n)+"\n"

        with open("fortran_stations.txt", "a") as myfile:
            myfile.write(stabulk)
            myfile.write(stabulk2)
            myfile.write(stabulk3)
            myfile.write(stabulk4)
            myfile.write(stabulk5)
            myfile.write(stabulk6)
            myfile.write(stabulk7)
            myfile.write(stabulk8)

        myfile.close()
    else:
        print("no waveforms or picks for station "+stast)

# Create the files of events, with phase arrival and window duration
#!/usr/bin/env python

#

## Import needed packages and functions
dataDir=cwd+"/data"
#eventDir=[]
allDataDir=listdir(dataDir)
cnt=0
for file1 in listdir(dataDir):
    if fnmatch.fnmatch(file1, '??????????????'):
        evDir=cwd+"/data/"+file1
        evDirNm=file1
#        print(evDir)
        outputName=cwd+"/data/"+evDirNm+"NAMESa.DAT"
        with open(outputName, "w") as myfile:
#        f = open(outputName,'w')
#    f.close()
            print(evDirNm, outputName)
            sacfilesEv=listdir(evDir)
            sacfilesEv=[f for f  in sacfilesEv if '.sac' in f]
            sacfiles_z=[f for f  in sacfilesEv if 'Z.sac' in f]
        
            for xx in range(len(sacfiles_z)):
                file_name=sacfiles_z[xx]
                sacfiles_all=[f for f  in sacfilesEv if file_name[0:20] in f]
                file_name_e= [f for f  in sacfiles_all if 'E.sac' in f or '1.sac' in f]
                file_name_n= [f for f  in sacfiles_all if 'N.sac' in f or '2.sac' in f]
                time_name=file_name[0:14]
                stast=file_name[15:20]
                ind1=[f for f, s in enumerate(stat_time) if stast in s]
                ind2=[f for f, s in enumerate(ev_time_file) if time_name in s]
                ind_all=set(ind1) & set(ind2)
                ind_all=list(map(int,ind_all))
                ind_all=[int(i) for i in ind_all]
                # match the picks and use the input from the setup parameters file to decide if you want P or S
                ppick=[p_pick_time[i] for i in ind_all]
                spick=[s_pick_time[i] for i in ind_all]
                win=[window[i] for i in ind_all]
                if (pick_type=='P'):
                    pick_use=ppick
                else:
                    pick_use=spick
            
                file_name_e=''.join(file_name_e)
                file_name_n=''.join(file_name_n)
                pick_use=''.join(pick_use)
                win=''.join(win)
                win=win[0:5]
                winFloat=float(win)
                if winFloat<minWindow:
                    win=str(minWindow)
                pick_use=pick_use[0:5]
                line_z=file_name+" "+pick_use+" "+win+"\n"
                line_e=file_name_e+"\n"
                line_n=file_name_n+"\n"
            
                myfile.write(line_z)
                myfile.write(line_e)
                myfile.write(line_n)
                print("Completed event directory: ",evDir)
        myfile.close()
    else:
        print("Not a directory, ",file1)

            
# clean up data directory
sacdir=cwd+'/data/'
sacfiles=listdir(sacdir)
sacfiles=[f for f  in sacfiles if '.sac' in f]
for filef in sacfiles:
	print(filef)
	filefFull=cwd+'/data/'+filef
	os.remove(filefFull)


