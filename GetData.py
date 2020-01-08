from obspy.clients.fdsn import Client
from obspy.core.event import Origin
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometer2degrees
from obspy.taup import TauPyModel
from obspy.io.sac.sactrace import SACTrace
from obspy import read_inventory, UTCDateTime

client = Client("IRIS")
model = TauPyModel(model='ak135')

import sys
import pdb

stn = 'YKW4'
ntwk = 'CN'
dir = '/Users/tianze/Research/WorkingDirectory/Yellowknife/Tianze/Downloads/YKW4'

pre_filt = None

def DeconvDecimate(st,f_prefilt):
	for tr in st:
		tr.detrend(type='constant')
		tr.detrend(type='linear')
		tr.remove_response(output='VEL',taper = True, taper_fraction=0.05, pre_filt= f_prefilt, water_level = 60.0)
	return st

# date = sys.argv[1]
# net_sta = sys.argv[2]
starttime = '1989-01-01';
endtime = '2016-01-01';

inv = client.get_stations(network='CN', station='YKW3',level='response')
info = inv.get_coordinates('CN.YKW3..BHZ','2001-01-01')
stla = info['latitude']
stlo = info['longitude']

cat = client.get_events(starttime=starttime, endtime=endtime, minmagnitude=5.5, catalog="ISC")

Event = cat.events
Event.reverse()

for evt in Event:
	org = evt.origins
	org = org[0]
	otime = org.time
	evlo = org.longitude
	evla = org.latitude
	evdp = org.depth

	if type(evdp) == float:
		evdp = evdp/1000
	else:
		continue

	geo = gps2dist_azimuth(stla, stlo, evla, evlo)
	dist = geo[0];
	baz = geo[1];
	gcarc = kilometer2degrees(dist/1000)

	if gcarc < 60 and gcarc > 30 and baz < 310 and baz > 290:
		#print([otime, gcarc, baz, evdp])
		arrival = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=gcarc, phase_list=['P'])
		arrival = arrival[0];
		ttime = arrival.time
		#print(ttime);
		atime = otime+ttime
		#atime_utc = UTCDateTime(atime)
		#print(atime_utc)
		#print(type(atime_utc))
		print(otime)
		try:
			st = client.get_waveforms(ntwk,stn,'*','BH*',atime-100,atime+1500,attach_response=True)
		except:
			print('No Data')
			continue
		ntr = len(st)
		print('Find '+str(ntr)+' Components')

		for tr in st:
			chan = tr.stats.channel
			print(chan)

		if ntr == 3:
			DeconvDecimate(st,pre_filt)

			yr = otime.year
			jd = otime.julday
			hr = otime.hour
			min = otime.minute
			sec = otime.second
			msec = otime.microsecond

			tr = st[0]
			chan = tr.stats.channel
			sacnm = dir+'/'+ntwk+'.'+stn+'.'+str(yr)+'-'+str(jd).zfill(3)+'-'+str(hr).zfill(2)+'-'+str(min).zfill(2)+'-'+str(sec).zfill(2)+'.'+chan+'.SAC'
			sac = SACTrace.from_obspy_trace(tr)
			sac.gcarc = gcarc
			sac.baz = baz
			sac.evlo = evlo
			sac.evla = evla
			sac.stlo = stlo
			sac.stla = stla
			sac.evdp = evdp
			sac.a = atime
			sac.o = otime
			sac.write(sacnm)
			print(sacnm)

			tr = st[1]
			chan = tr.stats.channel
			sacnm = dir+'/'+ntwk+'.'+stn+'.'+str(yr)+'-'+str(jd).zfill(3)+'-'+str(hr).zfill(2)+'-'+str(min).zfill(2)+'-'+str(sec).zfill(2)+'.'+chan+'.SAC'
			sac = SACTrace.from_obspy_trace(tr)
			sac.gcarc = gcarc
			sac.baz = baz
			sac.evlo = evlo
			sac.evla = evla
			sac.stlo = stlo
			sac.stla = stla
			sac.evdp = evdp
			sac.a = atime
			sac.o = otime
			sac.write(sacnm)
			print(sacnm)

			tr = st[2]
			chan = tr.stats.channel
			sacnm = dir+'/'+ntwk+'.'+stn+'.'+str(yr)+'-'+str(jd).zfill(3)+'-'+str(hr).zfill(2)+'-'+str(min).zfill(2)+'-'+str(sec).zfill(2)+'.'+chan+'.SAC'
			sac = SACTrace.from_obspy_trace(tr)
			sac.gcarc = gcarc
			sac.baz = baz
			sac.evlo = evlo
			sac.evla = evla
			sac.stlo = stlo
			sac.stla = stla
			sac.evdp = evdp
			sac.a = atime
			sac.o = otime
			sac.write(sacnm)
			print(sacnm)




#Net = net_sta.split('.')[0]
#Sta = net_sta.split('.')[1]
#pre_filt = (0.5,2,85,90)
#if Net=='BP':
#	chan = 'DP1'
#else:
#	chan = 'EHZ'
#try:
#	st = client.get_waveforms(Net,Sta,'*',chan,t-10,t+50,attach_response=True)
#	st.merge()
#	DeconvDecimate(st,pre_filt)
#	if Sta=='*':
#		for i in range(0,len(st)):
#			filename = Net+'_'+st[i].stats.station+'_'+date+'.SAC'
#			st[i].write(filename,format='SAC')
#	else:
#		filename = Net+'_'+Sta+'_'+date+'.SAC'
#		st.write(filename,format='SAC')
#except:
#	print('No data: %s.%s' %(Net,Sta))
