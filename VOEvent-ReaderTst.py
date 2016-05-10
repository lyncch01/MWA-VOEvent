import sys
import voeventparse
import datetime
import logging
from fourpiskytools.notify import Notifier
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

def handle_voevent(v):
	role = v.attrib['role']
	if role == 'observation':
		if is_grb(v):
			return = handle_grb(v)
		elif is_trans(v):
			return handle_trans(v)
	else:
		return handle_other(v)
		

def is_grb(v):
	ivorn = v.attrib['ivorn']
	if ivorn.find("ivo://nasa.gsfc.tan/vo-gcn.gsfc/SWIFT#BAT_GRB_Pos") == 0:
		return True
	else:
		return False

def is_trans(v):
	ivorn = v.attrib['ivorn']
	if ivorn.find("ivo://nasa.gsfc.tan/vo-gcn.gsfc/SWIFT#BAT_TRANS") == 0:
		return True
	else:
		return False

def handle_grb(v):
	srcRA = float(v.find(".//C1")) #degree
	srcDec = float(v.find(".//C2")) #degree
	date = (float(v.find(".//Param[@name='Burst_TJD']").attrib['value'])+ 2440000.5)*u.day
	Time = float(v.find(".//Param[@name='Burst_SOD']").attrib['value'])*u.second
	dateTime = date + Time
	fovVal, t0, t1 =  in_fov(srcRA, srcDec, dateTime)
	if fovVal == True:
		trigID = v.find(".//Param[@name = 'TrigID']").attrib['value']
		logging.basicConfig(filename='grbtst.log', level=logging.INFO)
		logger=logging.getLogger('grbtst')
		logger.handlers.append(logging.StreamHandler(sys.stdout))
		text = 'Swift packet recieved and identified as GRB \n'
		text1 = 'Burst occured at '+str(dateTime)+'\n'
		text2 = 'Coords are RA='+str(srcRA)+' deg DEC='+str(srcDec)+' deg \n'
		text3 = 'Source will be at an elevation > 45 degrees at MWA from '+t0+' to '+t1
		textTot = text + text1 + text2 + text3
		n = Notifier()
		n.send_notification(title='NEW SWIFT GRB!', text=textTot)
		return True, trigID, dateTime, srcRA, srcDec, t0, t1
	else:
		return False, 0., 0., 0., 0., 0., 0.
	
def handle_trans(v):
	pointSrc = v.find(".//Param[@name='Point_Source']").attrib['value']
	tgtInt = v.find(".//Param[@name='Target_of_Interest']").attrib['value']
	tgtFltCat = v.find(".//Param[@name='Target_in_Flt_Catalog']").attrib['value']
	if pointSrc == 'true' and tgtInt == 'true' and tgtFltCat =='true':
		srcRA = float(v.find(".//C1")) #degree
		srcDec = float(v.find(".//C2")) #degree
		date = (float(v.find(".//Param[@name='Burst_TJD']").attrib['value'])+ 2440000.5)*u.day
		Time = float(v.find(".//Param[@name='Burst_SOD']").attrib['value'])*u.second
		dateTime = date + Time
		fovVal, t0, t1 = in_fov(srcRA, srcDec, dateTime)
		if fovVal == True:
			trigID = v.find(".//Param[@name = 'TrigID']").attrib['value']
			logging.basicConfig(filename='grbtst.log', level=logging.INFO)
			logger=logging.getLogger('grbtst')
			logger.handlers.append(logging.StreamHandler(sys.stdout))
			text = 'Swift packet recieved and could be a stellar flare \n'
			text1 = 'Burst occured at '+str(dateTime)+'\n'
			text2 = 'Coords are RA='+str(srcRA)+' deg DEC='+str(srcDec)+' deg \n'
			text3 = 'Source will be at an elevation > 45 degrees at MWA from '+t0+' to '+t1
			textTot = text + text1 + text2 + text3
			n = Notifier()
			n.send_notification(title='NEW SWIFT TRANS EVENT!', text=textTot)
			return True, trigID, dateTime, srcRA, srcDec, t0, t1
		else:
			return False, 0., 0., 0., 0., 0., 0.
	else:
		return False, 0., 0., 0., 0., 0., 0.
			
def handle_other(v):
	return False, 0., 0., 0., 0., 0., 0.	
		

def in_fov(ra, dec, date):
	#Determines if the sources is in the sky covered by the MWA given its RA, Dec and date of observation
	
	utcoffset = 8.0*u.hour #Australian Western standard time
	startJD = Time(date, format='jd')
	startUTC = startJD.utc + utcoffset
	delta_time = np.linspace(0, 12, 100)*u.hour
	times = startUTC + delta_time
	srcAlt = []
	for i in times:
		srcAlt.append(calc_elev(ra, dec, i))
	srcAlt = np.array(srcAlt)
	
	#45 degrees lowest elevation (check with David!)
	indxAbove = np.where(srcAlt > 45.0)[0]
	srcAbove = srcAlt[indxAbove]
	timeAbove = times[indxAbove]
	if len(srcAbove) == 0:
		return False, 0, 0
	else:
		return True, timeAbove[0].iso, timeAbove[-1].iso
		
def calc_elev(ra, dec, time_int):
	#Calculates the elevation of the source with RA and Dec given a time time_int
	srcLoc = SkyCoord(ra*u.degree, dec*u.degree,frame='icrs')
	mwaLoc = EarthLocation(lat=-26.7031*u.degree, lon=116.671*u.degree, height=377.827*u.m)
	time = Time(time_int)
	srcAltAz = srcLoc.transform_to(AltAz(obstime=time, location=mwaLoc))
	return srcAltAz.alt.value