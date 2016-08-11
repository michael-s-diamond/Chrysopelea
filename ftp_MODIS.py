"""
Test FTP script
"""

import ftplib
import os
import datetime
import numpy as np
from pyhdf import SD
from pyhdf.SD import SDC

#Get today's date and current time
now = datetime.datetime.utcnow()

year = now.year
day = now.day
month = now.month
hour = now.hour
minute = now.minute

def julian_day(month, day, year):
    """
    Convert calendar day to julian day.
    
    Parameters
    ----------
    month : int
    Numeric month (e.g., January = 1).
    
    day : int
    Day of the month.
    
    year : int
    Full year (e.g., 2016).
    
    Modification history
    --------------------
    Written: Michael Diamond, 08/08/2016, Seattle, WA
    """
    if not 1 <= month <= 12: return 'Error: Month must be an integer between 1 and 12.'
    if not 1 <= day <= 31: return 'Error: Day must be an integer between 1 and 31'
    if day > 28:
        if month == 2:
            if year%4 != 0: return 'Error: February %s only has 28 days.' % year
            elif year%4 == 0:
                if day > 29: return 'Error: February %s only has 29 days.' % year
        if day > 30:
            valid_months = [1,3,4,7,8,10,12]
            if month not in valid_months: return 'Error: This month does not have 31 days.'
    #Non-leap year
    if year%4 != 0:
        if month == 1:
            offset = 0
        elif month == 2:
            offset = 32-1
        elif month == 3:
            offset = 60-1
        elif month == 4:
            offset = 91-1
        elif month == 5:
            offset = 121-1
        elif month == 6:
            offset = 152-1
        elif month == 7:
            offset = 182-1
        elif month == 8:
            offset = 213-1
        elif month == 9:
            offset = 244-1
        elif month == 10:
            offset = 274-1
        elif month == 11:
            offset = 305-1
        elif month == 12:
            offset = 335-1
    #Leap year
    elif year%4 == 0:
        if month == 1:
            offset = 0
        elif month == 2:
            offset = 32-1
        elif month == 3:
            offset = 61-1
        elif month == 4:
            offset = 92-1
        elif month == 5:
            offset = 122-1
        elif month == 6:
            offset = 153-1
        elif month == 7:
            offset = 183-1
        elif month == 8:
            offset = 214-1
        elif month == 9:
            offset = 245-1
        elif month == 10:
            offset = 275-1
        elif month == 11:
            offset = 306-1
        elif month == 12:
            offset = 336-1
    return offset + day

jday = julian_day(month, day, year)

"""
Get LANCE NRT data
"""
user = 'diamond2'
passwd = 'WA7Murray'

#
###Terra
#
os.chdir('/Users/michaeldiamond/Documents/ORACLES/flight_planning/test')
current_files = os.listdir('/Users/michaeldiamond/Documents/ORACLES/flight_planning/test')

host = 'nrt1.modaps.eosdis.nasa.gov'

path = '/allData/1/MOD06_L2/%s/%s/' % (year,jday)

ftp = ftplib.FTP(host,user,passwd)

ftp.cwd(path)
ftp.set_pasv(True)

directory = ftp.pwd()
print 'Accessing directory %s%s\n' % (host,directory)

files = ftp.nlst()
files.reverse()
for f in files:
    if 8 <= int(f[18:20]) <= 12 and f not in current_files:
        if f[-1] == 't':
            ftp.retrbinary('RETR %s' % f, open(f, 'wb').write)
            #Check geolocation
            fi = open(f,'r')
            read = fi.read()
            fi.close()
            r = read.split()
            i_e = r.index('EASTBOUNDINGCOORDINATE')+6
            e = float(r[i_e])
            i_w = r.index('WESTBOUNDINGCOORDINATE')+6
            w = float(r[i_w])
            i_n = r.index('NORTHBOUNDINGCOORDINATE')+6
            n = float(r[i_n])
            i_s = r.index('SOUTHBOUNDINGCOORDINATE')+6
            s = float(r[i_s])
            #Reject files outside of ORACLES study region
            too_north = s > -4.5
            too_south = n < -25.5
            bad_lat = np.logical_or(too_north,too_south)
            too_east = w > 15.5
            too_west = e < -15.5
            bad_lon = np.logical_or(too_east,too_west)
            if bad_lat or bad_lon: files.remove(f[0:34])
        else:
            #Get file
            print 'Getting file %s...' % f
            ftp.retrbinary('RETR %s' % f, open(f, 'wb').write)
            print 'Done!\n'

ftp.quit()


#
###Aqua
#

path = '/allData/1/MYD06_L2/%s/%s/' % (year,jday)

ftp = ftplib.FTP(host,user,passwd)

ftp.cwd(path)
ftp.set_pasv(True)

directory = ftp.pwd()
print 'Accessing directory %s%s\n' % (host,directory)

files = ftp.nlst()
files.reverse()
for f in files:
    if 8 <= int(f[18:20]) <= 12 and f not in current_files:
        if f[-1] == 't':
            ftp.retrbinary('RETR %s' % f, open(f, 'wb').write)
            #Check geolocation
            fi = open(f,'r')
            read = fi.read()
            fi.close()
            r = read.split()
            i_e = r.index('EASTBOUNDINGCOORDINATE')+6
            e = float(r[i_e])
            i_w = r.index('WESTBOUNDINGCOORDINATE')+6
            w = float(r[i_w])
            i_n = r.index('NORTHBOUNDINGCOORDINATE')+6
            n = float(r[i_n])
            i_s = r.index('SOUTHBOUNDINGCOORDINATE')+6
            s = float(r[i_s])
            #Reject files outside of ORACLES study region
            too_north = s > -4.5
            too_south = n < -25.5
            bad_lat = np.logical_or(too_north,too_south)
            too_east = w > 15.5
            too_west = e < -15.5
            bad_lon = np.logical_or(too_east,too_west)
            if bad_lat or bad_lon: files.remove(f[0:34])
        else:
            #Get file
            print 'Getting file %s...' % f
            ftp.retrbinary('RETR %s' % f, open(f, 'wb').write)
            print 'Done!\n'

ftp.quit()
