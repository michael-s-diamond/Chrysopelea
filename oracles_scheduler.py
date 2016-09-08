"""
Retrieve MODIS and SEVIRI data and plot at specified times

***Written for ORACLES NASA ESPO mission***

Modification history
--------------------
Written: Michael Diamond, Seattle, WA
Modified: Michael Diamond, 09/08/2016, Swakopmund, Namibia
    -Major streamlining of MODIS data
"""

import schedule
import time
import datetime
import os
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import ftp_MODIS
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import ftp_SEVIRI
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import daily_MODIS
os.chdir('/Users/michaeldiamond/')


#What would be needed to run mapmaker
def run_daily_MODIS():
    print 'Running daily_MODIS...'
    os.chdir('/Users/michaeldiamond')
    reload(daily_MODIS)
    print 'Done!\n'
    
schedule.every().day.at("15:00").do(run_daily_MODIS)

#Run ftp_MODIS every 5 min
def run_ftp_MODIS():
    print 'Running ftp_MODIS at %s:%s...' % (datetime.datetime.utcnow().hour,datetime.datetime.utcnow().minute)
    os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
    try: reload(ftp_MODIS)
    except: 
        print 'Something went wrong at %s' % (datetime.datetime.utcnow())
    os.chdir('/Users/michaeldiamond')
    print 'Done!\n'

schedule.every(5).minutes.do(run_ftp_MODIS)

#Run ftp_SEVIRI every 15 min
def run_ftp_SEVIRI():
    print 'Running ftp_SEVIRI at %s:%s...' % (datetime.datetime.utcnow().hour,datetime.datetime.utcnow().minute)
    os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
    try: reload(ftp_SEVIRI)
    except: 
        print 'Something went wrong at %s' % (datetime.datetime.utcnow())
    os.chdir('/Users/michaeldiamond')
    print 'Done!\n'

schedule.every(5).minutes.do(run_ftp_SEVIRI)

#Loop forever and ever and ever...
while True:
    schedule.run_pending()
    time.sleep(10)

