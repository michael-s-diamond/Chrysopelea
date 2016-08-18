"""
Download script for SEVIRI NRT data

*Created for use with ORACLES NASA ESPO mission*

Modification history
--------------------
Written: Michael Diamond, 8/15/2016, Seattle, WA
"""

import os
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import modipy as mod
import sevipy as sev
os.chdir('/Users/michaeldiamond/Documents/')
import datetime
import numpy as np
import matplotlib.pylab as plt

#Get today's date and current time
now = datetime.datetime.utcnow()
year = now.year
day = now.day
month = now.month
hour = now.hour
minute = now.minute
jday = mod.julian_day(month, day, year)

#Figure out what files are wanted
wanted = []
for hr in ['05','06','07','08','09','10','11','12','13','14','15','16','17','18','19']:
    for mn in ['00','30']:
        wanted.append('MET10.%s%s.%s%s.03km.C01.nc' % (year,jday,hr,mn))
        wanted.append('MET10.%s%s.%s%s.03km.C02.nc' % (year,jday,hr,mn))

"""
Get LARC SEVIRI data
"""
file_directory = '/Users/michaeldiamond/Documents/oracles_files/msg/%s' % jday
os.chdir(file_directory)
current_files = os.listdir(file_directory)
try: os.system('perl retrieve_raw_%s.pl' % jday)
except: print 'Error getting raw counts at %s' % now
try: os.system('perl retrieve_prod_%s.pl' % jday)
except: print 'Error getting products at %s' % now
new_files = os.listdir(file_directory)
#Rid list of repeat files
for f in current_files:
    if f in new_files: new_files.remove(f)

#Make plots
directory = '/Users/michaeldiamond/Documents/oracles/msg/%s' % jday
for f in new_files:
    if 5 <= f[:0] <= 19:
    
    #Read in file
    os.chdir(file_directory)
    cr = sev.CR(f,f[0:26]+'2.nc')
    #Move to image directory
    os.chdir(directory)
    #Color ratio for each file pair
    print 'Making CR plots for %s...' % f
    plt.figure(100)
    cr.merc()
    fig = plt.gcf()
    fig.set_size_inches(2*13.33,2*7.5)
    plt.savefig('%s_%s_%s_%s_CRS' % (cr.year,cr.month,cr.day,cr.time),dpi=300)
    print 'Done!\n'




