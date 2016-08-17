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
import requests


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
user = 'oracles'
passwd = 'oracles2016'
host = 'http://cloudsgate2.larc.nasa.gov'

print 'Checking for new SEVIRI data...'

fdir = '/Users/michaeldiamond/Documents/oracles_files/msg/%s' % jday
os.chdir(fdir)
current_files = os.listdir(fdir)
for f in wanted:
    if f in current_files: wanted.remove(f)
new_files = []
path = '/prod/exp/oracles/d2/sat-ncdf/msg/%s/%s/%s/' % (year,month,day)
url = host+path
files = []

url_worked = False
try:
    r = requests.get(url,auth=(user,passwd))
    for f in wanted:
        try:
            with open(f, 'wb').write as code:
                code.write(r.content)
            new_files.append(f)
        except: pass
    url_worked = True
    r.close()
except: print 'MSG url request at %s failed...\n' % now

#Make plots
directory = '/Users/michaeldiamond/Documents/oracles/msg/%s' % jday
for f in new_files[::2]:
    #Read in file
    os.chdir(fdir)
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

"""
Test for cloud products
"""
os.chdir('/Users/michaeldiamond/Downloads')
url = 'http://cloudsgate2.larc.nasa.gov/prod/exp/oracles/d2/prod-ncdf/msg/2016/08//16/'
f = 'MET10.2016229.1045.cldprod.03km.pix.nc.gz'
r = requests.get(url)
with open(f, 'wb').write as code:
    code.write(r.content)






