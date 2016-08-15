"""
FTP script for SEVIRI NRT data

*Created for use with ORACLES NASA ESPO mission*

Modification history
--------------------
Written: Michael Diamond, 8/15/2016, Seattle, WA
"""

import ftplib
import os
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import modipy as mod
import sevipy as sev
os.chdir('/Users/michaeldiamond/Documents/')
import datetime
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap

#Get today's date and current time
now = datetime.datetime.utcnow()
year = now.year
day = now.day
month = now.month
hour = now.hour
minute = now.minute
jday = mod.julian_day(month, day, year)

"""
Get LARC SEVIRI data
"""
user = 'oracles'
passwd = 'oracles2016'
host = 'cloudsgate2.larc.nasa.gov'

print 'Checking for new SEVIRI data...'
fdir = '/Users/michaeldiamond/Documents/oracles_files/msg/%s' % jday
os.chdir(fdir)
current_files = os.listdir(fdir)
new_files = []
path = '/prod/exp/oracles/d2/sat-ncdf/msg/%s/%s/%s/' % (year,month,day)
ftp = ftplib.FTP(host,user,passwd)
ftp.cwd(path)
ftp.set_pasv(True)

directory = ftp.pwd()
print 'Accessing directory %s%s\n' % (host,directory)

files = ftp.nlst()
files.reverse()
for f in files:
    if 5 <= int(f[14:16]) <= 19:
        C1 = f[-4] == 1
        C2 = f[-4] == 2 and f[-5] == 0
        if f not in current_files and np.logical_or(C1,C2):            
            print 'Getting file %s...' % f
            ftp.retrbinary('RETR %s' % f, open(f, 'wb').write)
            new_files.append(f)
            print 'Done!\n'

ftp.quit()

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












