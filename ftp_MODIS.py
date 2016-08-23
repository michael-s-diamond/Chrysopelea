"""
FTP script for MODIS NRT data

*Created for use with ORACLES NASA ESPO mission*

Modification history
--------------------
Written: Michael Diamond, 8/11/2016, Seattle, WA
Modified: Michael Diamond, 8/16/2016, Seattle, WA
    -Made script resilient to ftp failures
"""

import ftplib
import os
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import modipy as mod
os.chdir('/Users/michaeldiamond/')
from login import u, p
import datetime
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
from scipy.ndimage.interpolation import zoom

#Get today's date and current time
now = datetime.datetime.utcnow()
year = now.year
day = now.day
month = now.month
hour = now.hour
minute = now.minute
jday = mod.julian_day(month, day, year)

"""
Get LANCE NRT data
"""
user = u['MODIS']
passwd = p['MODIS']
host = 'nrt1.modaps.eosdis.nasa.gov'

#
###Terra cloud
#
print 'Checking for new Terra data...'
fdir = '/Users/michaeldiamond/Documents/oracles_files/terra/%s' % jday
os.chdir(fdir)
current_files = os.listdir(fdir)
new_files = []
path = '/allData/1/MOD06_L2/%s/%s/' % (year,jday)
files = []

#Try to access
ftp_worked = False
try:
    ftp = ftplib.FTP(host,user,passwd,timeout=60)
    ftp.cwd(path)
    ftp.set_pasv(True)

    directory = ftp.pwd()
    print 'Accessing directory %s%s\n' % (host,directory)

    files = ftp.nlst()
    ftp_worked = True
except: print 'Terra ftp at %s failed...\n' % now

files.reverse()
for f in files:
    if 8 <= int(f[18:20]) <= 12:
        if f not in current_files:
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
                new_files.append(f)
                print 'Done!\n'
        else:
            files.remove(f[0:34])

if ftp_worked: ftp.quit()


directory = '/Users/michaeldiamond/Documents/oracles/terra/%s' % jday
os.chdir(directory)
current_images = os.listdir(directory)
os.chdir(fdir)

#Add to new_files if things got overlooked before and no image was made
for f in current_files:
    if f[-1] == 'f':
        time = f[18:22]
        ref_im = '%s_%s_%s_%s_tri_ref.png' % (year,month,day,time)
        if ref_im not in current_images: new_files.append(f)
    else: pass

rsync = False
if len(new_files) > 0: rsync = True

#Make plots
for f in new_files:
    #Read in file
    os.chdir(fdir)
    cloud = mod.nrtMOD06(f)
    lon = zoom(cloud.lon,5.)
    lat = zoom(cloud.lat,5.)
    m = Basemap(llcrnrlon=-15.5,llcrnrlat=-25.5,urcrnrlon=15.5,urcrnrlat=-4.5,projection='merc',resolution='l')
    lon, lat = m(lon, lat)
    #Move to image directory
    os.chdir(directory)
    #Triplots for each file
    print 'Making triplots for %s...' % f
    plt.figure(3)
    print '...ref...'
    plt.clf()
    cloud.triplot(data='ref',full_res=True,num=3)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_%s_tri_ref' % (cloud.year,mod.month_num[cloud.month],cloud.day,cloud.time),dpi=150)
    print '...geo...'
    plt.clf()
    cloud.triplot(data='geo',full_res=False,num=3)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_%s_tri_geo' % (cloud.year,mod.month_num[cloud.month],cloud.day,cloud.time),dpi=150)
    print '...cot...'
    plt.clf()
    cloud.triplot(data='cot',full_res=True,num=3)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_%s_tri_cot' % (cloud.year,mod.month_num[cloud.month],cloud.day,cloud.time),dpi=150)
    print '...Nd...'
    plt.clf()
    cloud.triplot(data='Nd',full_res=True,num=3)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_%s_tri_Nd' % (cloud.year,mod.month_num[cloud.month],cloud.day,cloud.time),dpi=150)
    print 'Done!\n'
    #Now add tile to daily maps
    print 'Adding data to daily maps...'
    print '...delta ref...'
    plt.figure(7)
    d = cloud.delta_ref16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu_r',vmin=-6,vmax=6)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_delta' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print '...del ref...'
    plt.figure(14)
    d = cloud.del_ref16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu_r',vmin=-500,vmax=500)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_del' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print '...effective radius...'
    plt.figure(21)
    d = cloud.ref
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=4,vmax=24)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_ref' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print '...cloud optical thickness...'
    plt.figure(28)
    d = cloud.COT
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=0,vmax=32)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_cot' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print '...Nd...'
    plt.figure(35)
    d = cloud.Nd
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=0,vmax=1000)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_Nd' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print 'Done!\n'
    
#
###Aqua cloud
#
print 'Checking for new Aqua data...'
fdir = '/Users/michaeldiamond/Documents/oracles_files/aqua/%s' % jday
os.chdir(fdir)
current_files = os.listdir(fdir)
new_files = []
path = '/allData/1/MYD06_L2/%s/%s/' % (year,jday)
files = []

#Try to access
ftp_worked = False
try:
    ftp = ftplib.FTP(host,user,passwd,timeout=60)
    ftp.cwd(path)
    ftp.set_pasv(True)

    directory = ftp.pwd()
    print 'Accessing directory %s%s\n' % (host,directory)

    files = ftp.nlst()
    ftp_worked = True
except: print 'Aqua ftp at %s failed...\n' % now

files.reverse()
for f in files:
    if 12 <= int(f[18:20]) <= 15:
        if f not in current_files:
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
                new_files.append(f)
                print 'Done!\n'
        else:
            files.remove(f[0:34])

if ftp_worked: ftp.quit()

directory = '/Users/michaeldiamond/Documents/oracles/aqua/%s' % jday
os.chdir(directory)
current_images = os.listdir(directory)
os.chdir(fdir)

#Add to new_files if things got overlooked before and no image was made
for f in current_files:
    if f[-1] == 'f':
        time = f[18:22]
        ref_im = '%s_%s_%s_%s_tri_ref.png' % (year,month,day,time)
        if ref_im not in current_images: new_files.append(f)
    else: pass

if len(new_files) > 0: rsync = True

#Make plots
for f in new_files:
    #Read in file
    os.chdir(fdir)
    cloud = mod.nrtMOD06(f)
    lon = zoom(cloud.lon,5.)
    lat = zoom(cloud.lat,5.)
    m = Basemap(llcrnrlon=-15.5,llcrnrlat=-25.5,urcrnrlon=15.5,urcrnrlat=-4.5,projection='merc',resolution='l')
    lon, lat = m(lon, lat)
    #Move to image directory
    os.chdir(directory)
    #Triplots for each file
    print 'Making triplots for %s...' % f
    plt.figure(3)
    print '...ref...'
    plt.clf()
    cloud.triplot(data='ref',full_res=True,num=3)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_%s_tri_ref' % (cloud.year,mod.month_num[cloud.month],cloud.day,cloud.time),dpi=150)
    print '...geo...'
    plt.clf()
    cloud.triplot(data='geo',full_res=False,num=3)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_%s_tri_geo' % (cloud.year,mod.month_num[cloud.month],cloud.day,cloud.time),dpi=150)
    print '...cot...'
    plt.clf()
    cloud.triplot(data='cot',full_res=True,num=3)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_%s_tri_cot' % (cloud.year,mod.month_num[cloud.month],cloud.day,cloud.time),dpi=150)
    print '...Nd...'
    plt.clf()
    cloud.triplot(data='Nd',full_res=True,num=3)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_%s_tri_Nd' % (cloud.year,mod.month_num[cloud.month],cloud.day,cloud.time),dpi=150)
    print 'Done!\n'
    #Now add tile to daily maps
    print 'Adding data to daily maps...'
    print '...delta ref...'
    plt.figure(6)
    d = cloud.delta_ref16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu_r',vmin=-6,vmax=6)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_delta' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print '...del ref...'
    plt.figure(12)
    d = cloud.del_ref16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu_r',vmin=-500,vmax=500)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_del' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print '...effective radius...'
    plt.figure(18)
    d = cloud.ref
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=4,vmax=24)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_ref' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print '...cloud optical thickness...'
    plt.figure(24)
    d = cloud.COT
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=0,vmax=32)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_cot' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print '...Nd...'
    plt.figure(30)
    d = cloud.Nd
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=0,vmax=1000)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.plot([14.5247,0,-10],[-22.9390,-10,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_Nd' % (cloud.year,mod.month_num[cloud.month],cloud.day),dpi=150)
    print 'Done!\n'

if rsync:
    try: os.system('rsync -a /Users/michaeldiamond/Documents/oracles diamond2@olympus.atmos.washington.edu:~/public_html')
    except: print 'Rsync failed at %s' % now