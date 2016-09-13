"""
Process daily maps for MODIS at the end of each day

***Created for ORACLES NASA ESPO mission***

Modification history
--------------------
Written: Michael Diamond, 09/07/2016, Swakopmund, Namibia
"""

#Import libraries
import os
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import modipy as mod
reload(mod)
import datetime
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
from scipy.ndimage.interpolation import zoom
from matplotlib.colors import LogNorm

#Get today's date and current time
now = datetime.datetime.utcnow()
year = now.year
if now.month < 10: month = '0'+str(now.month)
else: month = str(now.month)
if now.day < 10: day = '0'+str(now.day)
else: day = str(now.day)
hour = now.hour
minute = now.minute
jday = mod.julian_day(now.month, now.day, year)

print 'Running MODIS daily mapmaker at %s' % now
plt.close("all")

"""
Terra
"""
print '\nTerra\n'
print 'Setting up maps...'
#Create maps
font = 'Arial'
size = 16
m = Basemap(llcrnrlon=-15.5,llcrnrlat=-25.5,urcrnrlon=15.5,urcrnrlat=-4.5,projection='merc',resolution='l')

#Delta ref
plt.figure(7)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummy have all the values needed
vmax = 6
vmed = 0
vmin = -6
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='RdYlBu_r')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[%sm]' % u"\u03BC", fontname=font,fontsize=size-2)
plt.title('%sref (2.1 %sm - 1.6 %sm) for %s/%s/%s from Terra' % (u"\u0394",u"\u03BC",u"\u03BC",month,day,year),\
fontname=font,fontsize=size)

#Del ref
plt.figure(14)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummy have all the values needed
vmax = 500
vmed = 0
vmin = -500
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='RdYlBu_r')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[per mil]', fontname=font,fontsize=size-2)
plt.title('%sref (2.1 %sm - 1.6 %sm) for %s/%s/%s from Terra' % (u"\u03B4",u"\u03BC",u"\u03BC",month,day,year),\
fontname=font,fontsize=size)

#Ref
plt.figure(21)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = 4*np.ones((2,2)) #Make the dummy have all the values needed
vmax = 24
vmed = 14
vmin = 4
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='viridis')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[%sm]' % u"\u03BC", fontname=font,fontsize=size-2)
plt.title('Effective radius (2.1 %sm) for %s/%s/%s from Terra' % (u"\u03BC",month,day,year),\
fontname=font,fontsize=size)

#COT
plt.figure(28)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummy have all the values needed
vmax = 32
vmed = 16
vmin = 0
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='viridis')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[unitless]', fontname=font,fontsize=size-2)
plt.title('Cloud optical thickness for %s/%s/%s from Terra' % (month,day,year),\
fontname=font,fontsize=size)

#Nd
plt.figure(35)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='steelblue')
m.fillcontinents(color='floralwhite',lake_color='steelblue',zorder=0)
dummy = np.ones((2,2)) #Make the dummy have all the values needed
vmax = 1000
vmed = 500
vmin = 1
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='cubehelix',norm = LogNorm(vmin=1, vmax=1000))
ticks = [1,10,100,1000]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[$\mathregular{cm^{-3}}$]', fontname=font,fontsize=size-2)
plt.title('Nd for %s/%s/%s from Terra' % (month,day,year),\
fontname=font,fontsize=size)

#Delta COT
plt.figure(42)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummy have all the values needed
vmax = 1
vmed = 0
vmin = -1
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='RdYlBu_r')
ticks = [vmin, (vmin+vmed)/2., vmed, (vmed+vmax)/2., vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[unitless]', fontname=font,fontsize=size-2)
plt.title('%sCOT for %s/%s/%s from Terra' % (u"\u0394",month,day,year),\
fontname=font,fontsize=size)

#Del COT
plt.figure(49)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummy have all the values needed
vmax = 100
vmed = 0
vmin = -100
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='RdYlBu_r')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[per mil]', fontname=font,fontsize=size-2)
plt.title('%sCOT for %s/%s/%s from Terra' % (u"\u03B4",month,day,year),\
fontname=font,fontsize=size)

#Delta Nd
plt.figure(56)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummy have all the values needed
vmax = 300
vmed = 0
vmin = -300
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='RdYlBu')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[per cc]', fontname=font,fontsize=size-2)
plt.title('%sNd for %s/%s/%s from Terra' % (u"\u0394",month,day,year),\
fontname=font,fontsize=size)

#Del Nd
plt.figure(63)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummary have all the values needed
vmax = 1000
vmed = 0
vmin = -1000
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='RdYlBu')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[per mil]', fontname=font,fontsize=size-2)
plt.title('%sNd for %s/%s/%s from Terra' % (u"\u03B4",month,day,year),\
fontname=font,fontsize=size)

#ACAOD
plt.figure(100)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummary have all the values needed
vmax = 3
vmed = 1.5
vmin = 0
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='inferno_r')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[unitless]', fontname=font,fontsize=size-2)
plt.title('ACAOD for %s/%s/%s from Terra' % (month,day,year),\
fontname=font,fontsize=size)

#ACAOD_ModAbsAero
plt.figure(107)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummary have all the values needed
vmax = 3
vmed = 1.5
vmin = 0
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='inferno_r')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[unitless]', fontname=font,fontsize=size-2)
plt.title('ACAOD_ModAbsAero for %s/%s/%s from Terra' % (month,day,year),\
fontname=font,fontsize=size)

print 'Done!\n'

#Set up Terra directories and get files
file_directory = '/Users/michaeldiamond/Documents/oracles_files/terra/%s' % jday
image_directory = '/Users/michaeldiamond/Documents/oracles/terra/%s' % jday
os.chdir(file_directory)
files = os.listdir(file_directory)
cloudfiles = []
for f in files:
    if f[-1] == 'f' and f[6] == 'L':
        cloudfiles.append(f)

#Make plots
for f in cloudfiles:   
    #Read in file
    os.chdir(file_directory)
    cloud = mod.nrtMOD06(f)
    time = cloud.time
    lon = zoom(cloud.lon,5.)
    lat = zoom(cloud.lat,5.)
    m = Basemap(llcrnrlon=-15.5,llcrnrlat=-25.5,urcrnrlon=15.5,urcrnrlat=-4.5,projection='merc',resolution='l')
    lon, lat = m(lon, lat)
    #Move to image directory
    os.chdir(image_directory)
    #Now add tile to daily maps
    print 'Adding data to daily maps...'
    print '...delta ref...'
    plt.figure(7)
    d = cloud.delta_ref16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu_r',vmin=-6,vmax=6)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_delta_ref' % (year,month,day),dpi=150)
    print '...del ref...'
    plt.figure(14)
    d = cloud.del_ref16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu_r',vmin=-500,vmax=500)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_del_ref' % (year,month,day),dpi=150)
    print '...effective radius...'
    plt.figure(21)
    d = cloud.ref
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=4,vmax=24)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_ref' % (year,month,day),dpi=150)
    print '...cloud optical thickness...'
    plt.figure(28)
    d = cloud.COT
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=0,vmax=32)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_cot' % (year,month,day),dpi=150)
    print '...Nd...'
    plt.figure(35)
    d = cloud.Nd
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='cubehelix',norm = LogNorm(vmin=1, vmax=1000))
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_Nd' % (year,month,day),dpi=150)
    print '...delta COT...'
    plt.figure(42)
    d = cloud.delta_COT16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu_r',vmin=-1,vmax=1)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_delta_cot' % (year,month,day),dpi=150)
    print '...del COT...'
    plt.figure(49)
    d = cloud.del_COT16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu_r',vmin=-100,vmax=100)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_del_cot' % (year,month,day),dpi=150)
    print '...delta Nd...'
    plt.figure(56)
    d = cloud.delta_Nd16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu',vmin=-300,vmax=300)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_delta_Nd' % (year,month,day),dpi=150)
    print '...del Nd...'
    plt.figure(63)
    d = cloud.del_Nd16
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='RdYlBu',vmin=-1000,vmax=1000)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_del_Nd' % (year,month,day),dpi=150)
    afile = 'MOD06ACAERO.A%s%s.%s.006.NRT.hdf' % (year,jday,time)
    if afile in files:
        #Now add tile to daily maps
        print '\nAdding aerosol data to daily maps...'
        print '...ACAOD map...'
        os.chdir(file_directory)
        aero = mod.nrtACAERO(afile)
        comp = mod.nrt_comp(f,afile)
        os.chdir(image_directory)
        plt.figure(100)
        d = aero.ds['Above_Cloud_AOD']
        plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='inferno_r',vmin=0,vmax=3)
        m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
        m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
        m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
        m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
        m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
        fig = plt.gcf()
        fig.set_size_inches(13.33,7.5)
        plt.savefig('%s_%s_%s_map_ACAOD' % (year,month,day),dpi=150)
        print '...ACAOD_ModAbsAero map...'
        plt.figure(107)
        d = aero.ds['Above_Cloud_AOD_ModAbsAero']
        plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='inferno_r',vmin=0,vmax=3)
        m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
        m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
        m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
        m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
        m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
        fig = plt.gcf()
        fig.set_size_inches(13.33,7.5)
        plt.savefig('%s_%s_%s_map_ACAOD_ModAbsAero' % (year,month,day),dpi=150)
        print 'Done!\n'
        print 'Making comparison plots...'
        for var in ['delta_ref16','delta_COT16','delta_Nd16','del_ref16','del_COT16','del_Nd16']:
            print '...%s...' % var
            plt.figure('comp')
            comp.compare(var)
            fig = plt.gcf()
            fig.set_size_inches(13.33,7.5)
            plt.savefig('%s_%s_%s_%s_comp_%s' % (year,month,day,time,var),dpi=100)
        print 'Done!\n'

plt.close("all")

"""
Aqua
"""
print '\nAqua\n'
print 'Setting up maps...'

#Ref
plt.figure(18)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = 4*np.ones((2,2)) #Make the dummy have all the values needed
vmax = 24
vmed = 14
vmin = 4
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='viridis')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[%sm]' % u"\u03BC", fontname=font,fontsize=size-2)
plt.title('Effective radius (2.1 %sm) for %s/%s/%s from Aqua' % (u"\u03BC",month,day,year),\
fontname=font,fontsize=size)

#COT
plt.figure(24)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummy have all the values needed
vmax = 32
vmed = 16
vmin = 0
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='viridis')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[unitless]', fontname=font,fontsize=size-2)
plt.title('Cloud optical thickness for %s/%s/%s from Aqua' % (month,day,year),\
fontname=font,fontsize=size)

#Nd
plt.figure(30)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='steelblue')
m.fillcontinents(color='floralwhite',lake_color='steelblue',zorder=0)
dummy = np.ones((2,2)) #Make the dummy have all the values needed
vmax = 1000
vmed = 500
vmin = 1
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='cubehelix',norm = LogNorm(vmin=1, vmax=1000))
ticks = [1,10,100,1000]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[$\mathregular{cm^{-3}}$]', fontname=font,fontsize=size-2)
plt.title('Nd for %s/%s/%s from Aqua' % (month,day,year),\
fontname=font,fontsize=size)

#ACAOD
plt.figure(101)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummy have all the values needed
vmax = 3
vmed = 1.5
vmin = 0
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='inferno_r')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[unitless]', fontname=font,fontsize=size-2)
plt.title('ACAOD for %s/%s/%s from Aqua' % (month,day,year),\
fontname=font,fontsize=size)

#ACAOD_ModAbsAero
plt.figure(106)
plt.clf()
m
m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents('k',zorder=0)
dummy = np.zeros((2,2)) #Make the dummary have all the values needed
vmax = 3
vmed = 1.5
vmin = 0
dummy[0,0] = vmax
dummy[-1,-1] = vmin
plt.pcolormesh(dummy,cmap='inferno_r')
ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
cbar = plt.colorbar(ticks = ticks)
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize=size-4)
cbar.set_label('[unitless]', fontname=font,fontsize=size-2)
plt.title('ACAOD_ModAbsAero for %s/%s/%s from Aqua' % (month,day,year),\
fontname=font,fontsize=size)

print 'Done!\n'

#Set up Aqua directories and get files
file_directory = '/Users/michaeldiamond/Documents/oracles_files/aqua/%s' % jday
image_directory = '/Users/michaeldiamond/Documents/oracles/aqua/%s' % jday
os.chdir(file_directory)
files = os.listdir(file_directory)
cloudfiles = []
for f in files:
    if f[-1] == 'f' and f[6] == 'L':
        cloudfiles.append(f)

#Make plots
for f in cloudfiles:   
    #Read in file
    os.chdir(file_directory)
    cloud = mod.nrtMOD06(f)
    time = cloud.time
    lon = zoom(cloud.lon,5.)
    lat = zoom(cloud.lat,5.)
    m = Basemap(llcrnrlon=-15.5,llcrnrlat=-25.5,urcrnrlon=15.5,urcrnrlat=-4.5,projection='merc',resolution='l')
    lon, lat = m(lon, lat)
    #Move to image directory
    os.chdir(image_directory)
    #Now add tile to daily maps
    print 'Adding data to daily maps...'
    print '...effective radius...'
    plt.figure(18)
    d = cloud.ref
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=4,vmax=24)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_ref' % (year,month,day),dpi=150)
    print '...cloud optical thickness...'
    plt.figure(24)
    d = cloud.COT
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='viridis',vmin=0,vmax=32)
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_cot' % (year,month,day),dpi=150)
    print '...Nd...'
    plt.figure(30)
    d = cloud.Nd
    plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='cubehelix',norm = LogNorm(vmin=1, vmax=1000))
    m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
    m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
    m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
    m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
    fig = plt.gcf()
    fig.set_size_inches(13.33,7.5)
    plt.savefig('%s_%s_%s_map_Nd' % (year,month,day),dpi=150)
    afile = 'MYD06ACAERO.A%s%s.%s.006.NRT.hdf' % (year,jday,time)
    if afile in files:
        #Now add tile to daily maps
        print '\nAdding aerosol data to daily maps...'
        print '...ACAOD map...'
        os.chdir(file_directory)
        aero = mod.nrtACAERO(afile)
        comp = mod.nrt_comp(f,afile)
        os.chdir(image_directory)
        plt.figure(101)
        d = aero.ds['Above_Cloud_AOD']
        plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='inferno_r',vmin=0,vmax=3)
        m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
        m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
        m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
        m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
        m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
        fig = plt.gcf()
        fig.set_size_inches(13.33,7.5)
        plt.savefig('%s_%s_%s_map_ACAOD' % (year,month,day),dpi=150)
        print '...ACAOD_ModAbsAero map...'
        plt.figure(106)
        d = aero.ds['Above_Cloud_AOD_ModAbsAero']
        plt.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],cmap='inferno_r',vmin=0,vmax=3)
        m.scatter(14.5247,-22.9390,s=250,c='orange',marker='D',latlon=True)
        m.scatter(-14.3559,-7.9467,s=375,c='c',marker='*',latlon=True)
        m.scatter(-5.7089,-15.9650,s=375,c='chartreuse',marker='*',latlon=True)
        m.plot([14.5247,13,0],[-22.9390,-23,-10],c='w',linewidth=5,linestyle='dashed',latlon=True)
        m.plot([14.5247,13,0],[-22.9390,-23,-10],c='k',linewidth=3,linestyle='dashed',latlon=True)
        fig = plt.gcf()
        fig.set_size_inches(13.33,7.5)
        plt.savefig('%s_%s_%s_map_ACAOD_ModAbsAero' % (year,month,day),dpi=150)
        print 'Done!\n'

plt.close("all")