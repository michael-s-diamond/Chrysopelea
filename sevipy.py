"""
Define basic functions for using SEVIRI data
"""

#Import libraries
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
import pysolar
import datetime

"""
General purpose functions
"""

#Convert Julian day to calendar day
def cal_day(julian_day,year):
    """
    Convert Julian day into calendar day.
    
    Parameters
    ----------
    julian_day : int
    Julian day.
    
    year : int
    Year.
    
    Return
    ------
    cal_day : tuple
    (month, day)
    """
    if julian_day <= 59:
        if julian_day <= 31:
            cal_day = (1,julian_day)
        elif julian_day > 31:
            cal_day = (2,julian_day-31)
    elif julian_day > 59:
        if year%4 != 0:
            if julian_day <= 90:
                cal_day = (3,julian_day-59)
            elif 90 < julian_day <= 120:
                cal_day = (4,julian_day-90)
            elif 120 < julian_day <= 151:
                cal_day = (5,julian_day-120)
            elif 151 < julian_day <= 181:
                cal_day = (6,julian_day-151)
            elif 181 < julian_day <= 212:
                cal_day = (7,julian_day-181)
            elif 212 < julian_day <= 243:
                cal_day = (8,julian_day-212)
            elif 243 < julian_day <= 273:
                cal_day = (9,julian_day-243)
            elif 273 < julian_day <= 304:
                cal_day = (10,julian_day-273)
            elif 304 < julian_day <= 334:
                cal_day = (11,julian_day-304)
            elif 334 < julian_day <= 365:
                cal_day = (12,julian_day-334)
            else:
                return 'Error: Non-leap year day must be between 1 and 365.'
        elif year%4 == 0:
            if julian_day == 60:
                cal_day = (2,29)
            elif 60 <= julian_day <= 91:
                cal_day = (3,julian_day-60)
            elif 91 < julian_day <= 121:
                cal_day = (4,julian_day-91)
            elif 121 < julian_day <= 152:
                cal_day = (5,julian_day-121)
            elif 152 < julian_day <= 182:
                cal_day = (6,julian_day-152)
            elif 182 < julian_day <= 213:
                cal_day = (7,julian_day-182)
            elif 213 < julian_day <= 244:
                cal_day = (8,julian_day-213)
            elif 244 < julian_day <= 274:
                cal_day = (9,julian_day-244)
            elif 274 < julian_day <= 305:
                cal_day = (10,julian_day-274)
            elif 305 < julian_day <= 335:
                cal_day = (11,julian_day-305)
            elif 335 < julian_day <= 366:
                cal_day = (12,julian_day-335)
            else:
                return 'Error: Leap year day must be between 1 and 366.'
    else:
        return 'Error: Inputs should be integers.'
    return cal_day

"""
Color ratios for ORACLES
"""

class CR(object):
    """
    Caculate C1:C2 color ratio from SEVERI data.
    
    Parameters
    ----------
    C1_file : string
    File name for channel 1 (600 nm).
    
    C2_file : string
    File name for channel 2 (800 nm).
    
    Return
    -------
    CR : array
    Color ratio of C1:C2.
    """
    
    def __init__(self,C1_file,C2_file):
        S_1 = 65.2296*10*(1/.56-1/.71)/(.71-.56) #Solar constant for 600 nm
        S_2 = 73.0127*10*(1/.74-1/.88)/(.88-.74) #Solar constant for 800 nm
        d = 1 #in AU
        #
        ###Load data
        #
        #Date
        self.day = int(C1_file[10:12+1])
        self.year = int(C1_file[6:9+1])
        self.time = C1_file[14:17+1]
        self.hour = int(self.time[0:2])
        self.minute = int(self.time[2:3])
        dsl = 1427 + (self.day - 153) #For 2016
        #Load channel 1 (600 nm)
        data_C1 = nc.Dataset(C1_file, 'r')
        count_C1 = data_C1['RAW'][:,:]
        g0_C1 = 0.5669
        g1_C1 = 1.23E-5
        self.lat = data_C1['Latitude'][:]/100.
        self.lon = data_C1['Longitude'][:]/100.
        #Load channel 2 (800 nm)
        data_C2 = nc.Dataset(C2_file, 'r')
        count_C2 = data_C2['RAW'][:,:]
        g0_C2 = 0.4529
        g1_C2 = 2.5E-6
        #
        ###Calculate solar zenith angle
        #
        sza = np.zeros([np.shape(self.lat)[0],np.shape(self.lon)[0]])
        for i in range(len(self.lat)):
            for j in range(len(self.lon)):
                date = datetime.datetime(self.year, cal_day(self.day, self.year)[0], \
                cal_day(self.day, self.year)[1],self.hour,self.minute,00)
                sza[i][j] = 90 - pysolar.solar.get_altitude_fast(self.lat[i],self.lon[j],date)
        self.sza = sza
        self.lon,self.lat = np.meshgrid(self.lon,self.lat)
        #
        ###Calculate radiances
        #
        self.Rad1 = (g0_C1+g1_C1*dsl)*(count_C1-51)
        self.Rad2 = (g0_C2+g1_C2*dsl)*(count_C2-51)
        #
        ###Calculate reflectances
        #
        self.R1 = (np.pi*self.Rad1*d**2)/(S_1*np.cos(sza*np.pi/180.))
        self.R2 = (np.pi*self.Rad2*d**2)/(S_2*np.cos(sza*np.pi/180.))
        self.CR = self.R2/self.R1
        #Close files
        data_C1.close()
        data_C2.close()
    
    def merc(self):
        """
        Plot the view as seen from the satellite
        """
        plt.clf()
        font = 16
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='i')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.CR,shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=.9,vmax=1.1)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 2:Channel 1 color ratio',fontsize=font-1)
        plt.title('Color ratio from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.show()
    
    def view(self):
        """
        Plot the view as seen from the satellite
        """
        plt.clf()
        m = Basemap(projection='nsper',lon_0=self.lon.mean(),lat_0=self.lat.mean(),resolution='l',satellite_height=36000*1000)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10))
        m.drawmeridians(np.arange(0,360,10))
        m.pcolormesh(self.lon,self.lat,self.CR,shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=.5,vmax=1.5)
        cbar = m.colorbar()
        cbar.set_label('Channel 1:Channel 2 color ratio')
        plt.title('View from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time))
        plt.show()
        
    def radmerc(self):
        """
        Plot the view as seen from the satellite
        """
        plt.clf()
        font = 16
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='i')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.Rad2/self.Rad1,shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=0.6,vmax=.8)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 2:Channel 1 color ratio',fontsize=font-1)
        plt.title('Radiance color ratio from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.show()
    
    def check(self):
        """
        Plot radiances and reflectances to make sure it all makes sense.
        """
        plt.clf()
        font = 12
        plt.subplot(2,2,1)
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='c')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.Rad1,shading='gouraud',cmap='viridis',latlon=True,vmin=0,vmax=300)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 1 radiance [W/m2/micron/sr]',fontsize=font-1)
        plt.title('600 nm radiance from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.subplot(2,2,2)
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='c')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.Rad2,shading='gouraud',cmap='plasma',latlon=True,vmin=0,vmax=300)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 2 radiance [W/m2/micron/sr]',fontsize=font-1)
        plt.title('800 nm radiance from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.subplot(2,2,3)
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='c')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.R1,shading='gouraud',cmap='YlGn',latlon=True,vmin=0,vmax=1)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 1 reflectance [unitless]',fontsize=font-1)
        plt.title('600 nm reflectance from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.subplot(2,2,4)
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='c')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.R2,shading='gouraud',cmap='YlOrRd',latlon=True,vmin=0,vmax=1)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 2 reflectance [unitless]',fontsize=font-1)
        plt.title('800 nm reflectance from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.show()
    
    def szaplot(self):
        """
        Plot the sza
        """
        plt.clf()
        font = 16
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='i')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.sza,shading='gouraud',cmap='magma_r',latlon=True,vmin=0,vmax=90)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Degrees',fontsize=font-1)
        plt.title('Solar zenith angle (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.show()

#Temporary
import os
os.chdir('/Users/michaeldiamond/Documents/ORACLES/SEVIRI')

C1f = 'MET10.2016216.1330.03km.C01.nc'
C2f = 'MET10.2016216.1330.03km.C02.nc'
cr = CR(C1f,C2f)

"""
Microphysics blend
"""
class blend(object):
    """
    Caculate C1:C2 color ratio from SEVERI data.
    
    Parameters
    ----------
    C1_file : string
    File name for channel 1 (600 nm).
    
    C2_file : string
    File name for channel 2 (800 nm).
    
    Return
    -------
    CR : array
    Color ratio of C1:C2.
    """
    
    def __init__(self,C1_file,C2_file):
        S_1 = 65.2296*10*(1/.56-1/.71)/(.71-.56) #Solar constant for 600 nm
        S_2 = 73.0127*10*(1/.74-1/.88)/(.88-.74) #Solar constant for 800 nm
        d = 1 #in AU
        #
        ###Load data
        #
        #Date
        self.day = int(C1_file[10:12+1])
        self.year = int(C1_file[6:9+1])
        self.time = C1_file[14:17+1]
        self.hour = int(self.time[0:2])
        self.minute = int(self.time[2:3])
        dsl = 1427 + (self.day - 153) #For 2016
        #Load channel 1 (600 nm)
        data_C1 = nc.Dataset(C1_file, 'r')
        count_C1 = data_C1['RAW'][:,:]
        g0_C1 = 0.5669
        g1_C1 = 1.23E-5
        self.lat = data_C1['Latitude'][:]/100.
        self.lon = data_C1['Longitude'][:]/100.
        #Load channel 2 (800 nm)
        data_C2 = nc.Dataset(C2_file, 'r')
        count_C2 = data_C2['RAW'][:,:]
        g0_C2 = 0.4529
        g1_C2 = 2.5E-6
        #
        ###Calculate solar zenith angle
        #
        sza = np.zeros([np.shape(self.lat)[0],np.shape(self.lon)[0]])
        for i in range(len(self.lat)):
            for j in range(len(self.lon)):
                date = datetime.datetime(self.year, cal_day(self.day, self.year)[0], \
                cal_day(self.day, self.year)[1],self.hour,self.minute,00)
                sza[i][j] = 90 - pysolar.solar.get_altitude_fast(self.lat[i],self.lon[j],date)
        self.sza = sza
        self.lon,self.lat = np.meshgrid(self.lon,self.lat)
        #
        ###Calculate radiances
        #
        self.Rad1 = (g0_C1+g1_C1*dsl)*(count_C1-51)
        self.Rad2 = (g0_C2+g1_C2*dsl)*(count_C2-51)
        #
        ###Calculate reflectances
        #
        self.R1 = (np.pi*self.Rad1*d**2)/(S_1*np.cos(sza*np.pi/180.))
        self.R2 = (np.pi*self.Rad2*d**2)/(S_2*np.cos(sza*np.pi/180.))
        self.CR = self.R2/self.R1
        #Close files
        data_C1.close()
        data_C2.close()
    
    def merc(self):
        """
        Plot the view as seen from the satellite
        """
        plt.clf()
        font = 16
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='i')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.CR,shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=.9,vmax=1.1)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 2:Channel 1 color ratio',fontsize=font-1)
        plt.title('Color ratio from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.show()
    
    def view(self):
        """
        Plot the view as seen from the satellite
        """
        plt.clf()
        m = Basemap(projection='nsper',lon_0=self.lon.mean(),lat_0=self.lat.mean(),resolution='l',satellite_height=36000*1000)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10))
        m.drawmeridians(np.arange(0,360,10))
        m.pcolormesh(self.lon,self.lat,self.CR,shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=.5,vmax=1.5)
        cbar = m.colorbar()
        cbar.set_label('Channel 1:Channel 2 color ratio')
        plt.title('View from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time))
        plt.show()
        
    def radmerc(self):
        """
        Plot the view as seen from the satellite
        """
        plt.clf()
        font = 16
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='i')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.Rad2/self.Rad1,shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=0.6,vmax=.8)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 2:Channel 1 color ratio',fontsize=font-1)
        plt.title('Radiance color ratio from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.show()
    
    def check(self):
        """
        Plot radiances and reflectances to make sure it all makes sense.
        """
        plt.clf()
        font = 12
        plt.subplot(2,2,1)
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='c')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.Rad1,shading='gouraud',cmap='viridis',latlon=True,vmin=0,vmax=300)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 1 radiance [W/m2/micron/sr]',fontsize=font-1)
        plt.title('600 nm radiance from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.subplot(2,2,2)
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='c')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.Rad2,shading='gouraud',cmap='plasma',latlon=True,vmin=0,vmax=300)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 2 radiance [W/m2/micron/sr]',fontsize=font-1)
        plt.title('800 nm radiance from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.subplot(2,2,3)
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='c')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.R1,shading='gouraud',cmap='YlGn',latlon=True,vmin=0,vmax=1)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 1 reflectance [unitless]',fontsize=font-1)
        plt.title('600 nm reflectance from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.subplot(2,2,4)
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='c')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.R2,shading='gouraud',cmap='YlOrRd',latlon=True,vmin=0,vmax=1)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Channel 2 reflectance [unitless]',fontsize=font-1)
        plt.title('800 nm reflectance from MSG SEVIRI (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.show()
    
    def szaplot(self):
        """
        Plot the sza
        """
        plt.clf()
        font = 16
        m = Basemap(llcrnrlon=self.lon.min(),llcrnrlat=self.lat.min(),urcrnrlon=self.lon.max(),\
        urcrnrlat=self.lat.max(),projection='merc',resolution='i')
        m.drawmapboundary(linewidth=1.5)        
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,5),labels=[1,0,0,0],fontsize=font-2)
        m.drawmeridians(np.arange(0,360,5),labels=[1,1,0,1],fontsize=font-2)
        m.pcolormesh(self.lon,self.lat,self.sza,shading='gouraud',cmap='magma_r',latlon=True,vmin=0,vmax=90)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=font-2) 
        cbar.set_label('Degrees',fontsize=font-1)
        plt.title('Solar zenith angle (%s/%s/%s, %s UTC)' % \
        (cal_day(int(self.day),int(self.year))[0],cal_day(int(self.day),int(self.year))[1],\
        self.year,self.time),fontsize=font+2)
        plt.show()
