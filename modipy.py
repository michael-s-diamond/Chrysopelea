"""
Define simple functions to analyze MODIS data using Python.

Modification history
--------------------
Written: Michael Diamond, 8/3-4/16, Seattle, WA
"""

#Import libraries
from pyhdf import SD
from pyhdf.SD import SDC
import numpy as np
import numpy.ma as ma
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap

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
    
    Returns
    -------
    cal_day : string
    Calendar month and day (no year).
    
    Modification history
    --------------------
    Written: Michael Diamond, 8/3/16, Seattle, WA
    """
    if julian_day <= 59:
        if julian_day <= 31:
            cal_day = 'January %s' % (julian_day)
        elif julian_day > 31:
            cal_day = 'February %s' % (julian_day-31)
    elif julian_day > 59:
        if year%4 != 0:
            if julian_day <= 90:
                cal_day = 'March %s' % (julian_day-59)
            elif 90 < julian_day <= 120:
                cal_day = 'April %s' % (julian_day-90)
            elif 120 < julian_day <= 151:
                cal_day = 'May %s' % (julian_day-120)
            elif 151 < julian_day <= 181:
                cal_day = 'June %s' % (julian_day-151)
            elif 181 < julian_day <= 212:
                cal_day = 'July %s' % (julian_day-181)
            elif 212 < julian_day <= 243:
                cal_day = 'August %s' % (julian_day-212)
            elif 243 < julian_day <= 273:
                cal_day = 'September %s' % (julian_day-243)
            elif 273 < julian_day <= 304:
                cal_day = 'October %s' % (julian_day-273)
            elif 304 < julian_day <= 334:
                cal_day = 'November %s' % (julian_day-304)
            elif 334 < julian_day <= 365:
                cal_day = 'December %s' % (julian_day-334)
            else:
                return 'Error: Non-leap year day must be between 1 and 365.'
        elif year%4 == 0:
            if julian_day == 60:
                cal_day = 'February 29'
            elif 60 <= julian_day <= 91:
                cal_day = 'March %s' % (julian_day-60)
            elif 91 < julian_day <= 121:
                cal_day = 'April %s' % (julian_day-91)
            elif 121 < julian_day <= 152:
                cal_day = 'May %s' % (julian_day-121)
            elif 152 < julian_day <= 182:
                cal_day = 'June %s' % (julian_day-152)
            elif 182 < julian_day <= 213:
                cal_day = 'July %s' % (julian_day-182)
            elif 213 < julian_day <= 244:
                cal_day = 'August %s' % (julian_day-213)
            elif 244 < julian_day <= 274:
                cal_day = 'September %s' % (julian_day-244)
            elif 274 < julian_day <= 305:
                cal_day = 'October %s' % (julian_day-274)
            elif 305 < julian_day <= 335:
                cal_day = 'November %s' % (julian_day-305)
            elif 335 < julian_day <= 366:
                cal_day = 'December %s' % (julian_day-335)
            else:
                return 'Error: Leap year day must be between 1 and 366.'
    else:
        return 'Error: Inputs should be integers.'
    return cal_day

#Average wavelength and wavelength spread for each spectral band
def avg_wavelength(band):
    """
    Average wavelength for given MODIS spectral band.
    
    Parameters
    ----------
    band : int
    MODIS band, from 1 to 36.
    
    Returns
    -------
    For band < 20: average wavelength in nm.
    For band >= 20: average wavelength in microns.
    
    Modification history
    --------------------
    Written: Michael Diamond, 8/3/16, Seattle, WA
    """
    if band == 1:
        return (670 + 620)/2. 
    elif band == 2:
        return (841 + 876)/2. 
    elif band == 3:
        return (459 + 479)/2. 
    elif band == 4:
        return (545 + 565)/2. 
    elif band == 5:
        return (1230 + 1250)/2. 
    elif band == 6:
        return (1628 + 1652)/2. 
    elif band == 7:
        return (2105 + 2155)/2. 
    elif band == 8:
        return (405 + 420)/2. 
    elif band == 9:
        return (438 + 448)/2. 
    elif band == 10:
        return (483 + 493)/2. 
    elif band == 11:
        return (526 + 536)/2. 
    elif band == 12:
        return (546 + 556)/2. 
    elif band == 13:
        return (662 + 672)/2. 
    elif band == 14:
        return (673 + 683)/2. 
    elif band == 15:
        return (743 + 753)/2. 
    elif band == 16:
        return (862 + 877)/2. 
    elif band == 17:
        return (890 + 920)/2. 
    elif band == 18:
        return (931 + 941)/2. 
    elif band == 19:
        return (915 + 965)/2. 
    elif band == 20:
        return (3.660 + 3.840)/2. 
    elif band == 21:
        return (3.929 + 3.989)/2. 
    elif band == 22:
        return (3.929 + 3.989)/2. 
    elif band == 23:
        return (4.020 + 4.080)/2. 
    elif band == 24:
        return (4.433 + 4.498)/2. 
    elif band == 25:
        return (4.482 + 4.549)/2. 
    elif band == 26:
        return (1.360 + 1.390)/2. 
    elif band == 27:
        return (6.535 + 6.895)/2. 
    elif band == 28:
        return (7.175 + 7.475)/2. 
    elif band == 29:
        return (8.400 + 8.700)/2. 
    elif band == 30:
        return (9.580 + 9.880)/2. 
    elif band == 31:
        return (10.780 + 11.280)/2. 
    elif band == 32:
        return (11.770 + 12.270)/2. 
    elif band == 33:
        return (13.185 + 13.485)/2. 
    elif band == 34:
        return (13.485 + 13.785)/2. 
    elif band == 35:
        return (13.785 + 14.085)/2. 
    elif band == 36:
        return (14.085 + 14.385)/2. 
    else:
        print 'Error: Band must be an integer between 1-36.'

def channel_width(band):
    """
    Wavelength spread for given MODIS spectral band.
    
    Parameters
    ----------
    band : int
    MODIS band, from 1 to 36.
    
    Returns
    -------
    For band < 20: spectral width in nm.
    For band >= 20: spectral width in microns.
    
    Modification history
    --------------------
    Written: Michael Diamond, 8/3/16, Seattle, WA
    """
    if band == 1:
        return np.abs(670 - 620) 
    elif band == 2:
        return np.abs(841 - 876) 
    elif band == 3:
        return np.abs(459 - 479) 
    elif band == 4:
        return np.abs(545 - 565) 
    elif band == 5:
        return np.abs(1230 - 1250) 
    elif band == 6:
        return np.abs(1628 - 1652) 
    elif band == 7:
        return np.abs(2105 - 2155) 
    elif band == 8:
        return np.abs(405 - 420) 
    elif band == 9:
        return np.abs(438 - 448) 
    elif band == 10:
        return np.abs(483 - 493) 
    elif band == 11:
        return np.abs(526 - 536) 
    elif band == 12:
        return np.abs(546 - 556) 
    elif band == 13:
        return np.abs(662 - 672) 
    elif band == 14:
        return np.abs(673 - 683) 
    elif band == 15:
        return np.abs(743 - 753) 
    elif band == 16:
        return np.abs(862 - 877) 
    elif band == 17:
        return np.abs(890 - 920) 
    elif band == 18:
        return np.abs(931 - 941) 
    elif band == 19:
        return np.abs(915 - 965) 
    elif band == 20:
        return np.abs(3.660 - 3.840) 
    elif band == 21:
        return np.abs(3.929 - 3.989) 
    elif band == 22:
        return np.abs(3.929 - 3.989) 
    elif band == 23:
        return np.abs(4.020 - 4.080) 
    elif band == 24:
        return np.abs(4.433 - 4.498) 
    elif band == 25:
        return np.abs(4.482 - 4.549) 
    elif band == 26:
        return np.abs(1.360 - 1.390) 
    elif band == 27:
        return np.abs(6.535 - 6.895) 
    elif band == 28:
        return np.abs(7.175 - 7.475) 
    elif band == 29:
        return np.abs(8.400 - 8.700) 
    elif band == 30:
        return np.abs(9.580 - 9.880) 
    elif band == 31:
        return np.abs(10.780 - 11.280) 
    elif band == 32:
        return np.abs(11.770 - 12.270) 
    elif band == 33:
        return np.abs(13.185 - 13.485) 
    elif band == 34:
        return np.abs(13.485 - 13.785) 
    elif band == 35:
        return np.abs(13.785 - 14.085) 
    elif band == 36:
        return np.abs(14.085 - 14.385) 
    else:
        print 'Error: Band must be an integer between 1-36.'

#Characteristic colorbar set by band
def colorbar(band):
    """
    Get colorbar based on MODIS spectral band.
    
    Parameters
    ----------
    band : int
    Spectral band number.
    
    Returns
    -------
    String of default matplotlib to use for a given band.
    
    Modification history
    --------------------
    Written: Michael Diamond, 08/04/2016, Seattle, WA
    """
    if band == 1 or 13 <= band <= 15:
        return 'OrRd'
    elif band == 2 or 5 <= band <= 7 or 16 <= band <= 19:
        return 'YlOrRd'
    elif band == 3 or 8 <= band <= 10:
        return 'PuBu'
    elif band == 4 or 11 <= band <= 12:
        return 'YlGn'
    else:
        return 'cubehelix'

"""
Class for MOD021KM/MYD021KM calibrated radiance/reflectance files from LAADS Web (https://ladsweb.nascom.nasa.gov/)
"""

class MOD021KM(object):
    """
    Create object to analyze a MOD021KM/MYD021KM calibrated radiance/reflectance file from LAADS Web (https://ladsweb.nascom.nasa.gov/).
    
    Parameters
    ----------
    filename : string
    Name of MOD021KM/MYD021KM file to analyze.
    
    Methods
    -------
    radiance: Get array of calibrated radiances.
    
    reflectance: Get array of calibrated reflectances.
    
    brightness_temperature: Get array of brightness temperature.
        Also, Tb() returns equivalent array.
    
    quick_plot: Fast plot of data. Intended as a check/first look at the data.
    
    Returns
    -------
    filename, month, time, satellite: string
    Name of MOD021KM/MYD021KM file used, month, time of retrieval, satellite (Terra or Aqua).
    
    jday, year, day : int
    Julian day, year, and calendar day.
    
    lon, lat : array
    Longitude and latitude (Note: 5 km x 5 km).
    
    Modification history
    --------------------
    Written (v.1.0): Michael Diamond, 08/04/16, Seattle, WA
    """
    
    def __init__(self,filename):
        self.filename= filename
        #Get geospatial information
        self.jday = int(filename[14:16+1]) #Julian day
        self.year = int(filename[10:14])
        self.month = cal_day(self.jday,self.year).split()[0]
        self.day = int(cal_day(self.jday,self.year).split()[1]) #Calendar day
        self.time = filename[18:19+1]+':'+filename[20:21+1]+' UTC'
        if filename[1] == 'O':
            self.satellite = 'Terra'
        elif filename[1] == 'Y':
            self.satellite = 'Aqua'
        h = SD.SD(self.filename, SDC.READ)
        self.lon = h.select('Longitude')[:,:]
        self.lat = h.select('Latitude')[:,:]
        SD.SD(self.filename).end()
        
    #Get radiances at different bands, account for scale and offsets
    def radiance(self,b, hi = False):
        """
        Get radiances from MODIS level 1b file.
        
        Parameters
        ----------
        filename : string
        Name of the .hdf level 1b file.
        
        b : int
        Band number, between 1-36.
        
        hi : boolean
        Hi- or lo-gain bands for bands 13 and 14. Set to True for hi-gain. Default is lo-gain.
        
        Returns
        -------
        data : masked array
        Array of radiance values in W/meters^2/micron/steradian. Invalid measurements are masked.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/4/16, Seattle, WA
        """
        if not type(b) == int:
            return 'Error: Band must be an integer between 1 and 36'
        if not 1 <= b <= 36:
            return 'Error: Band must be an integer between 1 and 36'
        h = SD.SD(self.filename, SDC.READ)
        if 1 <= b <= 2:
            ind = b-1
            rawdata = h.select('EV_250_Aggr1km_RefSB')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Apply scale and offset
            offset = attrs["radiance_offsets"][0][ind]
            scale = attrs["radiance_scales"][0][ind]
            data = (data - offset) * scale
            #Fix invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max, data < valid_min)
            data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
            return data
        elif 3 <= b <= 7:
            ind = b-3
            rawdata = h.select('EV_500_Aggr1km_RefSB')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Apply scale and offset
            offset = attrs["radiance_offsets"][0][ind]
            scale = attrs["radiance_scales"][0][ind]
            data = (data - offset) * scale
            #Fix invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            return data
        elif 8 <= b <= 19:
            if 8 <= b <= 12:
                ind = b-8
            elif b == 13:
                if hi == False:
                    ind = 5
                elif hi == True:
                    ind = 6
            elif b == 14:
                if hi == False:
                    ind = 7
                elif hi == True:
                    ind = 8
            elif 15 <= b <= 19:
                ind = b-6
            rawdata = h.select('EV_1KM_RefSB')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Apply offset and scale
            offset = attrs["radiance_offsets"][0][ind]
            scale = attrs["radiance_scales"][0][ind]
            data = (data - offset) * scale
            #Fix invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max, data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            return data
        elif 20 <= b <= 36 and b != 26:
            if 20 <= b <= 25:
                ind = b-20
            elif 27 <= b <= 36:
                ind = b-21
            rawdata = h.select('EV_1KM_Emissive')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Apply offset and scale
            offset = attrs["radiance_offsets"][0][ind]
            scale = attrs["radiance_scales"][0][ind]
            data = (data - offset) * scale
            #Fix invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max, data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            return data
        elif b == 26:
            ind = 0
            rawdata = h.select('EV_Band26')
            data = rawdata[:,:]
            attrs = rawdata.attributes(full=1)
            #Apply offset and scale
            offset = attrs["radiance_offsets"][0]
            scale = attrs["radiance_scales"][0]
            data = (data - offset) * scale
            #Fix invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max, data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            return data
        SD.SD(self.filename).end()

    #Get reflectances at different bands, account for scale and offsets
    def reflectance(self, b, hi = False):
        """
        Get reflectances from MODIS level 1b file.
        
        Parameters
        ----------
        filename : string
        Name of the .hdf level 1b file.
        
        b : int
        Band number, between 1-36.
        
        hi : boolean
        Hi- or lo-gain bands for bands 13 and 14. Set to True for hi-gain. Default is lo-gain.
        
        Returns
        -------
        data : masked array
        Array of relfectance values (unitless). Invalid measurements are masked.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/4/16, Seattle, WA
        """
        if not type(b) == int:
            return 'Error: Band must be an integer between 1 and 19 or 26'
        if not 1 <= b <= 36:
            return 'Error: Band must be an integer between 1 and 19 or 26'
        h = SD.SD(self.filename, SDC.READ)
        if 1 <= b <= 2:
            ind = b-1
            rawdata = h.select('EV_250_Aggr1km_RefSB')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Apply offset and scale
            offset = attrs["reflectance_offsets"][0][ind]
            scale = attrs["reflectance_scales"][0][ind]
            data = (data - offset) * scale
            #Fix invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            return data
        elif 3 <= b <= 7:
            ind = b-3
            rawdata = h.select('EV_500_Aggr1km_RefSB')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Apply offset and scale
            offset = attrs["reflectance_offsets"][0][ind]
            scale = attrs["reflectance_scales"][0][ind]
            data = (data - offset) * scale
            #Fix invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            return data
        elif 8 <= b <= 19:
            if 8 <= b <= 12:
                ind = b-8
            elif b == 13:
                if hi == False:
                    ind = 5
                elif hi == True:
                    ind = 6
            elif b == 14:
                if hi == False:
                    ind = 7
                elif hi == True:
                    ind = 8
            elif 15 <= b <= 19:
                ind = b-6
            rawdata = h.select('EV_1KM_RefSB')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Apply offset and scale
            offset = attrs["reflectance_offsets"][0][ind]
            scale = attrs["reflectance_scales"][0][ind]
            data = (data - offset) * scale
            #Fix invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            return data
        elif 20 <= b <= 36 and b != 26:
            return 'Error: Band must be an integer between 1 and 19 or 26'
        elif b == 26:
            ind = 0
            rawdata = h.select('EV_Band26')
            data = rawdata[:,:]
            attrs = rawdata.attributes(full=1)
            #Apply offset and scale
            offset = attrs["reflectance_offsets"][0]
            scale = attrs["reflectance_scales"][0]
            data = (data - offset) * scale
            #Fix invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            return data
        SD.SD(self.filename).end()

    #Calculate brightness temperature
    def Tb(self,band,hi=False):
        """
        Convert radiances to brightness temperatures. Same as brightness_temperature().
        
        Parameters
        ----------       
        band : int
        MODIS band, from 1 to 36.
        
        hi : boolean
        Hi- or lo-gain bands for bands 13 and 14. Default is lo-gain.
        
        Returns
        -------
        Tb : array
        Array of brightness temperature values in Kelvin.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/4/2016, Seattle, WA
        """
        h = 6.62607004.E-34 #m2 kg / s
        c = 299792458. #m/s
        k = 1.38064852.E-23 #m2 kg s-2 K-1
        B = self.radiance(band,hi)*10.**6 #In W/m2/m/sr
        if 1 <= band <= 19:
            wvl = avg_wavelength(band)/10.**9 #In m
        elif 20 <= band <= 36:
            wvl = avg_wavelength(band)/10.**6 #In m
        else:
            print 'Error: Band must be an integer between 1 and 36.'
        Tb = (h*c/(k*wvl)/(np.log(1+2.*h*c**2/(B*wvl**5)))) #In K
        return Tb
    
    def brightness_temperature(self,band,hi=False):
        """
        Convert radiances to brightness temperatures. Same as Tb().
        
        Parameters
        ----------       
        band : int
        MODIS band, from 1 to 36.
        
        hi : boolean
        Hi- or lo-gain bands for bands 13 and 14. Default is lo-gain.
        
        Returns
        -------
        Tb : array
        Array of brightness temperature values in Kelvin.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/4/2016, Seattle, WA
        """
        return self.Tb(band,hi)

    #Quickly plot data as check/first pass
    def quick_plot(self,band,hi=False,data='radiance',projection='merc'):
        """
        Quick plot data for a single MOD021KM or MYD021KM file. Intended as a check or first pass at data.
        
        Parameters
        ----------
        band : int
        MODIS band, from 1 to 36.
        
        hi : boolean
        Hi- or lo-gain bands for bands 13 and 14. Default is lo-gain.
        
        data : string
        Use 'rad' for radiances or 'ref' for reflectances. Default is radiance.
        
        projection : string
        """
        if not type(band) == int or band > 36:
            print 'Error: Band must be an integer from 1-36'
            return
        plt.figure()
        plt.clf()
        font = 'Arial'
        size = 20
        max_lon = self.lon.max()
        max_lat = self.lat.max()
        min_lon = self.lon.min()
        min_lat = self.lat.min()
        if projection == 'merc':
            m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
            urcrnrlat=max_lat+1,projection='merc',resolution='l')
        elif projection == 'global':
            m = Basemap(lon_0=0,projection='kav7',resolution='c')
        elif projection == 'satellite':
            m = Basemap(projection='nsper',lon_0=self.lon.mean(),lat_0=self.lat.mean(),\
            resolution='l',satellite_height=705000)
        else:
            print "Error: Projection must be 'merc', 'global', or 'satellite'."
            return
        m.drawcoastlines()
        m.drawcountries()
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        if data == 'radiance': 
            d = self.radiance(band,hi)[::5,::5]
            vmin = 0
            vmax = d.max()
        elif data == 'reflectance': 
            d = self.reflectance(band,hi)[::5,::5]
            vmin = 0
            vmax = 1
        elif data == 'brightness temperature' or data == 'Tb':
            d = self.Tb(band,hi)[::5,::5]
            vmin = d.min()
            vmax = d.max()
        m.pcolormesh(self.lon,self.lat,d[:np.shape(self.lon)[0],:np.shape(self.lat)[1]],\
        shading='gouraud',cmap=colorbar(band),latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=size-4) 
        label = {'radiance' : 'Radiance [W/m^2/micron/sr]', 'reflectance' : 'Reflectance [unitless]',\
        'brightness temperature' : 'Brightness temperature [K]', 'Tb' : 'Brightness temperature [K]'}
        cbar.set_label(label['%s' % data],fontname=font,fontsize=size-2)
        plt.title('Band %s %s for %s %s, %s from %s' % \
        (band,data,self.month,self.day,self.year,self.satellite), fontname=font,fontsize=size)

"""
Functions for M-D06_L2 files
"""

class cloud(object):
    """
    Calculate bias in MODIS effective radius data from M-D06_L2 cloud file.
    
    Parameters
    ----------
    cfile : string
    M-D06_L2 cloud file.
    
    Returns
    -------
    day, year, time, satellite : int, int, int, float
    Satellite properties.
    
    lat, lon : array, array
    Latitude and longitude for plotting.
    
    Re163 : array
    Effective radius calculated using the 1.63 micron channel.
    
    Re213 : array
    Effective radius calculated using the 2.13 micron channel.
    
    delta_Re : array
    Difference between the 2.13 micron and 1.63 micron channels.
    
    del_Re : array
    Difference between the 2.13 micron and 1.63 micron channels.
    
    R_Re : array
    Ratio of 2.13 micron over 1.63 micron channels.
    
    median : array
    Median effective radius bias.
    
    mean : array
    Mean effective radius bias.
    
    tau163 : array
    Cloud optical thickness calculated using the 1.63 micron channel.
    
    tau213 : array
    Cloud optical thickness calculated using the 2.13 micron channel.
    
    delta_tau : array
    Difference between the 2.13 micron and 1.63 micron channels.
    
    del_tau : array
    Difference between the 2.13 micron and 1.63 micron channels.
    
    R_tau : array
    Ratio of 2.13 micron over 1.63 micron channels.
    
    tau_median : array
    Median effective radius bias.
    
    tau_mean : array
    Mean effective radius bias.
    
    plot_Re : figure
    Figure showing the above arrays for the single file.
    
    plot_tau : figure
    Figure showing the above arrays for the single file.
    
    view_Re : figure
    View from satellite (effective radius bias).
    
    view_tau : figure
    View from satellite (cloud optical thickness bias).
    """
    
    def __init__(self,cfile):
        self.day = cfile[14:16+1]
        self.year = cfile[10:13+1]
    	self.time = cfile[18:21+1]
        if cfile[1] == 'Y':
            self.satellite = 'Aqua'
        elif cfile[1] == 'O':
            self.satellite = 'Terra'
        c = SD.SD(cfile, SDC.READ)
        #
        #1.63 micron band
        #
        self.Re163 = c.select('Cloud_Effective_Radius_16')[:]
        attrs_163 = c.select('Cloud_Effective_Radius_16').attributes(full=1)
        scale_163 = attrs_163['scale_factor'][0]
        offset_163 = attrs_163['add_offset'][0]
        self.tau163 = c.select('Cloud_Optical_Thickness_16')[:]
        attrs_tau163 = c.select('Cloud_Optical_Thickness_16').attributes(full=1)
        scale_tau163 = attrs_tau163['scale_factor'][0]
        offset_tau163 = attrs_tau163['add_offset'][0]
        #Fixer-upper
        self.Re163 = scale_163*(self.Re163 - offset_163)
        self.Re163 = ma.MaskedArray(self.Re163, mask = self.Re163 < 0, fill_value = 0)
        self.tau163 = scale_tau163*(self.tau163 - offset_tau163)
        self.tau163 = ma.MaskedArray(self.tau163, mask = self.tau163 < 0, fill_value = 0)
        #
        #2.13 micron band
        #
        self.Re213 = c.select('Cloud_Effective_Radius')[:]
        attrs_213 = c.select('Cloud_Effective_Radius').attributes(full=1)
        scale_213 = attrs_213['scale_factor'][0]
        offset_213 = attrs_213['add_offset'][0]
        self.tau213 = c.select('Cloud_Optical_Thickness')[:]
        attrs_tau213 = c.select('Cloud_Optical_Thickness').attributes(full=1)
        scale_tau213 = attrs_tau213['scale_factor'][0]
        offset_tau213 = attrs_tau213['add_offset'][0]
        #Fixer-upper
        self.Re213 = scale_213*(self.Re213 - offset_213)
        self.Re213 = ma.MaskedArray(self.Re213, mask = self.Re213 < 0, fill_value = 0)
        self.tau213 = scale_tau213*(self.tau213 - offset_tau213)
        self.tau213 = ma.MaskedArray(self.tau213, mask = self.tau213 < 0, fill_value = 0)
        
        #
        #3.7 micron band
        #
        self.Re37 = c.select('Cloud_Effective_Radius_37')[:]
        attrs_37 = c.select('Cloud_Effective_Radius_37').attributes(full=1)
        scale_37 = attrs_37['scale_factor'][0]
        offset_37 = attrs_37['add_offset'][0]
        self.tau37 = c.select('Cloud_Optical_Thickness_37')[:]
        attrs_tau37 = c.select('Cloud_Optical_Thickness_37').attributes(full=1)
        scale_tau37 = attrs_tau37['scale_factor'][0]
        offset_tau37 = attrs_tau37['add_offset'][0]
        #Fixer-upper
        self.Re37 = scale_37*(self.Re37 - offset_37)
        self.Re37 = ma.MaskedArray(self.Re37, mask = self.Re37 < 0, fill_value = 0)
        self.tau37 = scale_tau37*(self.tau37 - offset_tau37)
        self.tau37 = ma.MaskedArray(self.tau37, mask = self.tau37 < 0, fill_value = 0)
        
        #
        #Effective radius bias
        #
        self.delta_Re = self.Re213 - self.Re163
        self.del_Re = self.delta_Re/self.Re213
        self.R_Re = self.Re213/self.Re163
        #Create blocks with median values
        meds = np.zeros(np.shape(self.delta_Re))
        meds_del = np.zeros(np.shape(self.del_Re))
        meds_R = np.zeros(np.shape(self.R_Re))
        means = np.zeros(np.shape(self.delta_Re))
        for i in range(10):
            for j in range(10):
                #Pixels to use in this iteration
                pix = self.delta_Re[i*203:(i+1)*230,j*135:(j+1)*135]
                pix_del = self.del_Re[i*203:(i+1)*230,j*135:(j+1)*135]
                pix_R = self.R_Re[i*203:(i+1)*230,j*135:(j+1)*135]
                meds[i*203:(i+1)*230,j*135:(j+1)*135] = ma.median(pix)
                meds_del[i*203:(i+1)*230,j*135:(j+1)*135] = ma.median(pix_del)
                meds_R[i*203:(i+1)*230,j*135:(j+1)*135] = ma.median(pix_R)
                means[i*203:(i+1)*230,j*135:(j+1)*135] = ma.mean(pix)
        #Whatever's left
        finpix = self.delta_Re[10*203:,1350:]
        finpix_del = self.del_Re[10*203:,1350:]
        finpix_R = self.R_Re[10*203:,1350:]
        finlistpix = np.reshape(finpix,np.shape(finpix)[0]*np.shape(finpix)[1])
        meds[2030:,1350:] = ma.median(finlistpix)
        meds_del[2030:,1350:] = ma.median(finpix_del)
        meds_R[2030:,1350:] = ma.median(finpix_R)
        means[2030:,1350:] = ma.mean(finlistpix)
        meds = ma.MaskedArray(meds, mask = self.Re213 < 0, fill_value = 0)
        median = meds[::5,::5]
        meds_del = ma.MaskedArray(meds_del, mask = self.Re213 < 0, fill_value = 0)
        median_del = meds_del[::5,::5]
        meds_R = ma.MaskedArray(meds_R, mask = self.Re213 < 0, fill_value = 0)
        median_R = meds_R[::5,::5]
        means = ma.MaskedArray(means, mask = self.Re213 < 0, fill_value = 0)
        mean = means[::5,::5]
        self.median = median
        self.median_del = median_del
        self.median_R = median_R
        self.mean = mean
        #
        #Cloud optical thickness bias
        #
        self.delta_tau = self.tau213 - self.tau163
        self.del_tau = (self.tau213 - self.tau163)/self.tau213
        self.R_tau = self.tau213 / self.tau163
        #Create blocks with median values
        tau_meds = np.zeros(np.shape(self.delta_tau))
        tau_meds_del = np.zeros(np.shape(self.del_tau))
        tau_meds_R = np.zeros(np.shape(self.R_tau))
        tau_means = np.zeros(np.shape(self.delta_tau))
        for i in range(10):
            for j in range(10):
                #Pixels to use in this iteration
                pix = self.delta_tau[i*203:(i+1)*230,j*135:(j+1)*135]
                pix_del = self.del_tau[i*203:(i+1)*230,j*135:(j+1)*135]
                pix_R = self.R_tau[i*203:(i+1)*230,j*135:(j+1)*135]
                tau_meds[i*203:(i+1)*230,j*135:(j+1)*135] = ma.median(pix)
                tau_meds_del[i*203:(i+1)*230,j*135:(j+1)*135] = ma.median(pix_del)
                tau_meds_R[i*203:(i+1)*230,j*135:(j+1)*135] = ma.median(pix_R)
                tau_means[i*203:(i+1)*230,j*135:(j+1)*135] = ma.mean(pix)
        #Whatever's left
        finpix = self.delta_tau[10*203:,1350:]
        finpix_del = self.del_tau[10*203:,1350:]
        finpix_R = self.R_tau[10*203:,1350:]
        finlistpix = np.reshape(finpix,np.shape(finpix)[0]*np.shape(finpix)[1])
        tau_meds[2030:,1350:] = ma.median(finlistpix)
        tau_meds_del[2030:,1350:] = ma.median(finpix_del)
        tau_meds_R[2030:,1350:] = ma.median(finpix_R)
        tau_means[2030:,1350:] = ma.mean(finlistpix)
        tau_meds = ma.MaskedArray(tau_meds, mask = self.tau213 < 0, fill_value = 0)
        tau_median = tau_meds[::5,::5]
        tau_meds_del = ma.MaskedArray(tau_meds_del, mask = self.tau213 < 0, fill_value = 0)
        tau_median_del = tau_meds_del[::5,::5]
        tau_meds_R = ma.MaskedArray(tau_meds_R, mask = self.tau213 < 0, fill_value = 0)
        tau_median_R = tau_meds_R[::5,::5]
        tau_means = ma.MaskedArray(tau_means, mask = self.Re213 < 0, fill_value = 0)
        tau_mean = tau_means[::5,::5]
        self.tau_median = tau_median
        self.tau_median_del = tau_median_del
        self.tau_median_R = tau_median_R
        self.tau_mean = tau_mean
        #For plotting
        self.lat = c.select('Latitude')[:,:]
        self.lon = c.select('Longitude')[:,:]
        c.end()
    
    def plot_Re(self):
        """
        Make plot showing all the effective radius variables
        """
        plt.clf()
        lat = self.lat
        lon = self.lon
        Re_163 = self.Re163
        Re_213 = self.Re213
        delta_Re = self.delta_Re
        max_lon = lon.max()
        max_lat = lat.max()
        min_lon = lon.min()
        min_lat = lat.min()
        #Re 1.63 micron band
        plt.subplot(2,3,1)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        x, y = m(min_lon, max_lat-1)
        plt.text(x,y,'Julian day %s, %s' % (self.day, self.year),bbox=dict(facecolor='w', alpha=1))
        d163 = Re_163[::5,::5]
        m.pcolormesh(lon,lat,d163[:,:270],shading='gouraud',cmap='Spectral_r',latlon=True)
        cbar = m.colorbar()
        cbar.set_label('Effective radius (%sm)' % (u"\u03BC"))
        plt.title('A) Effective radius (1.63 %sm/860 nm)' % (u"\u03BC"))
        #Difference
        plt.subplot(2,3,2)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[0,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        x, y = m(min_lon, max_lat-1)
        plt.text(x,y,'Time: %s' % (self.time),bbox=dict(facecolor='w', alpha=1))
        dRe = delta_Re[::5,::5]
        m.pcolormesh(lon,lat,dRe[:,:270],shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=-7,vmax=7)
        cbar = m.colorbar()
        cbar.set_label('Difference (%sm)' % (u"\u03BC"))
        plt.title('B) Difference, R(1.63 %sm) - R(2.13 %sm)' % (u"\u03BC",u"\u03BC"))
        #Re 2.13 micron band
        plt.subplot(2,3,3)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[0,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        x, y = m(min_lon, max_lat-1)
        plt.text(x,y,'Satellite: %s' % (self.satellite),bbox=dict(facecolor='w', alpha=1))
        d213 = Re_213[::5,::5]
        m.pcolormesh(lon,lat,d213[:,:270],shading='gouraud',cmap='Spectral_r',latlon=True)
        cbar = m.colorbar()
        cbar.set_label('Effective radius (%sm)' % (u"\u03BC"))
        plt.title('C) Effective radius (2.13 %sm/860 nm)' % (u"\u03BC"))
        #Plot medians
        plt.subplot(2,2,3)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        m.pcolormesh(lon,lat,self.median[:np.shape(lat)[0],:np.shape(lat)[1]],shading='gouraud',cmap='YlGnBu_r',latlon=True,vmin=-5,vmax=0)
        cbar = m.colorbar()
        cbar.set_label('Bias (%sm)' % (u"\u03BC"))
        plt.title('D) Median bias')
        plt.show()
        #Plot means
        plt.subplot(2,2,4)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        m.pcolormesh(lon,lat,self.mean[:np.shape(lat)[0],:np.shape(lat)[1]],shading='gouraud',cmap='YlGnBu_r',latlon=True,vmin=-5,vmax=0)
        cbar = m.colorbar()
        cbar.set_label('Bias (%sm)' % (u"\u03BC"))
        plt.title('E) Mean bias')
        plt.show()
    
    def plot_tau(self):
        """
        Make plot showing all the cloud optical thickness variables
        """
        plt.clf()
        lat = self.lat
        lon = self.lon
        tau_163 = self.tau163
        tau_213 = self.tau213
        delta_tau = self.delta_tau
        max_lon = lon.max()
        max_lat = lat.max()
        min_lon = lon.min()
        min_lat = lat.min()
        #Tau 1.63 micron band
        plt.subplot(2,3,1)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        x, y = m(min_lon, max_lat-1)
        plt.text(x,y,'Julian day %s, %s' % (self.day, self.year),bbox=dict(facecolor='w', alpha=1))
        d163 = tau_163[::5,::5]
        m.pcolormesh(lon,lat,d163[:,:270],shading='gouraud',cmap='Spectral_r',latlon=True)
        cbar = m.colorbar()
        cbar.set_label('Cloud optical thickness')
        plt.title('A) 1.63 %sm/860 nm channels' % (u"\u03BC"))
        #Difference
        plt.subplot(2,3,2)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[0,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        x, y = m(min_lon, max_lat-1)
        plt.text(x,y,'Time: %s' % (self.time),bbox=dict(facecolor='w', alpha=1))
        dtau = delta_tau[::5,::5]
        m.pcolormesh(lon,lat,dtau[:,:270],shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=-1,vmax=1)
        cbar = m.colorbar()
        cbar.set_label('Difference')
        plt.title('B) Difference, 1.63 %sm - 2.13 %sm' % (u"\u03BC",u"\u03BC"))
        #Tau 2.13 micron band
        plt.subplot(2,3,3)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[0,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        x, y = m(min_lon, max_lat-1)
        plt.text(x,y,'Satellite: %s' % (self.satellite),bbox=dict(facecolor='w', alpha=1))
        d213 = tau_213[::5,::5]
        m.pcolormesh(lon,lat,d213[:,:270],shading='gouraud',cmap='Spectral_r',latlon=True)
        cbar = m.colorbar()
        cbar.set_label('Cloud optical thickness')
        plt.title('C) 2.13 %sm/860 nm channels' % (u"\u03BC"))
        #Plot medians
        plt.subplot(2,2,3)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        m.pcolormesh(lon,lat,self.tau_median[:np.shape(lat)[0],:np.shape(lat)[1]],shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=-.5,vmax=.5)
        cbar = m.colorbar()
        cbar.set_label('Bias')
        plt.title('D) Median bias')
        plt.show()
        #Plot means
        plt.subplot(2,2,4)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        m.pcolormesh(lon,lat,self.tau_mean[:np.shape(lat)[0],:np.shape(lat)[1]],shading='gouraud',cmap='RdYlBu_r',latlon=True,vmin=-.5,vmax=.5)
        cbar = m.colorbar()
        cbar.set_label('Bias')
        plt.title('E) Mean bias')
        plt.show()
            
    def view_Re(self):
        """
        Plot the view as seen from the satellite for Re bias
        """
        plt.clf()
        m = Basemap(projection='nsper',lon_0=self.lon.mean(),lat_0=self.lat.mean(),resolution='l',satellite_height=705000)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10))
        m.drawmeridians(np.arange(0,360,10))
        m.pcolormesh(self.lon,self.lat,self.median[:np.shape(self.lat)[0],:np.shape(self.lat)[1]],shading='gouraud',cmap='YlGnBu_r',latlon=True,vmin=-5,vmax=0)
        cbar = m.colorbar()
        cbar.set_label('Bias (%sm)' % (u"\u03BC"))
        plt.title('View from %s (%s, %s, at %s)' % \
        (self.satellite,cal_day(int(self.day),int(self.year)),self.year,self.time))
        plt.show()

