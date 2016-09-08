"""
Define simple functions to analyze MODIS data using Python.

Modification history
--------------------
Written (v.0.0): Michael Diamond, 8/3-4/2016, Seattle, WA
Modified (v.0.1): Michael Diamond, 8/8-9/2016, Seattle, WA
    -Fixed bugs in MOD021KM object
    -Created NRT MOD06 object for ORACLES campaign
"""

#Import libraries
from pyhdf import SD
from pyhdf.SD import SDC
import numpy as np
import numpy.ma as ma
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
from scipy.ndimage.interpolation import zoom
from matplotlib.colors import LogNorm

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

#Dictionary linking string month to numeric value
month_num = {'January' : 1, 'February' : 2, 'March' : 3, 'April' : 4, 'May' : 5, 'June' : 6,\
    'July' : 7, 'August' : 8, 'September' : 9, 'October' : 10, 'November' : 11, 'December' : 12}
#Dictionary linking string month number to name
month_name = {'1': 'January','2' : 'February','3' : 'March','4' : 'April','5' : 'May','6' : 'June',\
    '7' : 'July','8' : 'August','9' : 'September','10' : 'October','11' : 'November','12' : 'December'}

#Convert calendar day to Julian day
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
    String of default matplotlib cmap to use for a given band.
    
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
        return 'cubehelix_r'

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
    Written (v.1.0): Michael Diamond, 08/04/2016, Seattle, WA
    Modified (v.1.1): Michael Diamond, 08/08/2016, Seattle, WA
        -Fixed bug in brightness temperature
        -Changed aesthetics of plotting methods
        -Masks invalid data before applying scale and offset
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
        1 km x 1 km array of radiance values in W/meters^2/micron/steradian. Invalid measurements are masked.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/4/2016, Seattle, WA
        Modified: Michael Diamond, 8/8/2016, Seattle, WA
            -Mask invalid data before scaling
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
            #Mask invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max, data < valid_min)
            data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
            #Apply scale and offset
            offset = attrs["radiance_offsets"][0][ind]
            scale = attrs["radiance_scales"][0][ind]
            data = (data - offset) * scale
            return data
        elif 3 <= b <= 7:
            ind = b-3
            rawdata = h.select('EV_500_Aggr1km_RefSB')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Mask invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            #Apply scale and offset
            offset = attrs["radiance_offsets"][0][ind]
            scale = attrs["radiance_scales"][0][ind]
            data = (data - offset) * scale
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
            #Mask invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max, data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            #Apply offset and scale
            offset = attrs["radiance_offsets"][0][ind]
            scale = attrs["radiance_scales"][0][ind]
            data = (data - offset) * scale
            return data
        elif 20 <= b <= 36 and b != 26:
            if 20 <= b <= 25:
                ind = b-20
            elif 27 <= b <= 36:
                ind = b-21
            rawdata = h.select('EV_1KM_Emissive')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Mask invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max, data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            #Apply offset and scale
            offset = attrs["radiance_offsets"][0][ind]
            scale = attrs["radiance_scales"][0][ind]
            data = (data - offset) * scale
            return data
        elif b == 26:
            ind = 0
            rawdata = h.select('EV_Band26')
            data = rawdata[:,:]
            attrs = rawdata.attributes(full=1)
            #Mask invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max, data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            #Apply offset and scale
            offset = attrs["radiance_offsets"][0]
            scale = attrs["radiance_scales"][0]
            data = (data - offset) * scale
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
        1 km x 1 km array of relfectance values (unitless). Invalid measurements are masked.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/4/2016, Seattle, WA
        Modified: Michael Diamond, 8/8/2016, Seattle, WA
            -Mask invalid data before scaling
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
            #Mask invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            #Apply offset and scale
            offset = attrs["reflectance_offsets"][0][ind]
            scale = attrs["reflectance_scales"][0][ind]
            data = (data - offset) * scale
            return data
        elif 3 <= b <= 7:
            ind = b-3
            rawdata = h.select('EV_500_Aggr1km_RefSB')
            data = rawdata[ind,:,:]
            attrs = rawdata.attributes(full=1)
            #Mask invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            #Apply offset and scale
            offset = attrs["reflectance_offsets"][0][ind]
            scale = attrs["reflectance_scales"][0][ind]
            data = (data - offset) * scale
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
            #Mask invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            data = ma.MaskedArray(data,mask=invalid,fill_value=_FillValue)
            #Apply offset and scale
            offset = attrs["reflectance_offsets"][0][ind]
            scale = attrs["reflectance_scales"][0][ind]
            data = (data - offset) * scale
            return data
        elif 20 <= b <= 36 and b != 26:
            return 'Error: Band must be an integer between 1 and 19 or 26'
        elif b == 26:
            ind = 0
            rawdata = h.select('EV_Band26')
            data = rawdata[:,:]
            attrs = rawdata.attributes(full=1)
            #Mask invalid data
            valid_min = attrs["valid_range"][0][0]        
            valid_max = attrs["valid_range"][0][1]
            _FillValue = attrs["_FillValue"][0]
            invalid = np.logical_or(data > valid_max,data < valid_min)
            #Apply offset and scale
            offset = attrs["reflectance_offsets"][0]
            scale = attrs["reflectance_scales"][0]
            data = (data - offset) * scale
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
        1 km x 1 km array of brightness temperature values in Kelvin.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/4/2016, Seattle, WA
        Modified: Michael Diamond, 8/8/2016, Seattle, WA
            -Fixed syntax error in constants
        """
        h = 6.62607004*10**-34 #m2 kg / s
        c = 299792458. #m/s
        k = 1.38064852*10**-23 #m2 kg s-2 K-1
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
        1 km x 1 km array of brightness temperature values in Kelvin.
        
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
        Use 'merc' for mercator, 'global' for global plot (kav7), 'nsper' for Terra/Aqua's eye view.
        
        Returns
        -------
        Plot of data subsampled to 5 km x 5 km resolution.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/4/2016, Seattle, WA
        Modified: Michael Diamond, 8/8/2016, Seattle, WA
            -Aesthetic changes to plots
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
            m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
            m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        elif projection == 'global':
            m = Basemap(lon_0=0,projection='kav7',resolution='c')
            m.drawparallels(np.arange(-180,180,15),labels=[1,0,0,0])
            m.drawmeridians(np.arange(0,360,45),labels=[1,1,0,1])
        elif projection == 'satellite':
            m = Basemap(projection='nsper',lon_0=self.lon.mean(),lat_0=self.lat.mean(),\
            resolution='l',satellite_height=705000)
            m.bluemarble(alpha=.75)
        else:
            print "Error: Projection must be 'merc', 'global', or 'satellite'."
            return
        m.drawcoastlines()
        m.drawcountries()
        if data == 'radiance': 
            d = self.radiance(band,hi)[::5,::5]
            vmin = 0
            vmax = d.max()
            cmap = colorbar(band)
        elif data == 'reflectance': 
            d = self.reflectance(band,hi)[::5,::5]
            vmin = 0
            vmax = 1
            cmap = colorbar(band)+'_r'
        elif data == 'brightness temperature' or data == 'Tb':
            d = self.Tb(band,hi)[::5,::5]
            vmin = d.min()
            vmax = d.max()
            cmap = colorbar(band)
        m.pcolormesh(self.lon,self.lat,d[:np.shape(self.lon)[0],:np.shape(self.lat)[1]],\
        shading='gouraud',cmap=cmap,latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=size-4) 
        label = {'radiance' : 'Radiance [W/m^2/micron/sr]', 'reflectance' : 'Reflectance [unitless]',\
        'brightness temperature' : 'Brightness temperature [K]', 'Tb' : 'Brightness temperature [K]'}
        cbar.set_label(label['%s' % data],fontname=font,fontsize=size-2)
        plt.title('Band %s %s for %s %s, %s from %s' % \
        (band,data,self.month,self.day,self.year,self.satellite), fontname=font,fontsize=size)

"""
Functions for MOD06_L2/MYD06_L2 files downloaded from from LAADS Web (https://ladsweb.nascom.nasa.gov/).
"""

class MOD06(object):
    """
    Create object to analyze from MOD06_L2/MYD06_L2 cloud file downloaded from LAADS Web (https://ladsweb.nascom.nasa.gov/).
    
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

"""
Class for NRT MOD06_L2/MYD06_L2 files downloaded from from LANCE (https://lance.modaps.eosdis.nasa.gov/data_products/).
"""

class nrtMOD06(object):
    """
    Create object to analyze from MOD06_L2/MYD06_L2 cloud file downloaded from LANCE (https://lance.modaps.eosdis.nasa.gov/data_products/).
    
    ***Created for ORACLES campaign (https://espo.nasa.gov/home/oracles/)***
        
    Parameters
    ----------
    cfile : string
    M-D06_L2 cloud file.
    
    Methods
    -------
    get_ds: Get previously unaccessed dataset from file and save in dictionaries.
    
    quick_plot: Plot variables quickly as a check/first look at the data.
    
    tri_plot: Pre-defined groupings of three plots meant for use on the ORACLES field campaign.
    
    five_plot: Pre-defined groupings of five plots meant for use on the ORACLES field campaign.
        -value at each channel plus delta and del biases
    
    Returns
    -------
    jday, year, day : int
    Julian day, year, calendar day
    
    month, time, satellite : string
    Month, time of passage, satellite (Terra or Aqua)
    
    lat, lon : array, array
    5 km x 5 km latitude and longitude arrays.
    
    ds_name : dict
    Dictionary of named datasets available.
    
    ds : dict
    Dictionary of dataset arrays. Check ds_name for available parameters.
    
    units : dict
    Dictionary of units for each dataset array.
    
    Modification history
    --------------------
    Written: Michael Diamond, 08/08-09/2016, Seattle, WA
    Modified: Michael Diamond, 08/22/2016, Seattle, WA
    Modified: Michael Diamond, 08/29/2016, Swakopmund, Namibia
        -Collection 6 NRT processing now available!!!
        -Updates to make compatible with C6 (mostly to ref, also added COT differences)
    """
    
    def __init__(self,cfile):
        #Read in file
        self.file = cfile
        c = SD.SD(self.file, SDC.READ)
        #Dictionaries of all defined datasets in nrtMOD06 object
        ds_name = {} #Get full name from abbreviation
        ds = {} #Get dataset from abbreviation
        units = {} #To check units
        
        #
        ###Get geolocation data and date
        #
        self.jday = int(self.file[14:16+1])
        self.year = int(self.file[10:13+1])
    	self.month = cal_day(self.jday,self.year).split()[0]
    	self.day = int(cal_day(self.jday,self.year).split()[1])
    	self.time = self.file[18:21+1]
        if self.file[1] == 'Y':
            self.satellite = 'Aqua'
        elif self.file[1] == 'O':
            self.satellite = 'Terra'
        self.lon = c.select('Longitude')[:,:]
        ds_name['lon'] = 'Longitude'
        ds['lon'] = self.lon
        units['lon'] = 'degrees'
        self.lat = c.select('Latitude')[:,:]
        ds_name['lat'] = 'Latitude'
        ds['lat'] = self.lat
        units['lat'] = 'degrees'
        
        #
        ###Effective radius
        #
        data = c.select('Cloud_Effective_Radius')[:]
        attrs = c.select('Cloud_Effective_Radius').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.ref = scale*(data - offset)
        ds_name['ref'] = 'Cloud effective radius'
        ds['ref'] = self.ref
        units['ref'] = '%sm' % u"\u03BC"
        
        #
        ###Effective radius uncertainty
        #
        data = c.select('Cloud_Effective_Radius_Uncertainty')[:]
        attrs = c.select('Cloud_Effective_Radius_Uncertainty').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.ref_unc = scale*(data - offset)
        ds_name['ref_unc'] = 'Cloud effective radius uncertainty'
        ds['ref_unc'] = self.ref_unc
        units['ref_unc'] = '%'
        
        #
        ###Effective radius at 1.63 micron
        #
        data = c.select('Cloud_Effective_Radius_16')[:]
        attrs = c.select('Cloud_Effective_Radius_16').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.ref16 = scale*(data - offset)
        ds_name['ref16'] = 'Cloud effective radius (1.6 micron)'
        ds['ref16'] = self.ref16
        units['ref16'] = '%sm' % u"\u03BC"
        
        #
        ###Effective radius uncertainty at 1.63 micron
        #
        data = c.select('Cloud_Effective_Radius_Uncertainty_16')[:]
        attrs = c.select('Cloud_Effective_Radius_Uncertainty_16').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.ref16_unc = scale*(data - offset)
        ds_name['ref16_unc'] = 'Cloud effective radius uncertainty (1.6 micron)'
        ds['ref16_unc'] = self.ref16_unc
        units['ref16_unc'] = '%'
        
        #
        ###Effective radius at 3.7 micron
        #
        data = c.select('Cloud_Effective_Radius_37')[:]
        attrs = c.select('Cloud_Effective_Radius_37').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.ref37 = scale*(data - offset)
        ds_name['ref37'] = 'Cloud effective radius (3.7 micron)'
        ds['ref37'] = self.ref37
        units['ref37'] = '%sm' % u"\u03BC"
        
        #
        ###Effective radius uncertainty at 3.7 micron
        #
        data = c.select('Cloud_Effective_Radius_Uncertainty_37')[:]
        attrs = c.select('Cloud_Effective_Radius_Uncertainty_37').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.ref37_unc = scale*(data - offset)
        ds_name['ref37_unc'] = 'Cloud effective radius uncertainty (3.7 micron)'
        ds['ref37_unc'] = self.ref37_unc
        units['ref37_unc'] = '%'
        
        #
        ###Effective radius (1621)
        #
        data = c.select('Cloud_Effective_Radius_1621')[:]
        attrs = c.select('Cloud_Effective_Radius_1621').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.ref1621 = scale*(data - offset)
        ds_name['ref1621'] = 'Cloud effective radius (1621)'
        ds['ref1621'] = self.ref1621
        units['ref1621'] = '%sm' % u"\u03BC"
        
        #
        ###Effective radius uncertainty (1621)
        #
        data = c.select('Cloud_Effective_Radius_Uncertainty_1621')[:]
        attrs = c.select('Cloud_Effective_Radius_Uncertainty_1621').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.ref1621_unc = scale*(data - offset)
        ds_name['ref1621_unc'] = 'Cloud effective radius uncertainty (1621)'
        ds['ref1621_unc'] = self.ref1621_unc
        units['ref1621_unc'] = '%'
        
        #
        ###Cloud fraction
        #
        data = c.select('Cloud_Fraction')[:]
        attrs = c.select('Cloud_Fraction').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.cf = scale*(data - offset)
        ds_name['cf'] = 'Cloud fraction'
        ds['cf'] = self.cf
        units['cf'] = 'unitless'
        
        #
        ###Cloud optical thickness
        #
        data = c.select('Cloud_Optical_Thickness')[:]
        attrs = c.select('Cloud_Optical_Thickness').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.COT = scale*(data - offset)
        ds_name['COT'] = 'Cloud optical thickness'
        ds['COT'] = self.COT
        units['COT'] = 'unitless'
        
        #
        ###Cloud optical thickness uncertainty
        #
        data = c.select('Cloud_Optical_Thickness_Uncertainty')[:]
        attrs = c.select('Cloud_Optical_Thickness_Uncertainty').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.COT_unc = scale*(data - offset)
        ds_name['COT_unc'] = 'Cloud optical thickness uncertainty'
        ds['COT_unc'] = self.COT_unc
        units['COT_unc'] = '%'
        
        #
        ###Cloud optical thickness at 1.6 micron
        #
        data = c.select('Cloud_Optical_Thickness_16')[:]
        attrs = c.select('Cloud_Optical_Thickness_16').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.COT16 = scale*(data - offset)
        ds_name['COT16'] = 'Cloud optical thickness (1.6 micron)'
        ds['COT16'] = self.COT16
        units['COT16'] = 'unitless'
        
        #
        ###Cloud optical thickness uncertainty at 1.6 micron
        #
        data = c.select('Cloud_Optical_Thickness_Uncertainty_16')[:]
        attrs = c.select('Cloud_Optical_Thickness_Uncertainty_16').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.COT16_unc = scale*(data - offset)
        ds_name['COT16_unc'] = 'Cloud optical thickness uncertainty (1.6 micron)'
        ds['COT16_unc'] = self.COT16_unc
        units['COT16_unc'] = '%'
        
        #
        ###Cloud optical thickness at 3.7 micron
        #
        data = c.select('Cloud_Optical_Thickness_37')[:]
        attrs = c.select('Cloud_Optical_Thickness_37').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.COT37 = scale*(data - offset)
        ds_name['COT37'] = 'Cloud optical thickness (3.7 micron)'
        ds['COT37'] = self.COT37
        units['COT37'] = 'unitless'
        
        #
        ###Cloud optical thickness uncertainty at 3.7 micron
        #
        data = c.select('Cloud_Optical_Thickness_Uncertainty_37')[:]
        attrs = c.select('Cloud_Optical_Thickness_Uncertainty_37').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.COT37_unc = scale*(data - offset)
        ds_name['COT37_unc'] = 'Cloud optical thickness uncertainty (3.7 micron)'
        ds['COT37_unc'] = self.COT37_unc
        units['COT37_unc'] = '%'
        
        #
        ###Cloud optical thickness (1621)
        #
        data = c.select('Cloud_Optical_Thickness_1621')[:]
        attrs = c.select('Cloud_Optical_Thickness_1621').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.COT1621 = scale*(data - offset)
        ds_name['COT1621'] = 'Cloud optical thickness (1621)'
        ds['COT1621'] = self.COT1621
        units['COT1621'] = 'unitless'
        
        #
        ###Cloud optical thickness uncertainty (1621)
        #
        data = c.select('Cloud_Optical_Thickness_Uncertainty_1621')[:]
        attrs = c.select('Cloud_Optical_Thickness_Uncertainty_1621').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.COT1621_unc = scale*(data - offset)
        ds_name['COT1621_unc'] = 'Cloud optical thickness uncertainty (1621)'
        ds['COT1621_unc'] = self.COT1621_unc
        units['COT1621_unc'] = '%'
        
        #
        ###Cloud phase (optical properties)
        #
        data = c.select('Cloud_Phase_Optical_Properties')[:]
        attrs = c.select('Cloud_Phase_Optical_Properties').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.phase = scale*(data - offset)
        ds_name['phase'] = 'Cloud phase'
        ds['phase'] = self.phase
        units['phase'] = 'none'
        
        #
        ###Cloud top pressure
        #
        data = c.select('Cloud_Top_Pressure')[:]
        attrs = c.select('Cloud_Top_Pressure').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.CTP = scale*(data - offset)
        ds_name['CTP'] = 'Cloud top pressure'
        ds['CTP'] = self.CTP
        units['CTP'] = 'hPa'
        
        #
        ###Cloud top temperature
        #
        data = c.select('Cloud_Top_Temperature')[:]
        attrs = c.select('Cloud_Top_Temperature').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.CTT = scale*(data - offset)
        ds_name['CTT'] = 'Cloud top temperature'
        ds['CTT'] = self.CTT
        units['CTT'] = 'K'
        
        #
        ###Cloud water path
        #
        data = c.select('Cloud_Water_Path')[:]
        attrs = c.select('Cloud_Water_Path').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.CWP = scale*(data - offset)
        ds_name['CWP'] = 'Cloud water path'
        ds['CWP'] = self.CWP
        units['CWP'] = 'g/m^2'
        
        #
        ###Cloud water path uncertainty
        #
        data = c.select('Cloud_Water_Path_Uncertainty')[:]
        attrs = c.select('Cloud_Water_Path_Uncertainty').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.CWP_unc = scale*(data - offset)
        ds_name['CWP_unc'] = 'Cloud water path uncertainty'
        ds['CWP_unc'] = self.CWP_unc
        units['CWP_unc'] = '%'
        
        #
        ###Cloud water path (1621)
        #
        data = c.select('Cloud_Water_Path_1621')[:]
        attrs = c.select('Cloud_Water_Path_1621').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.CWP1621 = scale*(data - offset)
        ds_name['CWP1621'] = 'Cloud water path (1621)'
        ds['CWP1621'] = self.CWP1621
        units['CWP1621'] = 'g/m^2'
        
        #
        ###Cloud water path uncertainty (1621)
        #
        data = c.select('Cloud_Water_Path_Uncertainty_1621')[:]
        attrs = c.select('Cloud_Water_Path_Uncertainty_1621').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.CWP1621_unc = scale*(data - offset)
        ds_name['CWP1621_unc'] = 'Cloud water path uncertainty (1621)'
        ds['CWP1621_unc'] = self.CWP1621_unc
        units['CWP1621_unc'] = '%'
        
        #
        ###Difference between ref(2.1) - ref(1.6)
        #
        self.delta_ref16 = self.ref - self.ref16
        ds_name['delta_ref16'] = 'Effective radius difference (2.1-1.6 %sm)' % u"\u03BC"
        ds['delta_ref16'] = self.delta_ref16
        units['delta_ref16'] = '%sm' % u"\u03BC"
        #Normalized difference
        self.del_ref16 = 1000.*(self.delta_ref16)/self.ref
        ds_name['del_ref16'] = 'Normalized effective radius difference (2.1-1.6 %sm)' % u"\u03BC"
        ds['del_ref16'] = self.del_ref16
        units['del_ref16'] = 'per mil'
        
        #
        ###Difference between ref(2.1) - ref(3.7)
        #
        self.delta_ref37 = self.ref - self.ref37
        ds_name['delta_ref37'] = 'Effective radius difference (2.1-3.7 %sm)' % u"\u03BC"
        ds['delta_ref37'] = self.delta_ref37
        units['delta_ref37'] = '%sm' % u"\u03BC"
        #Normalized difference
        self.del_ref37 = 1000.*(self.delta_ref37)/self.ref
        ds_name['del_ref37'] = 'Normalized effective radius difference (2.1-3.7 %sm)' % u"\u03BC"
        ds['del_ref37'] = self.del_ref37
        units['del_ref37'] = 'per mil'
        
        #
        ###Difference between COT(2.1) - COT(1.6)
        #
        self.delta_COT16 = self.COT - self.COT16
        ds_name['delta_COT16'] = 'Cloud optical thickness difference (2.1-1.6 %sm)' % u"\u03BC"
        ds['delta_COT16'] = self.delta_COT16
        units['delta_COT16'] = 'unitless'
        #Normalized difference
        self.del_COT16 = 1000.*(self.delta_COT16)/self.COT
        ds_name['del_COT16'] = 'Normalized cloud optical thickness difference (2.1-1.6 %sm)' % u"\u03BC"
        ds['del_COT16'] = self.del_COT16
        units['del_COT16'] = 'per mil'
        
        #
        ###Difference between COT(2.1) - COT(3.7)
        #
        self.delta_COT37 = self.COT - self.COT37
        ds_name['delta_COT37'] = 'Cloud optical thickness difference (2.1-3.7 %sm)' % u"\u03BC"
        ds['delta_COT37'] = self.delta_COT37
        units['delta_COT37'] = 'unitless'
        #Normalized difference
        self.del_COT37 = 1000.*(self.delta_COT37)/self.COT
        ds_name['del_COT37'] = 'Normalized cloud optical thickness difference (2.1-3.7 %sm)' % u"\u03BC"
        ds['del_COT37'] = self.del_COT37
        units['del_COT37'] = 'per mil'
        
        #
        ###Sensor/view azimuth angle
        #
        data = c.select('Sensor_Azimuth')[:]
        attrs = c.select('Sensor_Azimuth').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.vaa = scale*(data - offset)
        ds_name['vaa'] = 'Sensor azimuth angle'
        ds['vaa'] = self.vaa
        units['vaa'] = 'degrees'
        
        #
        ###Sensor/view zenith angle
        #
        data = c.select('Sensor_Zenith')[:]
        attrs = c.select('Sensor_Zenith').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.vza = scale*(data - offset)
        ds_name['vza'] = 'Sensor zenith angle'
        ds['vza'] = self.vza
        units['vza'] = 'degrees'
        
        #
        ###Solar azimuth angle
        #
        data = c.select('Solar_Azimuth')[:]
        attrs = c.select('Solar_Azimuth').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.saa = scale*(data - offset)
        ds_name['saa'] = 'Solar azimuth angle'
        ds['saa'] = self.saa
        units['saa'] = 'degrees'
        
        #
        ###Solar zenith angle
        #
        data = c.select('Solar_Zenith')[:]
        attrs = c.select('Solar_Zenith').attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        self.sza = scale*(data - offset)
        ds_name['sza'] = 'Solar zenith angle'
        ds['sza'] = self.sza
        units['sza'] = 'degrees'
        
        #
        ###Relative azimuth angle
        #
        self.raa = self.vaa - self.saa
        ds_name['raa'] = 'Relative azimuth angle'
        ds['raa'] = self.raa
        units['raa'] = 'degrees'
        
        #
        ###Nd at 2.13 micron
        #
        k = .8
        gam_ad = 2.E-6 #kg/m^4
        frac_ad = 1
        gam_eff = gam_ad*frac_ad
        self.Nd = 10**.5/(4*np.pi*1000**.5)*gam_eff**.5*self.COT**.5/(self.ref*10**-6)**2.5/k/100.**3
        ds_name['Nd'] = 'Droplet concentration'
        ds['Nd'] = self.Nd
        units['Nd'] = '$\mathregular{cm^{-3}}$'
        
        #
        ###Nd at 1.6 micron
        #
        k = .8
        gam_ad = 2.E-6 #kg/m^4
        frac_ad = 1
        gam_eff = gam_ad*frac_ad
        self.Nd16 = 10**.5/(4*np.pi*1000**.5)*gam_eff**.5*self.COT16**.5/(self.ref16*10**-6)**2.5/k/100.**3
        ds_name['Nd16'] = 'Droplet concentration at 1.6 micron'
        ds['Nd16'] = self.Nd16
        units['Nd16'] = '$\mathregular{cm^{-3}}$'
        
        #
        ###Nd at 3.7 micron
        #
        k = .8
        gam_ad = 2.E-6 #kg/m^4
        frac_ad = 1
        gam_eff = gam_ad*frac_ad
        self.Nd37 = 10**.5/(4*np.pi*1000**.5)*gam_eff**.5*self.COT37**.5/(self.ref37*10**-6)**2.5/k/100.**3
        ds_name['Nd37'] = 'Droplet concentration at 3.7 micron'
        ds['Nd37'] = self.Nd37
        units['Nd37'] = '$\mathregular{cm^{-3}}$'
        
        #
        ###Difference between Nd(2.1) - Nd(1.6)
        #
        self.delta_Nd16 = self.Nd - self.Nd16
        ds_name['delta_Nd16'] = 'Cloud optical thickness difference (2.1-1.6 %sm)' % u"\u03BC"
        ds['delta_Nd16'] = self.delta_Nd16
        units['delta_Nd16'] = 'unitless'
        #Normalized difference
        self.del_Nd16 = 1000.*(self.delta_Nd16)/self.Nd
        ds_name['del_Nd16'] = 'Normalized cloud optical thickness difference (2.1-1.6 %sm)' % u"\u03BC"
        ds['del_Nd16'] = self.del_Nd16
        units['del_Nd16'] = 'per mil'
        
        #
        ###Difference between Nd(2.1) - Nd(3.7)
        #
        self.delta_Nd37 = self.Nd - self.Nd37
        ds_name['delta_Nd37'] = 'Cloud optical thickness difference (2.1-3.7 %sm)' % u"\u03BC"
        ds['delta_Nd37'] = self.delta_Nd37
        units['delta_Nd37'] = 'unitless'
        #Normalized difference
        self.del_Nd37 = 1000.*(self.delta_Nd37)/self.Nd
        ds_name['del_Nd37'] = 'Normalized cloud optical thickness difference (2.1-3.7 %sm)' % u"\u03BC"
        ds['del_Nd37'] = self.del_Nd37
        units['del_Nd37'] = 'per mil'
        
        #Finishing touches
        self.ds_name = ds_name
        self.ds = ds
        self.units = units
        c.end()
    
    def get_ds(self,dataset):
        """
        Get dataset from file that isn't already included.
        
        Parameters
        ----------
        dataset : string
        Name of dataset. 
        """
        c = SD.SD(self.file, SDC.READ)
        data = c.select('%s' % dataset)[:]
        attrs = c.select('%s' % dataset).attributes(full=1)
        scale = attrs['scale_factor'][0]
        offset = attrs['add_offset'][0]
        unit = attrs['units'][0]
        #Mask invalid data
        valid_min = attrs["valid_range"][0][0]        
        valid_max = attrs["valid_range"][0][1]
        _FillValue = attrs["_FillValue"][0]
        invalid = np.logical_or(data > valid_max, data < valid_min)
        data = ma.MaskedArray(data, mask=invalid, fill_value=_FillValue)
        #Fixer-upper
        data = scale*(data - offset)
        c.end()
        #Add to dictionaries for future use
        self.ds_name['%s' % dataset] = dataset
        self.ds['%s' % dataset] = data
        self.units['%s' % dataset] = unit
        return data
    
    #Quickly plot data as check/first pass
    def quick_plot(self,data='cf',projection='merc'):
        """
        Quick plot data for a single MOD06 or MYD06 NRT file. Intended as a check or first pass at data.
        
        Parameters
        ----------
        data : string
        Lookup key for "ds" dictionary. See dictionary for choices.
        
        projection : string
        Use 'merc' for mercator, 'global' for global plot (kav7), 'nsper' for Terra/Aqua's eye view.
        
        Returns
        -------
        Plot of data subsampled to 5 km x 5 km resolution.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/9/2016, Seattle, WA
        """
        if not type(data) == str:
            print 'Error: Data must be a valid key to the "ds" dictionary.'
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
            m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
            m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        elif projection == 'global':
            m = Basemap(lon_0=0,projection='kav7',resolution='c')
            m.drawparallels(np.arange(-180,180,15),labels=[1,0,0,0])
            m.drawmeridians(np.arange(0,360,45),labels=[1,1,0,1])
        elif projection == 'satellite':
            m = Basemap(projection='nsper',lon_0=self.lon.mean(),lat_0=self.lat.mean(),\
            resolution='l',satellite_height=705000)
            m.bluemarble(alpha=.75)
        else:
            print "Error: Projection must be 'merc', 'global', or 'satellite'."
            return
        m.drawcoastlines()
        m.drawcountries()
        d = self.ds['%s' % data]
        if np.shape(d)[0] != np.shape(self.lat)[0]: d = d[::5,::5]
        #Might make these customized in future, for now generic
        cmap = 'viridis'
        m.pcolormesh(self.lon,self.lat,d[:np.shape(self.lon)[0],:np.shape(self.lat)[1]],\
        shading='gouraud',cmap=cmap,latlon=True)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=size-4) 
        cbar.set_label('[%s]' % self.units['%s' % data],fontname=font,fontsize=size-2)
        plt.title('%s for %s %s, %s from %s' % \
        (self.ds_name['%s' % data],self.month,self.day,self.year,self.satellite), fontname=font,fontsize=size)
    
    def triplot(self,data='ref',projection='merc',full_res=False,num=None):
        """
        Plot three datasets for a single MOD06 or MYD06 NRT file.
        
        Parameters
        ----------
        data : string
        Choice of...
            'cot' for cloud optical thickness, cloud top pressure, and cloud water path
            'geo' for solar zenith angle, relative azimuth angle, and sensor zenith angle
            'ref' for standard effective radius (2.1 micron) and biases (1.6 micron)
            'Nd' for ref, Nd, and COT
            
        projection : string
        Use 'merc' for mercator. Other options possible in future update.
        
        full_res : boolean
        If True, plot at 1 km x 1 km resolution, interpolating as necessary. If False, plot at 5 km x 5 km resolution by subsampling if necessary.
        
        num : int or string
        Default None. If an integer or a string, calls that figure number/name to make plot.
 
        Returns
        -------
        Figure with three subplots.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/10-11/2016, Seattle, WA
        """
        plt.figure(num=num,figsize=(13.33,7.5))
        font = 'Arial'
        size = 16
        if full_res:
            lat = zoom(self.lat,5.)
            lon = zoom(self.lon,5.)
        else:
            lat = self.lat
            lon = self.lon
        max_lon = lon.max()
        max_lat = lat.max()
        min_lon = lon.min()
        min_lat = lat.min()
        
        plt.subplot(1,3,1)
        if projection == 'merc':
            m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
            urcrnrlat=max_lat+1,projection='merc',resolution='l')
            m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
            m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
            x, y = m(min_lon, max_lat-1)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        plt.text(x,y,'%s %s, %s' % (self.month, self.day, self.year),\
        bbox=dict(facecolor='w', alpha=1),fontname=font,fontsize=size-4)
        if data == 'ref': 
            key = 'delta_ref16'
            t = '%sref (1.6 %sm)' % (u"\u0394",u"\u03BC")
            c = 'RdYlBu_r'
            vmin = -6
            vmed = 0
            vmax = 6
        elif data == 'geo': 
            key = 'sza'
            t = 'Solar zenith angle'
            c = 'magma_r'
            vmin = 0.
            vmed = 45.
            vmax = 90.
        elif data == 'cot': 
            key = 'COT'
            t = 'Cloud optical thickness'
            c = 'viridis'
            vmin = 0
            vmed = 16
            vmax = 32
        elif data == 'Nd': 
            key = 'ref'
            t = 'Effective radius'
            c = 'viridis'
            vmin = 4
            vmed = 14
            vmax = 24
        d = self.ds['%s' % key]
        if not full_res and data != 'geo': d = d[::5,::5]
        if full_res and data == 'geo': d = zoom(d,5.)
        m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
        cmap=c,latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax],\
        location='bottom',pad=.35)
        cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
        cbar.ax.tick_params(labelsize=size-4)
        cbar.set_label('[%s]' % self.units['%s' % key],fontname=font,fontsize=size-2)
        plt.title('A) %s' % t,fontname=font,fontsize=size)
        
        plt.subplot(1,3,2)
        if projection == 'merc':
            m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
            urcrnrlat=max_lat+1,projection='merc',resolution='l')
            m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
            m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        plt.text(x,y,'Time: %s UTC' % (self.time),bbox=dict(facecolor='w', alpha=1),\
        fontname=font,fontsize=size-4)
        d = self.ref
        if data == 'ref': 
            key = 'ref'
            t = 'Effective radius (2.1 %sm)' % u"\u03BC"
            c = 'viridis'
            vmin = 4
            vmed = 14
            vmax = 24
        elif data == 'geo': 
            key = 'raa'
            t = 'Relative azimuth angle'
            c = 'RdYlBu_r'
            vmin = -180
            vmed = 0
            vmax = 180
        elif data == 'cot': 
            key = 'CTP'
            t = 'Cloud top pressure'
            c = 'plasma_r'
            vmin = 500
            vmed = 750
            vmax = 1000
        elif data == 'Nd': 
            key = 'Nd'
            t = 'Droplet concentration'
            c = 'viridis'
            vmin = 0
            vmed = 500
            vmax = 1000
        d = self.ds['%s' % key]
        if not full_res and data != 'geo' and data != 'cot': d = d[::5,::5]
        if not full_res and data == 'geo': d = d*self.ref[::5,::5][:np.shape(d)[0],:np.shape(d)[1]]/self.ref[::5,::5][:np.shape(d)[0],:np.shape(d)[1]] #Cloud mask?
        if full_res and data == 'geo': 
            d = zoom(d,5.)*self.ref[:np.shape(lon)[0],:np.shape(lat)[1]]/self.ref[:np.shape(lon)[0],:np.shape(lat)[1]] #Kludge
        if full_res and data == 'cot': 
            d = zoom(d,5.)*self.ref[:np.shape(lon)[0],:np.shape(lat)[1]]/self.ref[:np.shape(lon)[0],:np.shape(lat)[1]] #Kludge
        m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
        cmap=c,latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax],\
        location='bottom',pad=.35)
        cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
        cbar.ax.tick_params(labelsize=size-4)
        cbar.set_label('[%s]' % self.units['%s' % key],fontname=font,fontsize=size-2)
        plt.title('B) %s' % t,fontname=font,fontsize=size)
        
        plt.subplot(1,3,3)
        if projection == 'merc':
            m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
            urcrnrlat=max_lat+1,projection='merc',resolution='l')
            m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
            m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        plt.text(x,y,'Satellite: %s' % (self.satellite),bbox=dict(facecolor='w', alpha=1),\
        fontname=font,fontsize=size-4)
        d = self.ref
        if data == 'ref': 
            key = 'del_ref16'
            t = '%sref (1.6 %sm)' % (u"\u03B4",u"\u03BC")
            c = 'RdYlBu_r'
            vmin = -500
            vmed = 0
            vmax = 500
        elif data == 'geo': 
            key = 'vza'
            t = 'Sensor zenith angle'
            c = 'magma_r'
            vmin = 0.
            vmed = 45.
            vmax = 90.
        elif data == 'cot': 
            key = 'CWP'
            t = 'Cloud water path'
            c = 'viridis'
            vmin = 0
            vmed = 150
            vmax = 300
        elif data == 'Nd': 
            key = 'COT'
            t = 'Cloud optical thickness'
            c = 'viridis'
            vmin = 0
            vmed = 16
            vmax = 32
        d = self.ds['%s' % key]
        if not full_res and data != 'geo': d = d[::5,::5]
        if full_res and data == 'geo': d = zoom(d,5.)
        m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
        cmap=c,latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax],\
        location='bottom',pad=.35)
        cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
        cbar.ax.tick_params(labelsize=size-4)
        cbar.set_label('[%s]' % self.units['%s' % key],fontname=font,fontsize=size-2)
        plt.title('C) %s' % t,fontname=font,fontsize=size)
        
        plt.show()

    def five_plot(self,data='ref'):
        """
        Make plot showing variable at each wavelength and biases
        
        Parameters
        ----------
        data : string
        Choice of 'cot','ref', or 'Nd'. Default is effective radius.
        
        Modification history
        --------------------
        Written: Michael Diamond, 08/29/2016, Swakopmund, Namibia
        """
        plt.clf()
        font = 'Arial'
        size = 15
        if data == 'ref':
            dA = self.ref16
            dB = self.ref
            dC = self.ref37
            dD = self.delta_ref16
            dE = self.del_ref16
            vABC = (4,14,24)
            vD = (-6,0,6)
            vE = (-500,0,500)
            kABC = 'ref'
            kD = 'delta_ref16'
            kE = 'del_ref16'
        elif data == 'cot':
            dA = self.COT16
            dB = self.COT
            dC = self.COT37
            dD = self.delta_COT16
            dE = self.del_COT16
            vABC = (0,16,32)
            vD = (-1,0,1)
            vE = (-100,0,100)
            kABC = 'COT'
            kD = 'delta_COT16'
            kE = 'del_COT16'
        elif data == 'Nd':
            dA = self.Nd16
            dB = self.Nd
            dC = self.Nd37
            dD = self.delta_Nd16
            dE = self.del_Nd16
            vABC = (0,500,1000)
            vD = (-300,0,300)
            vE = (-1000,0,1000)
            kABC = 'Nd'
            kD = 'delta_Nd16'
            kE = 'del_Nd16'
        else:
            print 'Error: Invalid data input.'
            return
        #
        ###Variables at different wavelengths
        #
        key = kABC
        vmin = vABC[0]
        vmed = vABC[1]
        vmax = vABC[-1]
        lat = self.lat
        lon = self.lon
        max_lon = lon.max()
        max_lat = lat.max()
        min_lon = lon.min()
        min_lat = lat.min()
        #1.63 micron band
        plt.subplot(2,3,1)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        x, y = m(min_lon, max_lat-1)
        plt.text(x,y,'%s %s, %s' % (self.month,self.day, self.year),bbox=dict(facecolor='w', alpha=1),\
        fontname=font,fontsize=size-2)
        d = dA[::5,::5]
        if data == 'Nd': 
            m.drawmapboundary(fill_color='steelblue')
            m.fillcontinents(color='floralwhite',lake_color='steelblue',zorder=0)
            vmin = 1
            vmax = 1000
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='cubehelix',latlon=True,norm = LogNorm(vmin=vmin, vmax=vmax))
            cbar = m.colorbar(ticks=[1,10,100,1000])
            cbar.ax.set_xticklabels([1,10,100,1000])
            cbar.ax.tick_params(labelsize=size-4)
            cbar.set_label('[per cc]',fontname=font,fontsize=size-2)
        else:
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='viridis',latlon=True,vmin=vmin,vmax=vmax)
            cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
            cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
            cbar.ax.tick_params(labelsize=size-4)
            cbar.set_label('[%s]' % self.units['%s' % key],fontname=font,fontsize=size-2)
        plt.title('A) 1.63 %sm/860 nm channels' % (u"\u03BC"),fontname=font,fontsize=size)
        #2.13 micron band
        plt.subplot(2,3,2)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[0,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        x, y = m(min_lon, max_lat-1)
        plt.text(x,y,'Time: %s' % (self.time),bbox=dict(facecolor='w', alpha=1),\
        fontname=font,fontsize=size-2)
        d = dB[::5,::5]
        if data == 'Nd': 
            m.drawmapboundary(fill_color='steelblue')
            m.fillcontinents(color='floralwhite',lake_color='steelblue',zorder=0)
            vmin = 1
            vmax = 1000
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='cubehelix',latlon=True,norm = LogNorm(vmin=vmin, vmax=vmax))
            cbar = m.colorbar(ticks=[1,10,100,1000])
            cbar.ax.set_xticklabels([1,10,100,1000])
            cbar.ax.tick_params(labelsize=size-4)
            cbar.set_label('[per cc]',fontname=font,fontsize=size-2)
        else:
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='viridis',latlon=True,vmin=vmin,vmax=vmax)
            cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
            cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
            cbar.ax.tick_params(labelsize=size-4)
            cbar.set_label('[%s]' % self.units['%s' % key],fontname=font,fontsize=size-2)
        plt.title('B) Standard %s retrieval' % (key),fontname=font,fontsize=size)
        #3.7 micron band
        plt.subplot(2,3,3)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[0,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        x, y = m(min_lon, max_lat-1)
        plt.text(x,y,'Satellite: %s' % (self.satellite),bbox=dict(facecolor='w', alpha=1),\
        fontname=font,fontsize=size-2)
        d = dC[::5,::5]
        if data == 'Nd': 
            m.drawmapboundary(fill_color='steelblue')
            m.fillcontinents(color='floralwhite',lake_color='steelblue',zorder=0)
            vmin = 1
            vmax = 1000
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='cubehelix',latlon=True,norm = LogNorm(vmin=vmin, vmax=vmax))
            cbar = m.colorbar(ticks=[1,10,100,1000])
            cbar.ax.set_xticklabels([1,10,100,1000])
            cbar.ax.tick_params(labelsize=size-4)
            cbar.set_label('[per cc]',fontname=font,fontsize=size-2)
        else:
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='viridis',latlon=True,vmin=vmin,vmax=vmax)
            cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
            cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
            cbar.ax.tick_params(labelsize=size-4)
            cbar.set_label('[%s]' % self.units['%s' % key],fontname=font,fontsize=size-2)
        plt.title('C) 3.7 %sm/860 nm channels' % (u"\u03BC"),fontname=font,fontsize=size)
        #
        ###Bias plots
        #
        lon = zoom(lon,5.)
        lat = zoom(lat,5.)
        #Delta
        key = kD
        vmin = vD[0]
        vmed = vD[1]
        vmax = vD[-1]
        plt.subplot(2,2,3)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        d = dD
        if data == 'Nd':
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='RdYlBu',latlon=True,vmin=vmin,vmax=vmax)
        else:
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='RdYlBu_r',latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
        cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
        cbar.ax.tick_params(labelsize=size-4)
        cbar.set_label('[%s]' % self.units['%s' % key],fontname=font,fontsize=size-2)
        plt.title('D) %s (2.13 %sm - 1.63 %sm)' % (u"\u0394",u"\u03BC",u"\u03BC"),fontname=font,fontsize=size)
        #Del
        key = kE
        vmin = vE[0]
        vmed = vE[1]
        vmax = vE[-1]
        plt.subplot(2,2,4)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        d = dE
        if data == 'Nd':
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='RdYlBu',latlon=True,vmin=vmin,vmax=vmax)
        else:
            m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
            cmap='RdYlBu_r',latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
        cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax])
        cbar.ax.tick_params(labelsize=size-4)
        cbar.set_label('[%s]' % self.units['%s' % key],fontname=font,fontsize=size-2)
        plt.title('E) %s (2.13 %sm - 1.63 %sm)' % (u"\u03B4",u"\u03BC",u"\u03BC"),fontname=font,fontsize=size)
        plt.show()

"""
Class for NRT MOD06ACAERO/MYD06ACAERO files downloaded from from LANCE (https://lance.modaps.eosdis.nasa.gov/data_products/).
"""

class nrtACAERO(object):
    """
    Create object to analyze from MOD06ACAERO/MYD06ACAERO cloud file downloaded from LANCE (https://lance.modaps.eosdis.nasa.gov/data_products/).
    
    ***Created for ORACLES campaign (https://espo.nasa.gov/home/oracles/)***
    
    ***Still research-level product***
        
    Parameters
    ----------
    afile : string
    M-D06ACAERO cloud file.
    
    Methods
    -------
    get_ds: Get previously unaccessed dataset from file and save in dictionaries.
    
    quick_plot: Plot variables quickly as a check/first look at the data.
    
    tri_plot: Pre-defined groupings of three plots meant for use on the ORACLES field campaign.
    
    Returns
    -------
    jday, year, day : int
    Julian day, year, calendar day
    
    month, time, satellite : string
    Month, time of passage, satellite (Terra or Aqua)
    
    lat, lon : array, array
    5 km x 5 km latitude and longitude arrays.
    
    ds_name : dict
    Dictionary of named datasets available.
    
    ds : dict
    Dictionary of dataset arrays. Check ds_name for available parameters.
    
    Modification history
    --------------------
    Written: Michael Diamond, 08/30/2016, Swakopmund, Namibia
    """
    
    def __init__(self,cfile):
        #Read in file
        self.file = cfile
        a = SD.SD(self.file, SDC.READ)
        #Dictionaries of all defined datasets
        ds_name = {} #Get full name from abbreviation
        ds = {} #Get dataset from abbreviation
        
        #
        ###Get geolocation data and date
        #
        self.jday = int(self.file[17:20])
        self.year = int(self.file[13:17])
    	self.month = cal_day(self.jday,self.year).split()[0]
    	self.day = int(cal_day(self.jday,self.year).split()[1])
    	self.time = self.file[21:25]
        if self.file[1] == 'Y':
            self.satellite = 'Aqua'
        elif self.file[1] == 'O':
            self.satellite = 'Terra'
        self.lon = a.select('Longitude')[:,:]
        ds_name['lon'] = 'Longitude'
        ds['lon'] = self.lon
        self.lat = a.select('Latitude')[:,:]
        ds_name['lat'] = 'Latitude'
        ds['lat'] = self.lat
        
        datasets = ['Above_Cloud_AOD','Above_Cloud_AOD_ModAbsAero','Clear_Sky_AOD',\
        'Cloud_Effective_Radius','Cloud_Effective_Radius_ModAbsAero',\
        'Cloud_Optical_Thickness','Cloud_Optical_Thickness_ModAbsAero',\
        'Above_Cloud_AOD_Uncertainty','Above_Cloud_AOD_ModAbsAero_Uncertainty']
        
        for dset in datasets:
            data = a.select('%s' % dset)[:]
            attrs = a.select('%s' % dset).attributes(full=1)
            scale = attrs['scale_factor'][0]
            offset = attrs['add_offset'][0]
            #_FillValue = attrs["_FillValue"][0] #Why needed? No max/min valid
            ds_name['%s' % dset] = dset
            ds['%s' % dset] = scale*(data - offset)
        
        #
        ###Nd
        #
        k = .8
        gam_ad = 2.E-6 #kg/m^4
        frac_ad = 1
        gam_eff = gam_ad*frac_ad
        ds['Nd'] = 10**.5/(4*np.pi*1000**.5)*gam_eff**.5*ds['Cloud_Optical_Thickness']**.5/(ds['Cloud_Effective_Radius']*10**-6)**2.5/k/100.**3
        ds_name['Nd'] = 'Nd'
        ds['Nd_ModAbsAero'] = 10**.5/(4*np.pi*1000**.5)*gam_eff**.5*ds['Cloud_Optical_Thickness_ModAbsAero']**.5/(ds['Cloud_Effective_Radius_ModAbsAero']*10**-6)**2.5/k/100.**3
        ds_name['Nd_ModAbsAero'] = 'Nd_ModAbsAero'       
        
        #
        ###Masking invalid AOD data
        #
        ds['Above_Cloud_AOD'] = ma.MaskedArray(ds['Above_Cloud_AOD'],ds['Above_Cloud_AOD_Uncertainty']>100,fill_value=0).filled()
        ds['Above_Cloud_AOD'] = ma.MaskedArray(ds['Above_Cloud_AOD'],ds['Cloud_Optical_Thickness']<4,fill_value=0)
        ds['Above_Cloud_AOD_ModAbsAero'] = ma.MaskedArray(ds['Above_Cloud_AOD_ModAbsAero'],ds['Above_Cloud_AOD_ModAbsAero_Uncertainty']>100,fill_value=0).filled()
        ds['Above_Cloud_AOD_ModAbsAero'] = ma.MaskedArray(ds['Above_Cloud_AOD_ModAbsAero'],ds['Cloud_Optical_Thickness_ModAbsAero']<4,fill_value=0)
        
        #Finishing touches
        self.ds_name = ds_name
        self.ds = ds
        a.end()
    
    #Quickly plot data as check/first pass
    def quick_plot(self,data='Above_Cloud_AOD',projection='merc'):
        """
        Quick plot data for a single MOD06 or MYD06 NRT file. Intended as a check or first pass at data.
        
        Parameters
        ----------
        data : string
        Lookup key for "ds" dictionary. See dictionary for choices.
        
        projection : string
        Use 'merc' for mercator, 'global' for global plot (kav7), 'nsper' for Terra/Aqua's eye view.
        
        Returns
        -------
        Plot of data subsampled to 5 km x 5 km resolution.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/9/2016, Seattle, WA
        """
        if not type(data) == str:
            print 'Error: Data must be a valid key to the "ds" dictionary.'
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
            m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0])
            m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])
        elif projection == 'global':
            m = Basemap(lon_0=0,projection='kav7',resolution='c')
            m.drawparallels(np.arange(-180,180,15),labels=[1,0,0,0])
            m.drawmeridians(np.arange(0,360,45),labels=[1,1,0,1])
        elif projection == 'satellite':
            m = Basemap(projection='nsper',lon_0=self.lon.mean(),lat_0=self.lat.mean(),\
            resolution='l',satellite_height=705000)
            m.bluemarble(alpha=.75)
        else:
            print "Error: Projection must be 'merc', 'global', or 'satellite'."
            return
        m.drawcoastlines()
        m.drawcountries()
        d = self.ds['%s' % data]
        #Might make these customized in future, for now generic
        cmap = 'viridis'
        m.pcolormesh(self.lon,self.lat,d[:np.shape(self.lon)[0],:np.shape(self.lat)[1]],\
        shading='gouraud',cmap=cmap,latlon=True)
        cbar = m.colorbar()
        cbar.ax.tick_params(labelsize=size-4) 
        plt.title('%s for %s %s, %s from %s' % \
        (self.ds_name['%s' % data],self.month,self.day,self.year,self.satellite), fontname=font,fontsize=size)
    
    def AOD_plot(self):
        """
        Plot both ACAODs and clear sky AOD.
        
        Returns
        -------
        Figure with three subplots.
        
        Modification history
        --------------------
        Written: Michael Diamond, 8/30/2016, Swakopmund, Namibia
        """
        plt.clf()
        font = 'Arial'
        size = 16
        lat = self.lat
        lon = self.lon
        max_lon = lon.max()
        max_lat = lat.max()
        min_lon = lon.min()
        min_lat = lat.min()
        plt.subplot(1,3,1)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
        x, y = m(min_lon, max_lat-1)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        plt.text(x,y,'%s %s, %s' % (self.month, self.day, self.year),\
        bbox=dict(facecolor='w', alpha=1),fontname=font,fontsize=size-4)
        key = 'Above_Cloud_AOD'
        t = 'ACAOD'
        c = 'inferno_r'
        vmin = 0
        vmed = 2.5
        vmax = 5
        d = self.ds['%s' % key]
        m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
        cmap=c,latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2., vmed, (vmed+vmax)/2., vmax],\
        location='bottom',pad=.35)
        cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2., vmed, (vmed+vmax)/2., vmax])
        cbar.ax.tick_params(labelsize=size-4)
        plt.title('A) %s' % t,fontname=font,fontsize=size)
        
        plt.subplot(1,3,2)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
        x, y = m(min_lon, max_lat-1)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        plt.text(x,y,'Time: %s' % (self.time),bbox=dict(facecolor='w', alpha=1),\
        fontname=font,fontsize=size-4)
        key = 'Above_Cloud_AOD_ModAbsAero'
        t = 'ACAOD (ModAbsAero)'
        c = 'inferno_r'
        vmin = 0
        vmed = 2.5
        vmax = 5
        d = self.ds['%s' % key]
        m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
        cmap=c,latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2., vmed, (vmed+vmax)/2., vmax],\
        location='bottom',pad=.35)
        cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2., vmed, (vmed+vmax)/2., vmax])
        cbar.ax.tick_params(labelsize=size-4)
        plt.title('B) %s' % t,fontname=font,fontsize=size)
        
        plt.subplot(1,3,3)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
        x, y = m(min_lon, max_lat-1)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        plt.text(x,y,'Satellite: %s' % (self.satellite),bbox=dict(facecolor='w', alpha=1),\
        fontname=font,fontsize=size-4)
        key = 'Clear_Sky_AOD'
        t = 'Clear Sky AOD'
        c = 'inferno_r'
        vmin = 0
        vmed = 2.5
        vmax = 5
        d = self.ds['%s' % key]
        m.pcolormesh(lon,lat,d[:np.shape(lon)[0],:np.shape(lat)[1]],\
        cmap=c,latlon=True,vmin=vmin,vmax=vmax)
        cbar = m.colorbar(ticks=[vmin, (vmin+vmed)/2., vmed, (vmed+vmax)/2., vmax],\
        location='bottom',pad=.35)
        cbar.ax.set_xticklabels([vmin, (vmin+vmed)/2., vmed, (vmed+vmax)/2., vmax])
        cbar.ax.tick_params(labelsize=size-4)
        plt.title('C) %s' % t,fontname=font,fontsize=size)
        
        plt.show()

class nrt_comp(object):
    """
    Compare Meyer et al. 2015 ACAOD with nrtMOD06 biases
    
    ***Written for use on ORACLES field campaign (https://espo.nasa.gov/home/oracles/)***
    
    Parameters
    ----------
    cfile : string
    NRT MOD06_L2 file.
    
    afile : string
    NRT MOD06ACAERO file.
    
    *Files should be from same time*
    
    Modification history
    --------------------
    Written: Michael Diamond, 08/30/2016, Swakopmund, Namibia
    Modified: Michael Diamond, 08/31/2016, P3 flying over SEA
        -Finished plots
    """
    
    def __init__(self,cfile,afile):
        #Read in files
        c = nrtMOD06(cfile)
        a = nrtACAERO(afile)
        self.lon = a.lon
        self.lat = a.lat
        self.time = c.time
        self.month = c.month
        self.day = c.day
        self.jday = c.jday
        self.year = c.year
        self.satellite = c.satellite
        Delta = u"\u0394"
        Del = u"\u03B4"
        micron = u"\u03BC"+'m'
        
        #Make dictionaries
        self.ds = {}
        self.name = {}
        self.cmap = {}
        self.v = {}
        self.units = {}
        self.keys = []
        
        #Get ACAOD data
        key = 'ACAOD'
        ACAOD = a.ds['Above_Cloud_AOD']
        ACAOD = ma.MaskedArray(ACAOD,a.ds['Above_Cloud_AOD_Uncertainty'] > 100,fill_value=0).filled()
        invalid = a.ds['Cloud_Optical_Thickness']<4
        ACAOD = ma.MaskedArray(ACAOD, invalid)
        self.ds['%s' % key] = ACAOD
        self.name['%s' % key] = 'ACAOD'
        self.cmap['%s' % key] = 'inferno_r'
        self.v['%s' % key] = (0,1.5,3)
        self.units['%s' % key] = 'unitless'
        self.keys.append(key)
        
        #Get delta ref data
        key = 'delta_ref16'
        self.ds['%s' % key] = c.ds['%s' % key]
        self.name['%s' % key] = '%sref' % Delta
        self.cmap['%s' % key] = 'RdYlBu_r'
        self.v['%s' % key] = (-6,0,6)
        self.units['%s' % key] = '%s' % micron
        self.keys.append(key)
        
        #Get del ref data
        key = 'del_ref16'
        self.ds['%s' % key] = c.ds['%s' % key]
        self.name['%s' % key] = '%sref' % Del
        self.cmap['%s' % key] = 'RdYlBu_r'
        self.v['%s' % key] = (-500,0,500)
        self.units['%s' % key] = 'per mil'
        self.keys.append(key)
        
        #Get delta COT data
        key = 'delta_COT16'
        self.ds['%s' % key] = c.ds['%s' % key]
        self.name['%s' % key] = '%sCOT' % Delta
        self.cmap['%s' % key] = 'RdYlBu_r'
        self.v['%s' % key] = (-1.,0,1.)
        self.units['%s' % key] = 'unitless'
        self.keys.append(key)
        
        #Get del COT data
        key = 'del_COT16'
        self.ds['%s' % key] = c.ds['%s' % key]
        self.name['%s' % key] = '%sCOT' % Del
        self.cmap['%s' % key] = 'RdYlBu_r'
        self.v['%s' % key] = (-100,0,100)
        self.units['%s' % key] = 'per mil'
        self.keys.append(key)
        
        #Get delta Nd data
        key = 'delta_Nd16'
        self.ds['%s' % key] = c.ds['%s' % key]
        self.name['%s' % key] = '%sNd' % Delta
        self.cmap['%s' % key] = 'RdYlBu'
        self.v['%s' % key] = (-300,0,300)
        self.units['%s' % key] = '%s' % micron
        self.keys.append(key)
        
        #Get del Nd data
        key = 'del_Nd16'
        self.ds['%s' % key] = c.ds['%s' % key]
        self.name['%s' % key] = '%sNd' % Del
        self.cmap['%s' % key] = 'RdYlBu'
        self.v['%s' % key] = (-1000,0,1000)
        self.units['%s' % key] = 'per mil'
        self.keys.append(key)        
        
    def compare(self,key='delta_ref'):
        """
        Compare ACAOD with bias given by the key.
        
        Parameters
        ----------
        key : string
        Name of bias to plot.
        
        Modification history
        --------------------
        Written: Michael Diamond, 08/30/2016, Swakopmund, Namibia
        """
        micron = u"\u03BC"+'m'
        plt.clf()
        font = 'Arial'
        size = 16
        max_lon = self.lon.max()
        max_lat = self.lat.max()
        min_lon = self.lon.min()
        min_lat = self.lat.min()
        #ACAOD
        vmin = self.v['ACAOD'][0]
        vmed = self.v['ACAOD'][1]
        vmax = self.v['ACAOD'][-1]
        cmap = self.cmap['ACAOD']
        plt.subplot(1,2,1)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.pcolormesh(self.lon,self.lat,self.ds['ACAOD'],cmap=cmap,vmin=vmin,vmax=vmax,latlon=True)
        ticks = [0,.5,1,1.5,2,2.5,3]
        cbar = m.colorbar(ticks=ticks,extend='max')
        cbar.ax.set_xticklabels(ticks)
        cbar.ax.tick_params(labelsize=size-4)
        plt.title('A) ACAOD on %s %s, %s, from %s' % (self.month,self.day,self.year,self.satellite),\
        fontname=font,fontsize=size)
        #Bias
        bias = self.ds[key]
        vmin = self.v[key][0]
        vmed = self.v[key][1]
        vmax = self.v[key][-1]
        cmap = self.cmap[key]
        plt.subplot(1,2,2)
        m = Basemap(llcrnrlon=min_lon-1,llcrnrlat=min_lat-1,urcrnrlon=max_lon+1,\
        urcrnrlat=max_lat+1,projection='merc',resolution='l')
        m.drawparallels(np.arange(-180,180,10),labels=[1,0,0,0],fontsize=size-2,fontname=font)
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1],fontsize=size-2,fontname=font)
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents('k',zorder=0)
        m.pcolormesh(self.lon,self.lat,bias,cmap=cmap,vmin=vmin,vmax=vmax,latlon=True)
        ticks = [vmin, (vmin+vmed)/2, vmed, (vmed+vmax)/2, vmax]
        cbar = m.colorbar(ticks=ticks,extend='both')
        cbar.ax.set_xticklabels(ticks)
        cbar.ax.tick_params(labelsize=size-4)
        cbar.set_label('[%s]' % self.units['%s' % key],fontname=font,fontsize=size-2)
        plt.title('B) %s (2.13 %s - 1.63 %s)' % (self.name[key],micron,micron), \
        fontname=font,fontsize=size)
        plt.show()