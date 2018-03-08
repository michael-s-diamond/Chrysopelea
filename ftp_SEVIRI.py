"""
Download script for SEVIRI NRT data

*Created for use with ORACLES NASA ESPO mission*

Modification history
--------------------
Written: Michael Diamond, 08/15/2016, Seattle, WA
Modified: Michael Diamond, 09/08/2016, Swakopmund, Namibia
    -Cloud heights and thicknesses in feet
"""

import os
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import modipy as mod
import sevipy as sev
reload(sev)
os.chdir('/Users/michaeldiamond/Documents/')
import datetime
import matplotlib.pylab as plt

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

"""
Get LARC SEVIRI data
"""
file_directory = '/Users/michaeldiamond/Documents/oracles_files/msg/%s' % jday
os.chdir(file_directory)
current_files = os.listdir(file_directory)
try: os.system('perl retrieve_prod_%s.pl' % jday)
except: print 'Error getting products at %s' % now
new_files = os.listdir(file_directory)
#Rid list of repeat files
for f in current_files:
    if f in new_files: new_files.remove(f)
for f in new_files:
    if len(f) > 35 and f[0:35] in current_files: new_files.remove(f)
    if len(f) > 38 and f[0:38] in current_files: new_files.remove(f)
#Add files if images were not made properly
directory = '/Users/michaeldiamond/Documents/oracles/msg/%s' % jday
os.chdir(directory)
current_images = os.listdir(directory)
os.chdir(file_directory)

#Add to new_files if things got overlooked before and no image was made
for f in current_files:
    if f[0] == 'M' and 5 <= int(f[14:16]) <= 19:
        time = f[14:18]
        if f[-4] == '1':
            ref_im = '%s_%s_%s_%s_CRS.png' % (year,month,day,time)
            if ref_im not in current_images: 
                new_files.append(f)
        elif f[19] == 'c':
            ref_im = '%s_%s_%s_%s_DZ.png' % (year,month,day,time)
            if ref_im not in current_images: new_files.append(f)
        elif f[19] == 'a':
            ref_im = '%s_%s_%s_%s_AOD.png' % (year,month,day,time)
            if ref_im not in current_images: new_files.append(f)
        else: pass
    else: pass

#Make plots
for f in new_files:
    if f[0] == 'M' and 5 <= int(f[14:16]) <= 19:
        if f[-4] == '1':
            if f[0:26]+'2.nc' in new_files or f[0:26]+'2.nc' in current_files:
                #Read in file
                os.chdir(file_directory)
                cr = sev.CR(f,f[0:26]+'2.nc')
                #Move to image directory
                os.chdir(directory)
                #Color ratio for each file pair
                print 'Making CR plots for %s...' % f
                plt.figure(1)
                cr.merc()
                fig = plt.gcf()
                fig.set_size_inches(13.33,7.5)
                plt.savefig('%s_%s_%s_%s_CRS' % (year,month,day,cr.time),dpi=150)
                print 'Done!\n'
            else:
                os.chdir(file_directory)
                os.system('rm %s' % f)
        elif f[19] == 'c':
            #Read in file
            os.chdir(file_directory)
            try:
                os.system('gunzip %s' % f)
                fc = f[0:38]
                cloud = sev.cloud(fc)
                #Move to image directory
                os.chdir(directory)
                #Effective radius for each file pair
                print 'Making plots for %s...' % fc
                for var in ['Re','Nd','Tau','Pbot','Ptop','Ztf','Zbf','DZ']:
                    print '...%s...' % var
                    plt.figure(1)
                    cloud.plot(var)
                    fig = plt.gcf()
                    fig.set_size_inches(13.33,7.5)
                    plt.savefig('%s_%s_%s_%s_%s' % (year,month,day,cloud.time,var),dpi=150)
                print 'Done!\n'
            except: os.system('rm %s' % f)
        elif f[19] == 'a':
            ##Read in file
            os.chdir(file_directory)
            try:
                os.system('gunzip %s' % f)
                fa = f[0:35]
                aero = sev.aero(fa)
                #Move to image directory
                os.chdir(directory)
                #Effective radius for each file pair
                print 'Making plots for %s...' % fa
                for var in ['AOD','ATYP']:
                    print '...%s...' % var
                    plt.figure(1)
                    aero.plot(var)
                    fig = plt.gcf()
                    fig.set_size_inches(13.33,7.5)
                    plt.savefig('%s_%s_%s_%s_%s' % (year,month,day,aero.time,var),dpi=150)
                print 'Done!\n'
            except: os.system('rm %s' % f)
        else: pass
    else: pass

plt.close("all")
