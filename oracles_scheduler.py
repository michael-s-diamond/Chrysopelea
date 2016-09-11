"""
Retrieve MODIS and SEVIRI data and plot at specified times

***Written for ORACLES NASA ESPO mission***

Modification history
--------------------
Written: Michael Diamond, 08/10-23/2016, Seattle, WA
Modified: Michael Diamond, 09/08/2016, Swakopmund, Namibia
    -Major streamlining of MODIS data
    -Get rid of mapmaker; now entirely in Chrysopelea functions
"""

#Import libraries
import schedule
import time
import datetime
import os
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
from modipy import julian_day
import ftp_MODIS
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import ftp_SEVIRI
os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
import daily_MODIS
os.chdir('/Users/michaeldiamond/')

#Set up directories and perl scripts each night
def daily_reset():
    #Get time
    now = datetime.datetime.utcnow()
    jday = julian_day(now.month,now.day,now.year)
    if now.month < 10: month = '0'+str(now.month)
    else: month = str(now.month)
    if now.day < 10: day = '0'+str(now.day)
    else: day = str(now.day)
    
    #Make directories
    print 'Creating directories for %s/%s/%s...' % (now.month,now.day,now.year)
    file_directory = '/Users/michaeldiamond/Documents/oracles_files'
    image_directory = '/Users/michaeldiamond/Documents/oracles'
    #Create file directories
    os.chdir(file_directory+'/terra')
    os.system('mkdir ./%s' % jday)
    os.chdir(file_directory+'/aqua')
    os.system('mkdir ./%s' % jday)
    os.chdir(file_directory+'/msg')
    os.system('mkdir ./%s' % jday)
    #Create image directories
    os.chdir(image_directory+'/terra')
    os.system('mkdir ./%s' % jday)
    os.chdir(image_directory+'/aqua')
    os.system('mkdir ./%s' % jday)
    os.chdir(image_directory+'/msg')
    os.system('mkdir ./%s' % jday)
    
    #Make perl script for downloading files
    os.chdir(file_directory+'/msg/%s' % jday)
    output = open('retrieve_raw_%s.pl' % jday,'wb')
    output.write('#!/usr/bin/perl\n\n')
    output.write('use strict;\nuse warnings;\n\n')
    output.write("my $year = '%s';\nmy $month = '%s';\nmy $day = '%s';\n\n" % (now.year,month,day))
    output.write("my $username = 'oracles';\nmy $password = 'oracles2016';\n\n")
    output.write('my $host = "cloudsgate2.larc.nasa.gov";\nmy $dir = "prod/exp/oracles/d2/sat-ncdf/msg";\n\n')
    output.write('# first, get the day directory\nmy $cmd = "curl -u $username:$password http://$host/$dir/$year/$month/$day/";\n')
    output.write("my @results = `$cmd`;\n\n")
    output.write("my @download_files = ();\n\n")
    output.write("foreach (@results) {\n")
    output.write("\tif (m/\<a href\=\'(.+)(MET10.+\.C01.nc)\'\stitle/) {\n")
    output.write('\t\tpush @download_files, { full => "$1$2", fn => $2 };\n')
    output.write('\t}\n}\n\n')
    output.write("foreach (@results) {\n")
    output.write("\tif (m/\<a href\=\'(.+)(MET10.+\.C02.nc)\'\stitle/) {\n")
    output.write('\t\tpush @download_files, { full => "$1$2", fn => $2 };\n')
    output.write('\t}\n}\n\n')
    output.write('# now, go get all the files...\nmy $count = 1;\n')
    output.write("foreach my $target_file (@download_files) {\n")
    output.write('\tif (-e "./$target_file->{fn}") {\n')
    output.write('\t} else {\n')
    output.write('\t\t$cmd = "curl -u $username:$password -o $target_file->{fn} http://$host$target_file->{full} ";\n')
    output.write("\t\tmy $result = `$cmd`;\n\n")
    output.write('\t}\n\t$count++;\n}')
    output.close()
    output = open('retrieve_prod_%s.pl' % jday,'wb')
    output.write('#!/usr/bin/perl\n\n')
    output.write('use strict;\nuse warnings;\n\n')
    output.write("my $year = '%s';\nmy $month = '%s';\nmy $day = '%s';\n\n" % (now.year,month,day))
    output.write("my $username = 'oracles';\nmy $password = 'oracles2016';\n\n")
    output.write('my $host = "cloudsgate2.larc.nasa.gov";\nmy $dir = "prod/exp/oracles/d2/prod-ncdf/msg";\n\n')
    output.write('# first, get the day directory\nmy $cmd = "curl -u $username:$password http://$host/$dir/$year/$month/$day/";\n')
    output.write("my @results = `$cmd`;\n\n")
    output.write("my @download_files = ();\n\n")
    output.write("foreach (@results) {\n")
    output.write("\tif (m/\<a href\=\'(.+)(MET10.+\.nc.gz)\'\stitle/) {\n")
    output.write('\t\tpush @download_files, { full => "$1$2", fn => $2 };\n')
    output.write('\t}\n}\n\n')
    output.write('# now, go get all the files...\nmy $count = 1;\n')
    output.write("foreach my $target_file (@download_files) {\n")
    output.write('\tif (-e "./$target_file->{fn}") {\n')
    output.write('\t} else {\n')
    output.write('\t\t$cmd = "curl -u $username:$password -o $target_file->{fn} http://$host$target_file->{full} ";\n')
    output.write("\t\tmy $result = `$cmd`;\n\n")
    output.write('\t}\n\t$count++;\n}')
    output.close()
    
    os.chdir('/Users/michaeldiamond')
    print 'Done!\n'

schedule.every().day.at("17:01").do(daily_reset)

#Run daily_MODIS at 2200 UTC every night
def run_daily_MODIS():
    print 'Running daily_MODIS...'
    os.chdir('/Users/michaeldiamond/GitHub/Chrysopelea')
    try: reload(daily_MODIS)
    except: 
        print 'Something went wrong at %s' % (datetime.datetime.utcnow())
    os.chdir('/Users/michaeldiamond')
    print 'Done!\n'
    
schedule.every().day.at("13:00").do(run_daily_MODIS)

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

#Run ftp_SEVIRI every 5 min
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
