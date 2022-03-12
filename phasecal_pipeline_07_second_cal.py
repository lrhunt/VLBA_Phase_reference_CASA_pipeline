import os
import pyfits as fits
import numpy as np
import time
from distutils.spawn import find_executable
try:  
    from casatasks.private import tec_maps
except:
    from recipes import tec_maps
import random
import csv
import shutil
import itertools
import copy
import UF_PIPELINE_UTILS as ufpu
import pickle

log_directory=os.path.join(os.getcwd(),'cal_logs')
if not os.path.exists(log_directory):
     os.mkdir(log_directory)
     
casalog.setlogfile(os.path.join(log_directory,'second_phase_bandpass_cal.log'))

with open('cal_script_param_objs.pkl','rb') as fl:
       cal_script_param_dict=pickle.load(fl)
       initial_ms=pickle.load(fl)

for band in cal_script_param_dict:
    exprmnt=cal_script_param_dict[band].cal_file_root.split('/')[-1]
    if not os.path.exists(os.path.join(cal_script_param_dict[band].second_cal_dir,exprmnt+'.tsys')):
        shutil.copytree(cal_script_param_dict[band].cal_file_root+'.tsys',os.path.join(cal_script_param_dict[band].second_cal_dir,exprmnt+'.tsys'))
        shutil.copytree(cal_script_param_dict[band].cal_file_root+'.gc',os.path.join(cal_script_param_dict[band].second_cal_dir,exprmnt+'.gc'))
        shutil.copytree(cal_script_param_dict[band].cal_file_root+'.accor',os.path.join(cal_script_param_dict[band].second_cal_dir,exprmnt+'.accor'))
        if ufpu.test_tec_maps():
            shutil.copytree(cal_script_param_dict[band].cal_file_root+'.tec',os.path.join(cal_script_param_dict[band].second_cal_dir,exprmnt+'.tec'))
    rmtables(tablenames=cal_script_param_dict[band].ms_file)
    rmtables(tablenames=initial_ms.ms_file)
    if os.path.exists(cal_script_param_dict[band].ms_file+'.flagversions'):
        shutil.rmtree(cal_script_param_dict[band].ms_file+'.flagversions')
    if os.path.exists(initial_ms.ms_file+'.flagversions'):
        shutil.rmtree(initial_ms.ms_file+'.flagversions')     
    cal_script_param_dict[band].change_cal_file_root()
    cal_script_param_dict[band].apply_cal_tables_list=[]
    cal_script_param_dict[band].apply_cal_interp_list=[]
    cal_script_param_dict[band].apply_cal_spwmap_list=[]
    cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.tsys')
    cal_script_param_dict[band].apply_cal_interp_list.append('nearest')
    cal_script_param_dict[band].apply_cal_spwmap_list.append([])
    cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.gc')
    cal_script_param_dict[band].apply_cal_interp_list.append('nearest')
    cal_script_param_dict[band].apply_cal_spwmap_list.append([])
    cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.accor')
    cal_script_param_dict[band].apply_cal_interp_list.append('nearest')
    cal_script_param_dict[band].apply_cal_spwmap_list.append([])
    if ufpu.test_tec_maps():
        cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.tec')
        cal_script_param_dict[band].apply_cal_interp_list.append('nearest')
        cal_script_param_dict[band].apply_cal_spwmap_list.append([])


if ('casa' in locals()):
    hdu=fits.open(initial_ms.idifitsfile)
    ant_dict=ufpu.make_ant_dict(hdu)
    time_range = ''
    if not time_range:
        t = time.strptime("2000y01d00h00m00s", "%Yy%jd%Hh%Mm%Ss")
        btime = time.mktime(t) + 40587.0 * 86400
        t = time.strptime("2100y01d00h00m00s", "%Yy%jd%Hh%Mm%Ss")
        etime = time.mktime(t) + 40587.0 * 86400
    row_dict={}
    with open(initial_ms.cal_file_root+'.txt','w') as f:
        for i in range(0,len(hdu['gain_curve'].data)):
            for j in range(0,hdu['gain_curve'].header['no_band']):
                f.write(str(ufpu.find_band_letter(hdu['FREQUENCY'].header['REF_FREQ']+hdu['FREQUENCY'].data['BANDFREQ'][0][j])))
                f.write(' ')
                f.write(str(hdu['FREQUENCY'].header['REF_FREQ']+hdu['FREQUENCY'].data['BANDFREQ'][0][j]))
                f.write(' ')
                f.write(str(hdu['FREQUENCY'].header['REF_FREQ']+hdu['FREQUENCY'].data['BANDFREQ'][0][j]+hdu['frequency'].data['total_bandwidth'][0][0]))
                f.write(' ')
                f.write(str(btime))
                f.write(' ')
                f.write(str(etime))
                f.write(' ')
                f.write(str(ant_dict[hdu['gain_curve'].data['antenna_no'][i]]))
                f.write(' ')
                poly=ufpu.transform_poly(np.split(hdu['gain_curve'].data['gain_1'][i],hdu['gain_curve'].header['no_band'])[j])
                f.write(' ')
                f.write(str(poly[0]*np.sqrt(hdu['gain_curve'].data[i]['sens_1'][j])))
                f.write(' ')
                f.write(str(poly[1]*np.sqrt(hdu['gain_curve'].data[i]['sens_1'][j])))
                f.write(' ')
                f.write(str(poly[2]*np.sqrt(hdu['gain_curve'].data[i]['sens_1'][j])))
                f.write(' ')
                f.write(str(poly[3]*np.sqrt(hdu['gain_curve'].data[i]['sens_1'][j])))
                f.write(' ')
                f.write(str(poly[0]*np.sqrt(hdu['gain_curve'].data[i]['sens_1'][j])))
                f.write(' ')
                f.write(str(poly[1]*np.sqrt(hdu['gain_curve'].data[i]['sens_1'][j])))
                f.write(' ')
                f.write(str(poly[2]*np.sqrt(hdu['gain_curve'].data[i]['sens_1'][j])))
                f.write(' ')
                f.write(str(poly[3]*np.sqrt(hdu['gain_curve'].data[i]['sens_1'][j])))
                f.write('\n')

    hdu.close()
    # Defining column names for the CASA gaincurve table

    columnnames = ["BANDNAME","BFREQ","EFREQ","BTIME","ETIME","ANTENNA","GAIN"]

    # Defining datatypes for columns so that CASA can read the txt file

    datatypes = ["A","D","D","D","D","A","R4,2"]

    # creating CASA table to use for gaincurve to do initial amplitude calibration

    tb.fromascii(initial_ms.cal_file_root+'.gc.tbl', asciifile=initial_ms.cal_file_root+'.txt', sep=' ',columnnames=columnnames, datatypes=datatypes)

    #import idifits file as ms, put ms in ms_dir

importfitsidi(fitsidifile=initial_ms.idifitsfile,
              vis=initial_ms.ms_file,
              scanreindexgap_s=15.0)

#changing antenna table to get the diameters of the vlba antennas.
#This might be something that is possible in the future just by reading the idifits file?


diams= initial_ms.num_antennas*[25.0]
tb.open(initial_ms.ms_file+'/ANTENNA', nomodify=False)
tb.putcol('DISH_DIAMETER', diams)
tb.close()

#create ionosphere correction
#USNO Network is down, cannot access ionosphere calibration files. This is being commented out
#07/16/2020

#tec_maps.create(vis=initial_ms.ms_file,imname='iono')
#gencal(vis=initial_ms.ms_file,caltable=cal_file_root+'.tec',caltype='tecim',infile='iono.IGS_TEC.im')

#Flagging the first and last second of each scan, for some reason fringefit randomly fails if 
#you don't do this

flagdata(vis=initial_ms.ms_file,mode='quack',quackinterval=2.0,quackmode='endb')
flagdata(vis=initial_ms.ms_file,mode='quack',quackinterval=2.0,quackmode='beg')


tb.open(initial_ms.ms_file, nomodify=False)
scan_table=tb.getcol('SCAN_NUMBER')
new_arrayid=np.zeros(len(scan_table))

for scan in initial_ms.scn_low_bsln:
     new_arrayid[np.where(scan_table==scan)[0]]=1

tb.putcol('ARRAY_ID', new_arrayid)
tb.close()




# 5. Flag first and last second so fringefit works (quack)

# Making the required strings/files for stat gaincal
for band in cal_script_param_dict:
    spw_str='{}~{}'.format(min(cal_script_param_dict[band].spw_freq_dict[band]),max(cal_script_param_dict[band].spw_freq_dict[band]))
    split(vis=initial_ms.ms_file,outputvis=cal_script_param_dict[band].ms_file,datacolumn='data',spw=spw_str)
    for fl in os.listdir(cal_script_param_dict[band].gaincal_csv_dir):
        if 'unflagged' in fl:
            gaincal_csv_file=os.path.join(cal_script_param_dict[band].gaincal_csv_dir,fl)
            with open(gaincal_csv_file,'r') as csv_file:                                   
                reader=csv.reader(csv_file)        
                stat_gaincal_dict={rows[0]:rows[1] for rows in reader}
    stat_gaincal_ant_string=''
    stat_gaincal_prmtr_list=[]
    for ant_strng in stat_gaincal_dict:
          stat_gaincal_dict[ant_strng]=[float(f) for f in stat_gaincal_dict[ant_strng].split('[')[-1].split(']')[0].split(',')]
          stat_gaincal_ant_string=stat_gaincal_ant_string+ant_strng+','
          stat_gaincal_prmtr_list.append(np.median(stat_gaincal_dict[ant_strng]))
    stat_gaincal_ant_string = stat_gaincal_ant_string[:-1]
    gencal(vis=cal_script_param_dict[band].ms_file,
           caltable=cal_script_param_dict[band].cal_file_root+'.amp',
           caltype='amp',
           antenna=stat_gaincal_ant_string,
           parameter=stat_gaincal_prmtr_list)
    cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.amp')
    cal_script_param_dict[band].apply_cal_interp_list.append('nearest')
    cal_script_param_dict[band].apply_cal_spwmap_list.append([])
    if find_executable('aoflagger'):
     	print('running aoflagger!')
     	flag_command='aoflagger {file_name}'.format(file_name=cal_script_param_dict[band].ms_file)
     	print(flag_command)
     	os.system(flag_command)
    else:
        flagdata(vis=cal_script_param_dict[band].ms_file,mode='tfcrop')
    flagdata(vis=initial_ms.ms_file,
             autocorr=True,
             mode='manual')
    fringefit(vis=cal_script_param_dict[band].ms_file,
              caltable=cal_script_param_dict[band].cal_file_root+'.sbd',
              scan=str(cal_script_param_dict[band].phasecal_scan),
              solint='inf',
              zerorates=True,
              refant=cal_script_param_dict[band].refants,
              minsnr=5.0,
              gaintable=cal_script_param_dict[band].apply_cal_tables_list,
              interp=cal_script_param_dict[band].apply_cal_interp_list,
              parang=True,
              globalsolve=False)
    cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.sbd')
    cal_script_param_dict[band].apply_cal_interp_list.append('nearest')
    cal_script_param_dict[band].apply_cal_spwmap_list.append([])
    cal_script_param_dict[band].find_bandpass_scan(refants=cal_script_param_dict[band].refants,
                                                     gaincal_corr=cal_script_param_dict[band].apply_cal_tables_list)
    ufpu.write_cal_dict(pflag_dict={},
                        snr_dict=cal_script_param_dict[band].bandpass_snr_dict,
                        outfile=os.path.join(cal_script_param_dict[band].sn_csv_dir,'Bandpass_SNR_Dictionary_2_{}.csv'.format(band)))           
    top_bp_snr_scan=list(cal_script_param_dict[band].bandpass_snr_dict)[0]
    for scn in cal_script_param_dict[band].bandpass_snr_dict:
        if cal_script_param_dict[band].bandpass_snr_dict[scn]>cal_script_param_dict[band].bandpass_snr_dict[top_bp_snr_scan]:
            top_bp_snr_scan=scn
    cal_script_param_dict[band].bandpass_scan=top_bp_snr_scan
    bandpass(vis=cal_script_param_dict[band].ms_file,
             caltable=cal_script_param_dict[band].cal_file_root+'.bpass',
             scan=str(cal_script_param_dict[band].bandpass_scan),
             gaintable=cal_script_param_dict[band].apply_cal_tables_list,
             interp=cal_script_param_dict[band].apply_cal_interp_list,
             solnorm=True,
             solint='inf',
             refant=cal_script_param_dict[band].refants,
             bandtype='B',
             parang=True)
    cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.bpass')
    cal_script_param_dict[band].apply_cal_interp_list.append('linear,linear')
    cal_script_param_dict[band].apply_cal_spwmap_list.append([])    
    flagdata(vis=cal_script_param_dict[band].ms_file,
              mode='manual',
              spw=cal_script_param_dict[band].flag_spw_string)
    with open(os.path.join(cal_script_param_dict[band].supplement_cal_dir,'bandpass_calibration_scan_2_{}.txt'.format(band)),'w') as bpsscn:
        bpsscn.write('{}_band={}'.format(band,top_bp_snr_scan))
    fringefit(vis=cal_script_param_dict[band].ms_file,
              scan=cal_script_param_dict[band].initial_fringefit_scans_string,
              caltable=cal_script_param_dict[band].cal_file_root+'.mbd',
              solint='inf',
              refant=cal_script_param_dict[band].refants,
              minsnr=5.0,
              combine='spw',
              gaintable=cal_script_param_dict[band].apply_cal_tables_list,
              interp=cal_script_param_dict[band].apply_cal_interp_list,
              parang=True,
              globalsolve=False)
    cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.mbd')
    cal_script_param_dict[band].apply_cal_interp_list.append('linear')
    cal_script_param_dict[band].apply_cal_spwmap_list.append(cal_script_param_dict[band].nspw*[0])
    applycal(vis=cal_script_param_dict[band].ms_file,
             gaintable=cal_script_param_dict[band].apply_cal_tables_list,
             interp=cal_script_param_dict[band].apply_cal_interp_list,
             spwmap=cal_script_param_dict[band].apply_cal_spwmap_list,
             parang=True)
    cal_script_param_dict[band].phasecal_split_ms_files(band=band)
     
with open('cal_script_param_objs.pkl','wb') as fl:
     pickle.dump(cal_script_param_dict,fl,protocol=pickle.HIGHEST_PROTOCOL)
     pickle.dump(initial_ms,fl,protocol=pickle.HIGHEST_PROTOCOL)

