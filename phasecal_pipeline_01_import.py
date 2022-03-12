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

casalog.setlogfile(os.path.join(log_directory,'import.log'))

initial_ms=ufpu.cal_script_params()

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

initial_ms.get_metadata_twntfr_hr(phaseref=True)

try:
    tec_maps.create(vis=initial_ms.ms_file,doplot=True,imname='iono')
    gencal(vis=initial_ms.ms_file,caltable=initial_ms.cal_file_root+'.tec',caltype='tecim',infile='iono.IGS_TEC.im')
    initial_ms.apply_cal_tables_list.append(initial_ms.cal_file_root+'.tec')
    initial_ms.apply_cal_interp_list.append('nearest')
    initial_ms.apply_cal_spwmap_list.append([])
except:
    print('Tec Maps did not work')

# Run listobs so I can find where sources failed

listobsfile=os.path.join(initial_ms.supplement_cal_dir,initial_ms.experiment_name+'.listobs')
listobs(vis=initial_ms.ms_file)
listobs(vis=initial_ms.ms_file,listfile=listobsfile)

flagdata(vis=initial_ms.ms_file,mode='quack',quackinterval=1.0,quackmode='endb')
flagdata(vis=initial_ms.ms_file,mode='quack',quackinterval=1.0,quackmode='beg')

tb.open(initial_ms.ms_file, nomodify=False)
scan_table=tb.getcol('SCAN_NUMBER')
new_arrayid=np.zeros(len(scan_table))

for scan in initial_ms.scn_low_bsln:
     new_arrayid[np.where(scan_table==scan)[0]]=1

tb.putcol('ARRAY_ID', new_arrayid)
tb.close()

#run aoflagger on the data!
cal_script_param_dict={}
if len(initial_ms.spw_freq_dict)==1:
     cal_script_param_dict[list(initial_ms.spw_freq_dict)[0]]=copy.deepcopy(initial_ms)
else:
     for sp in initial_ms.spw_freq_dict:
          cal_script_param_dict[sp]=cal_script_params(sp)
          spw_str='{}~{}'.format(min(initial_ms.spw_freq_dict[sp]),max(initial_ms.spw_freq_dict[sp]))
          split(vis=initial_ms.ms_file,outputvis=cal_script_param_dict[sp].ms_file,datacolumn='data',spw=spw_str)
          cal_script_param_dict[sp].get_metadata_twntfr_hr(phaseref=True)
          if ufpu.test_tec_maps():
               shutil.copytree(initial_ms.cal_file_root+'.tec',os.path.join(cal_script_param_dict[sp].cal_file_root+'.tec'))

with open('cal_script_param_objs.pkl','wb') as fl:
     pickle.dump(cal_script_param_dict,fl,protocol=pickle.HIGHEST_PROTOCOL)
     pickle.dump(initial_ms,fl,protocol=pickle.HIGHEST_PROTOCOL)
