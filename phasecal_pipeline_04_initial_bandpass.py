import os
import numpy as np
import time
from distutils.spawn import find_executable
import random
import csv
import shutil
import itertools
import UF_PIPELINE_UTILS as ufpu
import pickle

log_directory=os.path.join(os.getcwd(),'cal_logs')
if not os.path.exists(log_directory):
     os.mkdir(log_directory)
     
casalog.setlogfile(os.path.join(log_directory,'initial_bandpass_cal.log'))

with open('cal_script_param_objs.pkl','rb') as fl:
       cal_script_param_dict=pickle.load(fl)
       initial_ms=pickle.load(fl)
#   e. Find best bandpass calibrator scan
for band in cal_script_param_dict:
    if cal_script_param_dict[band].cal_file_root+'.bpass' in cal_script_param_dict[band].apply_cal_tables_list:
        cal_script_param_dict[band].apply_cal_tables_list.remove(cal_script_param_dict[band].cal_file_root+'.bpass')
        cal_script_param_dict[band].apply_cal_interp_list.remove('linear,linear')
        cal_script_param_dict[band].apply_cal_spwmap_list.remove([])
    if cal_script_param_dict[band].cal_file_root+'.mbd' in cal_script_param_dict[band].apply_cal_tables_list:
        cal_script_param_dict[band].apply_cal_tables_list.remove(cal_script_param_dict[band].cal_file_root+'.mbd')
        cal_script_param_dict[band].apply_cal_interp_list.remove('linear')
        cal_script_param_dict[band].apply_cal_spwmap_list.remove(cal_script_param_dict[band].nspw*[0])
    cal_script_param_dict[band].find_bandpass_scan(refants=cal_script_param_dict[band].refants)
    ufpu.write_cal_dict(pflag_dict={},
                    snr_dict=cal_script_param_dict[band].bandpass_snr_dict,
                    outfile=os.path.join(cal_script_param_dict[band].sn_csv_dir,'Bandpass_SNR_Dictionary_1_{}.csv'.format(band)))       
    top_bp_snr_scan=list(cal_script_param_dict[band].bandpass_snr_dict)[0]
    for scn in cal_script_param_dict[band].bandpass_snr_dict:
        if cal_script_param_dict[band].bandpass_snr_dict[scn]>cal_script_param_dict[band].bandpass_snr_dict[top_bp_snr_scan]:
            top_bp_snr_scan=scn
#   f. Bandpass calibration
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

# 10. Phase calibaration on calibrator sources
    fringefit(vis=cal_script_param_dict[band].ms_file,
              caltable=cal_script_param_dict[band].cal_file_root+'.mbd',
              field=cal_script_param_dict[band].cal_sources_string,
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
             field=cal_script_param_dict[band].cal_sources_string,
             gaintable=cal_script_param_dict[band].apply_cal_tables_list,
             interp=cal_script_param_dict[band].apply_cal_interp_list,
             spwmap=cal_script_param_dict[band].apply_cal_spwmap_list,
             gainfield=cal_script_param_dict[band].cal_sources_string,
             parang=True)

with open('cal_script_param_objs.pkl','wb') as fl:
     pickle.dump(cal_script_param_dict,fl,protocol=pickle.HIGHEST_PROTOCOL)
     pickle.dump(initial_ms,fl,pickle.HIGHEST_PROTOCOL)

