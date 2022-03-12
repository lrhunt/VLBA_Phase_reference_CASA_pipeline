import os
import numpy as np
import pyfits as fits
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
import UF_PIPELINE_UTILS as ufpu
import pickle

log_directory=os.path.join(os.getcwd(),'cal_logs')
if not os.path.exists(log_directory):
    os.mkdir(log_directory)

casalog.setlogfile(os.path.join(log_directory,'ty_gc_accor.log'))

with open('cal_script_param_objs.pkl','rb') as fl:
    cal_script_param_dict=pickle.load(fl)
    initial_ms=pickle.load(fl)

for band in cal_script_param_dict:
    gencal(vis=cal_script_param_dict[band].ms_file,
           caltable=cal_script_param_dict[band].cal_file_root+'.tsys',
           caltype='tsys',
           uniform=False)
    if cal_script_param_dict[band].cal_file_root+'.tsys' not in cal_script_param_dict[band].apply_cal_tables_list:
        cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.tsys')
        cal_script_param_dict[band].apply_cal_interp_list.append('nearest,nearest')
        cal_script_param_dict[band].apply_cal_spwmap_list.append([])
    if 'casa' not in locals():
        if (casatools.version()[0]>6) & (casatools.version()[1]>3):
            gencal(vis=cal_script_param_dict[band].ms_file,
                   caltable=cal_script_param_dict[band].cal_file_root+'.gc',
                   caltype='gc')
        else:
            gencal(vis=cal_script_param_dict[band].ms_file,
                   caltable=cal_script_param_dict[band].cal_file_root+'.gc',
                   caltype='gc',
                   infile=initial_ms.cal_file_root+'.gc.tbl')
    else:
        gencal(vis=cal_script_param_dict[band].ms_file,
               caltable=cal_script_param_dict[band].cal_file_root+'.gc',
               caltype='gc',
               infile=initial_ms.cal_file_root+'.gc.tbl')
    if cal_script_param_dict[band].cal_file_root+'.gc' not in cal_script_param_dict[band].apply_cal_tables_list:
        cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.gc')
        cal_script_param_dict[band].apply_cal_interp_list.append('nearest')
        cal_script_param_dict[band].apply_cal_spwmap_list.append([])
    accor(vis=cal_script_param_dict[band].ms_file,
          caltable=cal_script_param_dict[band].cal_file_root+'.accor',
          solint='20s')
    if cal_script_param_dict[band].cal_file_root+'.accor' not in cal_script_param_dict[band].apply_cal_tables_list:
        cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.accor')
        cal_script_param_dict[band].apply_cal_interp_list.append('nearest')
        cal_script_param_dict[band].apply_cal_spwmap_list.append([])
    smoothcal(vis=cal_script_param_dict[band].ms_file,
              tablein=cal_script_param_dict[band].cal_file_root+'.accor',
              smoothtime=1800.0)
    if find_executable('aoflagger'):
	    print('running aoflagger!')
	    flag_command='aoflagger {file_name}'.format(file_name=cal_script_param_dict[band].ms_file)
	    os.system(flag_command)
    else:
        flagdata(vis=cal_script_param_dict[band].ms_file,mode='tfcrop')
    flagdata(vis=cal_script_param_dict[band].ms_file,
             autocorr=True,
             mode='manual')


with open('cal_script_param_objs.pkl','wb') as fl:
     pickle.dump(cal_script_param_dict,fl,protocol=pickle.HIGHEST_PROTOCOL)
     pickle.dump(initial_ms,fl,pickle.HIGHEST_PROTOCOL)
