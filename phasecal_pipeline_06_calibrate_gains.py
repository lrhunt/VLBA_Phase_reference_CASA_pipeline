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
     
casalog.setlogfile(os.path.join(log_directory,'making_gaincorrection_table.log'))

with open('cal_script_param_objs.pkl','rb') as fl:
       cal_script_param_dict=pickle.load(fl)
       initial_ms=pickle.load(fl)

amp_dirs = {}
all_gaincal_vals={}
unflagged_gaincal_vals={}
print('Program Check_Antenna_Gains has started')
for band in cal_script_param_dict:
     amp_dirs[band]=[]
     all_gaincal_vals[band]={}
     unflagged_gaincal_vals[band]={}
     for path,dr,fil in os.walk(cal_script_param_dict[band].exp_im_dir):
          if (path.endswith('.ap')) & ('gain_cal' in path):
               amp_dirs[band].append(path)

     print('Found .ap calibration tables for {} band'.format(band))


     print('Making gaincal plot for experiment {expr}'.format(expr=cal_script_param_dict[band].experiment_name))               
     for aptab in amp_dirs[band]:
          tb.open(aptab)
          antenna_array=tb.getcol('ANTENNA1')
          gaincal_val_array=np.real(tb.getcol('CPARAM')[0][0])
          flag_array=tb.getcol('FLAG')[0][0]
          for i in range(0,len(antenna_array)):
               try:
                    all_gaincal_vals[band][antenna_array[i]].append(gaincal_val_array[i])
               except:
                    all_gaincal_vals[band][antenna_array[i]]=[gaincal_val_array[i]]
               if not flag_array[i]:
                    try:
                         unflagged_gaincal_vals[band][antenna_array[i]].append(gaincal_val_array[i])
                    except:
                         unflagged_gaincal_vals[band][antenna_array[i]]=[gaincal_val_array[i]]
     tb.close()
     all_gaincal_outfile = os.path.join(cal_script_param_dict[band].gaincal_csv_dir,'all_gain_vals_{}.csv'.format(band))
     unflagged_gaincal_outfile = os.path.join(cal_script_param_dict[band].gaincal_csv_dir,'unflagged_gain_vals_{}.csv'.format(band))
     csv.field_size_limit(sys.maxsize)
     with open(all_gaincal_outfile,'w') as csv_file:
          writer=csv.writer(csv_file)
          for key,value in all_gaincal_vals[band].items():
               writer.writerow([key,value])
     with open(unflagged_gaincal_outfile,'w') as csv_file:
          writer=csv.writer(csv_file)
          for key,value in unflagged_gaincal_vals[band].items():
               writer.writerow([key,value])

with open('cal_script_param_objs.pkl','wb') as fl:
     pickle.dump(cal_script_param_dict,fl,protocol=pickle.HIGHEST_PROTOCOL)
     pickle.dump(initial_ms,fl,pickle.HIGHEST_PROTOCOL)
