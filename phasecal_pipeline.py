import sys
import os

for arg in sys.argv:
    if '/' in arg:
        pipe_path=arg

if pipe_path.startswith('/'):
    path_list=pipe_path[1:].split('/')
else:
    path_list=pipe_path.split('/')

path_string='/'
for strng in path_list[:-1]:
    path_string=os.path.join(path_string,strng)


execfile(os.path.join(path_string,'phasecal_pipeline_01_import.py'))
execfile(os.path.join(path_string,'phasecal_pipeline_02_ty_gc_accor_cal.py'))
execfile(os.path.join(path_string,'phasecal_pipeline_03_initial_phasecal.py'))
execfile(os.path.join(path_string,'phasecal_pipeline_04_initial_bandpass.py'))
execfile(os.path.join(path_string,'phasecal_pipeline_05_selfcal_calibrators.py'))
execfile(os.path.join(path_string,'phasecal_pipeline_06_calibrate_gains.py'))
execfile(os.path.join(path_string,'phasecal_pipeline_07_second_cal.py'))
execfile(os.path.join(path_string,'phasecal_pipeline_08_final_imaging_and_split_sources.py'))
