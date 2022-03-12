import os
import numpy as np
import time
from distutils.spawn import find_executable
import csv
import UF_PIPELINE_UTILS as ufpu
import pickle

log_directory=os.path.join(os.getcwd(),'cal_logs')
if not os.path.exists(log_directory):
     os.mkdir(log_directory)
     
casalog.setlogfile(os.path.join(log_directory,'initial_phase_cal.log'))

with open('cal_script_param_objs.pkl','rb') as fl:
       cal_script_param_dict=pickle.load(fl)
       initial_ms=pickle.load(fl)
#**********************************************************
# X-Band: Find calibrator scans and reference antenna lists
#**********************************************************

 
rfnt='LA'
highest_pflags_dict={}
highest_snr_dict={}
best_scan_dict={}
for band in cal_script_param_dict:
    if cal_script_param_dict[band].cal_file_root+'.sbd' in cal_script_param_dict[band].apply_cal_tables_list:
        cal_script_param_dict[band].apply_cal_tables_list.remove(cal_script_param_dict[band].cal_file_root+'.sbd')
        cal_script_param_dict[band].apply_cal_interp_list.remove('nearest')
        cal_script_param_dict[band].apply_cal_spwmap_list.remove([])
    highest_pflags_dict[band]={}
    highest_snr_dict[band]={}
    best_scan_dict[band]={}
    rfnt='LA'
    cal_script_param_dict[band].find_strong_cal_source(refants=rfnt,slint='inf')
    ufpu.write_cal_dict(pflag_dict=cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'],
                        snr_dict=cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'],
                        outfile=os.path.join(cal_script_param_dict[band].sn_csv_dir,'Phase_Calibration_1_{}_{}.csv'.format(band,rfnt)))
    highest_pflags_dict[band][rfnt]=1.0
    highest_snr_dict[band][rfnt]=1.0
    best_scan_dict[band][rfnt]=0
    for scan in cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags']:
        if cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'][scan]<=highest_pflags_dict[band][rfnt]:
            if cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'][scan]>highest_snr_dict[band][rfnt]:
                best_scan_dict[band][rfnt]=scan
                highest_pflags_dict[band][rfnt]=cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'][scan]
                highest_snr_dict[band][rfnt]=cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'][scan]
    calscan_indx=np.where(np.array(cal_script_param_dict[band].cal_scan_list)==best_scan_dict[band][rfnt])[0][0]
    if calscan_indx<3: 
        short_calscan_list=cal_script_param_dict[band].cal_scan_list[calscan_indx:calscan_indx+6] 
    elif calscan_indx>len(cal_script_param_dict[band].cal_scan_list)-3: 
        short_calscan_list=cal_script_param_dict[band].cal_scan_list[calscan_indx-6:calscan_indx] 
    else: 
        short_calscan_list=cal_script_param_dict[band].cal_scan_list[calscan_indx-3:calscan_indx+3] 
    short_calscan_list=cal_script_param_dict[band].cal_scan_list[calscan_indx-3:calscan_indx+3]
    if 'FD' in cal_script_param_dict[band].inner_antennas_available:
        rfnt='FD'
        cal_script_param_dict[band].find_strong_cal_source(refants=rfnt,cal_scan_list=short_calscan_list,slint='inf')
        ufpu.write_cal_dict(pflag_dict=cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'],
                       snr_dict=cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'],
                       outfile=os.path.join(cal_script_param_dict[band].sn_csv_dir,'Phase_Calibration_1_{}_{}.csv'.format(band,rfnt)))
        highest_pflags_dict[band][rfnt]=1.0
        highest_snr_dict[band][rfnt]=1.0
        best_scan_dict[band][rfnt]=0
        for scan in cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags']:
               if cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'][scan]<=highest_pflags_dict[band][rfnt]:
                    if cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'][scan]>highest_snr_dict[band][rfnt]:
                        best_scan_dict[band][rfnt]=scan
                        highest_pflags_dict[band][rfnt]=cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'][scan]
                        highest_snr_dict[band][rfnt]=cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'][scan]
    if 'PT' in cal_script_param_dict[band].inner_antennas_available:
        rfnt='PT'
        cal_script_param_dict[band].find_strong_cal_source(refants=rfnt,cal_scan_list=short_calscan_list,slint='inf')
        ufpu.write_cal_dict(pflag_dict=cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'],
                       snr_dict=cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'],
                       outfile=os.path.join(cal_script_param_dict[band].sn_csv_dir,'Phase_Calibration_1_{}_{}.csv'.format(band,rfnt)))
        highest_pflags_dict[band][rfnt]=1.0
        highest_snr_dict[band][rfnt]=1.0
        best_scan_dict[band][rfnt]=0
        for scan in cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags']:
               if cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'][scan]<=highest_pflags_dict[band][rfnt]:
                    if cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'][scan]>highest_snr_dict[band][rfnt]:
                        best_scan_dict[band][rfnt]=scan
                        highest_pflags_dict[band][rfnt]=cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'][scan]
                        highest_snr_dict[band][rfnt]=cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'][scan]
    if 'KP' in cal_script_param_dict[band].inner_antennas_available:
        rfnt='KP'
        cal_script_param_dict[band].find_strong_cal_source(refants=rfnt,cal_scan_list=short_calscan_list,slint='inf')
        ufpu.write_cal_dict(pflag_dict=cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'],
                       snr_dict=cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'],
                       outfile=os.path.join(cal_script_param_dict[band].sn_csv_dir,'Phase_Calibration_1_{}_{}.csv'.format(band,rfnt)))
        highest_pflags_dict[band][rfnt]=1.0
        highest_snr_dict[band][rfnt]=1.0
        best_scan_dict[band][rfnt]=0
        for scan in cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags']:
               if cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'][scan]<=highest_pflags_dict[band][rfnt]:
                    if cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'][scan]>highest_snr_dict[band][rfnt]:
                        best_scan_dict[band][rfnt]=scan
                        highest_pflags_dict[band][rfnt]=cal_script_param_dict[band].refant_selection_dict[rfnt]['pflags'][scan]
                        highest_snr_dict[band][rfnt]=cal_script_param_dict[band].refant_selection_dict[rfnt]['snr'][scan]               
    for ant in best_scan_dict[band]:
        casalog.post('Best Scan {} = {}\n'.format(ant,best_scan_dict[band][ant]))
    lowest_pflags_overall=1.0
    highest_snr_overall=1.0
    for rfnt in highest_snr_dict[band]:
        if highest_pflags_dict[band][rfnt]<=lowest_pflags_overall:
            if highest_snr_dict[band][rfnt]>highest_snr_overall:
                refant=rfnt
                lowest_pflags_overall=highest_pflags_dict[band][rfnt]
                highest_snr_overall=highest_snr_dict[band][rfnt]
                calscan=best_scan_dict[band][rfnt]
    cal_script_param_dict[band].phasecal_scan=calscan


    refants=ufpu.make_refant_list(inner_ants=np.array(cal_script_param_dict[band].inner_antennas_available),
                             outer_ants=np.array(cal_script_param_dict[band].outer_antennas_available),
                             initial_ant=refant)

    cal_script_param_dict[band].save_refants(refants)
    fringefit(vis=cal_script_param_dict[band].ms_file,
              caltable=cal_script_param_dict[band].cal_file_root+'.sbd',
              scan=str(cal_script_param_dict[band].phasecal_scan),
              gaintable=cal_script_param_dict[band].apply_cal_tables_list,
              interp=cal_script_param_dict[band].apply_cal_interp_list,
              solint='inf',
              zerorates=True,
              refant=cal_script_param_dict[band].refants,
              minsnr=5.0,
              parang=True,
              globalsolve=False)
    cal_script_param_dict[band].apply_cal_tables_list.append(cal_script_param_dict[band].cal_file_root+'.sbd')
    cal_script_param_dict[band].apply_cal_interp_list.append('nearest')
    cal_script_param_dict[band].apply_cal_spwmap_list.append([])
with open('cal_script_param_objs.pkl','wb') as fl:
     pickle.dump(cal_script_param_dict,fl,protocol=pickle.HIGHEST_PROTOCOL)
     pickle.dump(initial_ms,fl,pickle.HIGHEST_PROTOCOL)

          
# Write diagnostic files (SBD cal scan (s/x), refants (s/x), cal scan list)

#with open(os.path.join(supplement_cal_dir,'cal_scan_list.txt'),'w') as csl:
#     for scn in cal_scan_list:
#          if scn==cal_scan_list[-1]:
#               csl.write('{}'.format(scn))
#          else:
#               csl.write('{},'.format(scn))
               
#with open(os.path.join(supplement_cal_dir,'refant_list.txt'),'w') as rfnts:
#     rfnts.write('x-band={}'.format(refants_x))
#     rfnts.write('s-band={}'.format(refants_s))
     
#with open(os.path.join(supplement_cal_dir,'sbd_calibration_scan.txt'),'w') as sbdcal:
#     sbdcal.write('x-band={}'.format(calscan_x))
#     sbdcal.write('s-band={}'.format(calscan_s))

