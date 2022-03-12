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
     
casalog.setlogfile(os.path.join(log_directory,'final_imaging.log'))

with open('cal_script_param_objs.pkl','rb') as fl:
       cal_script_param_dict=pickle.load(fl)
       initial_ms=pickle.load(fl)

for band in cal_script_param_dict:
    for msfl in cal_script_param_dict[band].split_ms_fils:
        src_nm=msfl.split('/')[-1].split('.ms')[0]
        src_im_root=os.path.join(cal_script_param_dict[band].final_ims,src_nm)
        cal_source=src_nm.split('_')[1]
        target_source=src_nm.split('_')[0]
        if not os.path.exists(src_im_root):
            os.mkdir(src_im_root)
        casalog.setlogfile(os.path.join(src_im_root,src_nm+'.log'))
        src_im_cal_dir=os.path.join(src_im_root,'cal_dir')
        src_im_cal_root=os.path.join(src_im_cal_dir,src_nm)
        if not os.path.exists(src_im_cal_dir):
            os.mkdir(src_im_cal_dir)
        src_im_img_dir=os.path.join(src_im_root,'img_dir')
        src_im_img_root=os.path.join(src_im_img_dir,src_nm)
        if not os.path.exists(src_im_img_dir):
            os.mkdir(src_im_img_dir)
        current_ms=os.path.join(msfl)
        listobs(vis=current_ms)
        imsz=512
        high_noisethreshold_1,low_noisethreshold_1,no_mask_created_1=ufpu.image_loop(ms_file=current_ms,
                                                                                im_name=src_im_img_root+'_1',
                                                                                prllel=False,im_size=imsz,
                                                                                fieldname=cal_source)
        if no_mask_created_1:
            imsz=2048
            high_noise_threshold_1,low_noise_threshold_1,no_mask_created_11=ufpu.image_loop(ms_file=current_ms,
                                                                                       im_name=src_im_img_root+'_big_im_1',
                                                                                       prllel=False,
                                                                                       im_size=imsz,
                                                                                       fieldname=cal_source)
        else:
            no_mask_created_11=False
        if no_mask_created_11:
            imsz=512
            high_noise_threshold_12,low_noise_threshold_12,mask_created_12=ufpu.image_loop(ms_file=current_ms,
                                                                                      im_name=src_im_img_root+'_user_msk_1',
                                                                                      prllel=False,
                                                                                      mask_type='user',
                                                                                      n_iter=50,
                                                                                      clean_to=1.0,fieldname=cal_source)
        pflags_for_solint={}
        #Need try/except in here somewhere in case all data gets flagged and no cal table is made. 
        for tme in range(20,180,20):
            p_solint=str(tme)+'s'
            gaincal(vis=current_ms,
                    field=cal_source,
                    caltable=src_im_cal_root+'_temp.p',
                    solint=p_solint, 
                    refant=cal_script_param_dict[band].refants, 
                    minblperant=3, 
                    gaintype='G', 
                    calmode='p')
            try:
                tb.open(src_im_cal_root+'_temp.p'.format(tme))
                total_flags=0
                total_datapoints=0
                for i in range(len(tb.getcol('FLAG'))):
                    for j in range(len(tb.getcol('FLAG')[i])):
                        total_flags+=float(sum(tb.getcol('FLAG')[i][j]))
                        total_datapoints+=float(len(tb.getcol('FLAG')[i][j]))
                pflags_for_solint[tme]=total_flags/total_datapoints
                tb.close()
            except:
                pflags_for_solint[tme]=float(200-tme)/200
            rmtables(src_im_cal_root+'_temp.p')
            with open(os.path.join(src_im_root,'solint_percentage_{}.txt'.format(src_nm)),'a') as fl:
                fl.write('{},{}\n'.format(tme,pflags_for_solint[tme]))
        tmp_solint=500
        tmp_percent=1
        for slnt in pflags_for_solint:
            if (pflags_for_solint[slnt]<=tmp_percent) & (slnt<tmp_solint):
                tmp_solint=slnt
                tmp_percent=pflags_for_solint[slnt]
        p_solint=str(tmp_solint)+'s'
        gaincal(vis=current_ms,
                field=cal_source, 
                caltable=src_im_cal_root+'.p',
                solint=p_solint, 
                refant=cal_script_param_dict[band].refants,
                minblperant=3, 
                gaintype='G', 
                calmode='p')
        applycal(vis=current_ms, 
                 field=cal_source,
                 gaintable=[src_im_cal_root+'.p'],
                 gainfield=[cal_source],
                 interp=['linear'])
        if no_mask_created_11:
            high_noisethreshold_2,low_noisethreshold_2,no_mask_created_21=ufpu.image_loop(ms_file=current_ms,
                                                                                     im_name=src_im_img_root+'_2',
                                                                                     clean_to=1.0,
                                                                                     start_high_mask_threshold=high_noisethreshold_1,
                                                                                     low_mask_threshold=low_noisethreshold_1,
                                                                                     prllel=False,
                                                                                     im_size=imsz,
                                                                                     fieldname=cal_source)
        else:
            high_noisethreshold_2,low_noisethreshold_2,no_mask_created_21=ufpu.image_loop(ms_file=current_ms,
                                                                                     im_name=src_im_img_root+'_2',
                                                                                     clean_to=5.0,
                                                                                     start_high_mask_threshold=high_noisethreshold_1,
                                                                                     low_mask_threshold=low_noisethreshold_1,
                                                                                     prllel=False,
                                                                                     im_size=imsz,
                                                                                     fieldname=cal_source)
        ap_solint=str(3*tmp_solint)+'s'
        gaincal(vis=current_ms,
                caltable=src_im_cal_root+'.ap',
                field=cal_source,
                solint=ap_solint,
                refant=cal_script_param_dict[band].refants,
                minblperant=4,
                solnorm=True,
                gaintype='G',
                calmode='ap',
                gaintable=[src_im_cal_root+'.p'],
                gainfield=[cal_source])
        applycal(vis=current_ms,
                gaintable=[src_im_cal_root+'.p',src_im_cal_root+'.ap'],
                gainfield=[cal_source,cal_source],
                interp=['linear','linear'])
        if no_mask_created_21:
            x,y,z=ufpu.image_loop(ms_file=current_ms,
                             im_name=src_im_img_root+'_3',
                             clean_to=1.0,
                             start_high_mask_threshold=5.0,
                             end_high_mask_threshold=4.0,
                             low_mask_threshold=3.0,
                             prllel=False,
                             sdlobthresh=0.7,
                             im_size=imsz,
                             fieldname=cal_source)
        else:
            x,y,z=ufpu.image_loop(ms_file=current_ms,
                             im_name=src_im_img_root+'_3',
                             clean_to=3.0,
                             start_high_mask_threshold=5.0,
                             end_high_mask_threshold=4.0,
                             low_mask_threshold=3.0,
                             prllel=False,
                             sdlobthresh=0.9,
                             im_size=imsz,
                             fieldname=cal_source)
        applycal(vis=current_ms,
                 field=target_source,
                 gaintable=[src_im_cal_root+'.p',src_im_cal_root+'.ap'],
                 gainfield=[cal_source,cal_source],
                 interp=['linear','linear'],
                 applymode='calonly',
                 parang=False) 
        if os.path.exists(current_ms.split(target_source)[0]+target_source+'.ms'):
            rmtables(current_ms.split(target_source)[0]+target_source+'.ms')
        split(vis=current_ms,
              outputvis=current_ms.split(target_source)[0]+target_source+'.ms',
              field=target_source,
              datacolumn='corrected')
        
with open('cal_script_param_objs.pkl','wb') as fl:
     pickle.dump(cal_script_param_dict,fl,protocol=pickle.HIGHEST_PROTOCOL)
     pickle.dump(initial_ms,fl,protocol=pickle.HIGHEST_PROTOCOL)

