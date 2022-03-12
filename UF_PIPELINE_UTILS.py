#*********************************************************************************************************************************************
# This script defines many of the functions used throughout the rest of the CASA UF imaging script. It also imports some of the packages that
# we use throughout the rest of the script. 
#*********************************************************************************************************************************************
# Date          Ver.     Editor          Changes
# 08/23/21      1.0      Lucas Hunt      This is the initial implementation. This version includes many functions. Will list and outline below
#
#                                        find_start_end_times(listobs_file)
#                                          - This function reads through a listobs output file and creates 3 dictionaries that contain the 
#                                            start and end times for each scan. Returns start_date,start_time_dict,end_time_dict
#                                        find_bandpass_scan(ms_file,cal_file_root,cal_time_dict,refants='LA,PT,FD,KP')
#                                          - This function will take 
#                                        calibrator_dict(ms_file,cal_scan_list,num_scans=50,refants='LA,PT,FD,KP')
#                                          - uses fringefit on the scans in scan list to determine the signal to noise on each baseline for 
#                                            each spectral window. Since we have information on signal to noise, we can use this to find a scan
#                                            for single band delay and for bandpass calibration. 
#                                        write_cal_dict(pflag_dict={},snr_dict={},outfile='Printed_Dict.csv')
#                                          - This will write the above dictionary to a file to look at later if you'd like
#                                        make_refant_list(inner_ants,outer_ants,initial_ant)
#                                          - This will take an initial reference antenna (usually determined as having highest S/N using the 
#                                            two functions above this), and randomly assign an order for reference antennas based, first using
#                                            the inner antennas, and then the outer antennas.
#                                        find_bandpass_scan(ms_fl,cl_fle_rt,cl_scns,refants='LA,PT,FD,KP')
#                                          - This takes a measurement set, a root string (so that it can put the bp table in a specific spot)
#                                            a list of scans, and a reference antenna list. It then carries out the bandpass calibration on the
#                                            list of scans to try to find the scan with the highest signal to noise ratio
#                                        get_unimaged(img_dir='',splt='')
#                                          - This will look in the image directory and the split measurement set directory to determine where the
#                                            imaging script left off. 
#                                        calc_pixel_size(msin)
#                                          - This function will take in the measurement set, and an array of baseline lengths (typed out by hand
#                                            in this file) to determine the pixel size
#                                        image_loop(ms_file,
#                                                   im_name,
#                                                   im_size=512,
#                                                   robust_weight=2,
#                                                   clean_to=8.0,
#                                                   start_high_mask_threshold=15.0,
#                                                   end_high_mask_threshold=5.0,
#                                                   low_mask_threshold=10.0,
#                                                   prllel=False,
#                                                   n_iter=4000,
#                                                   mask_type='auto-multithresh',
#                                                   sdlobthresh=0.7)
#                                          - This task has many arguments, but usually the defaults are fine. It will take a measurement set
#                                            and try to make an image that matches the parameters in the arguments. If it can't do that it will
#                                            decrement the noisethreshold parameter of tclean's auto-multithresh until it is able to make an
#                                            image. It will then return the final value for noisethreshold that made an image, the value for 
#                                            lownoisethreshold (looking for faint emission) and whether a mask was made or not.
#                                        are_lists_the_same(lista,listb)
#                                          - Simple function to compare two lists and whether lista and listb have the same length, and if all of
#                                            the elements of a are in b. (if they have the same length, and all elements in a are in b, then the
#                                            lists should be the same.
#*********************************************************************************************************************************************
# 09/15/2021    1.1      Lucas hunt      Include a parameter class in order to save information about each measurement set (important directories
#                                        Antenna names, reference antennas, calibration scans, etc.
#*********************************************************************************************************************************************
#*********************************************************************************************************************************************
# 09/15/2021    1.2      Lucas hunt      Add a couple of operators to find the frequency bands that are in the observation (i.e. S/X observations)
#                                        and create a dict that gives the spectral windows that are associated with each frequency band. This
#                                        should make it easier to split the data automatically if you have a different spw setup (i.e. RDV vs.
#                                        UF/UG/UH observations. 
#*********************************************************************************************************************************************
#*********************************************************************************************************************************************
# 12/1/2021     1.3      Lucas hunt      Implementing pickle instead of having to edit the init section of the class every time I want to add
#                                        a parameter. 
#*********************************************************************************************************************************************

import os
import numpy as np
import pyfits as fits
import time
from distutils.spawn import find_executable
import random
import csv
import shutil
import itertools
try:
#     #from casa import tbtool,msmdtool,tclean,fringefit,bandpass,rmtables,calstat,imstat,partition,split
     import casa as ctsks
     import casac
     msmd=ctsks.msmdtool()
     tb=ctsks.tbtool()
except:
     import casatools as ctls
     import casatasks as ctsks
     msmd=ctls.msmetadata()
     tb=ctls.table()  
from io import BytesIO as StringIO
import contextlib
import sys
#import casa




class cal_script_params:
     def __init__(self,band=''):
          self.wd=os.getcwd()
          self.experiment_name=self.wd.split('/')[-1]
          self.__cal_script_param_exists=False
          if os.path.exists(os.path.join(self.wd,'cal_script_param')):
               cal_script_dir_files=os.listdir(os.path.join(self.wd,'cal_script_param'))
               for fl in cal_script_dir_files:
                    if (not band) & (len(fl.split('.txt')[0].split('_'))==1):
                         self.__cal_script_param_exists=True 
                    if len(fl.split('.txt')[0].split('_'))>1:
                         if fl.split('.txt')[0].split('_')[-1].upper()==band.upper():
                              self.__cal_script_param_exists=True
          files=os.listdir(self.wd)
          for fl in files:
               if fl.endswith('.idifits'):
                    self.idifitsfile=os.path.join(self.wd,fl)
          if not band:
               self.ms_fl=self.experiment_name+'.ms'
          else:
               self.ms_fl=self.experiment_name+'_{}.ms'.format(band.upper())
          self.__make_important_directories()
          self.ms_file=os.path.join(self.exp_ms_dir,self.ms_fl)
          self.cal_file_root=os.path.join(self.first_cal_dir,self.ms_fl.split('.ms')[0])
          self.scans=[]
          self.scan_low_bsln=[]
          self.cal_scan_list=[]
          self.possible_cal_scans={}
          self.cal_sources_string_array=[]
          self.cal_sources_string=''
          self.cal_sources_array=[]
          self.all_antennas=[]
          self.missing_inner_antennas=[]
          self.missing_outer_antennas=[]
          self.inner_antennas_available=[]
          self.outer_antennas_available=[]
          self.num_antennas=0
          self.flag_spw_string=''
          self.fields=[]
          self.pflags_for_scan={}
          self.snr_sum_for_scan={}
          self.bandpass_snr_dict={}
          self.refant_selection_dict={}
          self.phasecal_scan=0
          self.nspw=0
          self.nchan=0
          self.bandpass_scan=0
          self.split_cal_ms_fils=[]
          self.split_ms_fils=[]
          self.spw_freq_dict={}
          hdu=fits.open(self.idifitsfile)
          freq_array=hdu['FREQUENCY'].data[0]['BANDFREQ']+hdu['FREQUENCY'].header['REF_FREQ']
          for i in range(len(freq_array)):
               band_key=find_band_letter(freq_array[i])
               if band_key not in self.spw_freq_dict:
                    self.spw_freq_dict[band_key]=[i]
               elif i not in self.spw_freq_dict[band_key]:
                    self.spw_freq_dict[band_key].append(i)
          hdu.close()
          self.need_split=False
          self.initial_fringefit_scans_string=''
          self.initial_fringefit_scans=[]
          self.final_fringefit_scans=[]
          self.final_fringefit_scans_string=''
          self.bad_scans=[]
          self.ant_number_name_dict={}
          self.refants=''
          self.apply_cal_tables_list=[]
          self.cal_fields=[]
          self.source_fields=[]
          self.source_cal_dict={}
          self.apply_cal_interp_list=[]
          self.apply_cal_spwmap_list=[]
#********************************************************************************************************       
     def __make_important_directories(self):
          self.wd=os.getcwd()
          self.cal_dir=os.path.join(self.wd,'cal_tables')
          if not os.path.exists(self.cal_dir):
               os.mkdir(self.cal_dir)
          
          self.experiment_cal_dir=os.path.join(self.cal_dir,self.ms_fl.split('.ms')[0])
          if not os.path.exists(self.experiment_cal_dir):
               os.mkdir(self.experiment_cal_dir)

          self.first_cal_dir=os.path.join(self.experiment_cal_dir,'first_cal')
          if not os.path.exists(self.first_cal_dir):
               os.mkdir(self.first_cal_dir)
          
          self.second_cal_dir=os.path.join(self.experiment_cal_dir,'second_cal')
          if not os.path.exists(self.second_cal_dir):
               os.mkdir(self.second_cal_dir)                       

          self.supplement_cal_dir=os.path.join(self.experiment_cal_dir,'supplement')
          if not os.path.exists(self.supplement_cal_dir):
               os.mkdir(self.supplement_cal_dir)
               
          self.gaincal_csv_dir=os.path.join(self.supplement_cal_dir,'antenna_gaincal_csvs')
          if not os.path.exists(self.gaincal_csv_dir):
               os.mkdir(self.gaincal_csv_dir)

          self.sn_csv_dir=os.path.join(self.supplement_cal_dir,'sn_csvs')
          if not os.path.exists(self.sn_csv_dir):
               os.mkdir(self.sn_csv_dir)
      
          self.im_dir=os.path.join(self.wd,'images')
          if not os.path.exists(self.im_dir):
               os.mkdir(self.im_dir)

          self.exp_im_dir=os.path.join(self.im_dir,self.ms_fl.split('.ms')[0])
          if not os.path.exists(self.exp_im_dir):
               os.mkdir(self.exp_im_dir)
               
          self.gain_cal_ims=os.path.join(self.exp_im_dir,'gain_cal')
          if not os.path.exists(self.gain_cal_ims):
               os.mkdir(self.gain_cal_ims)
         
          self.final_ims=os.path.join(self.exp_im_dir,'final_images')
          if not os.path.exists(self.final_ims):
               os.mkdir(self.final_ims)          
          
          self.ms_dir=os.path.join(self.wd,'ms_files')
          if not os.path.exists(self.ms_dir):
               os.mkdir(self.ms_dir)

          self.exp_ms_dir=os.path.join(self.ms_dir,self.ms_fl.split('.ms')[0])
          if not os.path.exists(self.exp_ms_dir):
               os.mkdir(self.exp_ms_dir)

          self.cal_ms_dir=os.path.join(self.exp_ms_dir,'cal_ms')
          if not os.path.exists(self.cal_ms_dir):
               os.mkdir(self.cal_ms_dir)
               
          self.split_ms_dir=os.path.join(self.exp_ms_dir,'split_ms')
          if not os.path.exists(self.split_ms_dir):
               os.mkdir(self.split_ms_dir)
          
          self.cal_script_params_dir=os.path.join(self.wd,'cal_script_param')
          if not os.path.exists(self.cal_script_params_dir):
               os.mkdir(self.cal_script_params_dir)
#********************************************************************************************************    
     def save_phasecal_scan(self,pcscn):
          self.phasecal_scan=pcscn
#********************************************************************************************************    
     def save_refants(self,refant):
         self.refants=refant
#********************************************************************************************************   
     def change_cal_file_root(self,new_CFR=''):
          if ('first_cal' in self.cal_file_root) & (not new_CFR):
               self.cal_file_root=os.path.join(self.second_cal_dir,self.ms_file.split('/')[-1].split('.')[0])
          elif new_CFR:
               self.cal_file_root=new_CFR
          else:
              print('already in second cal directory. Not going to change')
#********************************************************************************************************
     def get_metadata_twntfr_hr(self,phaseref=False):
          msmd.open(self.ms_file)
          
          for antid in msmd.antennaids():
               self.ant_number_name_dict[antid]=msmd.antennanames()[antid]
          
          self.scans=msmd.scannumbers()
          self.scn_low_bsln=[]
          if not phaseref:
               self.cal_scan_list=[]
               self.possible_cal_scans={}
               for scan in self.scans:
                    ants=msmd.antennasforscan(scan)
                    if len(ants)<5:
                         self.scn_low_bsln.append(scan)
                    if str(msmd.antennasforscan(scan))==str(msmd.antennaids()):
                         self.cal_scan_list.append(scan)
                    self.possible_cal_scans[scan]=msmd.antennasforscan(scan)
          else:
               fields=ctsks.vishead(vis=self.ms_file,mode='get',hdkey='field')[0]
               fld_codes=ctsks.vishead(vis=self.ms_file,mode='get',hdkey='fld_code')[0]
               self.cal_fields=[]
               self.source_fields=[]
               for val in range(0,len(fields)):
                    if fld_codes[val]=='V':
                         self.cal_fields.append(fields[val])
                    else:
                         self.source_fields.append(fields[val])
               if (not self.cal_fields) | (not self.source_fields):
                    self.cal_fields=[]
                    self.source_fields=[]
                    for i in range(len(fields)):
                         if i%2==0:
                              self.source_fields.append(fields[i])
                         else:
                              self.cal_fields.append(fields[i])
               self.source_cal_dict={}
               for fld in self.source_fields:
                    self.source_cal_dict[fld]=fields[msmd.fieldsforscan(msmd.scansforfield(fld)[0]+1)][0]
               self.possible_cal_scans={}
               self.cal_scan_list=[]
               for fld in self.cal_fields:
                    for scan in msmd.scansforfield(fld):
                         self.possible_cal_scans[scan]=msmd.antennasforscan(scan)
          
          num_ants_pcs=0
          if not self.cal_scan_list:
               for scan in list(self.possible_cal_scans):
                    if len(self.possible_cal_scans[scan])==num_ants_pcs:
                         self.cal_scan_list.append(scan)
                    if len(self.possible_cal_scans[scan])>num_ants_pcs:
                         self.cal_scan_list=[]
                         num_ants_pcs=len(self.possible_cal_scans[scan])
                         self.cal_scan_list.append(scan)

          tmp_cal_scan=self.cal_scan_list
          self.cal_scan_list=[]
          for i in range(len(tmp_cal_scan)):
               if i<65:
                    self.cal_scan_list.append(tmp_cal_scan[i])
          
          self.cal_sources_string_array=[]
          self.cal_sources_string=''
          self.cal_sources_array=[]
          if self.cal_fields:
               for fld in self.cal_fields:
                   self.cal_sources_string_array.append(fld)
                   self.cal_sources_array.append(msmd.fieldsforname(fld)[0])
                   self.cal_sources_string=self.cal_sources_string+fld+','
          else:
               for scn in self.cal_scan_list:
                    if msmd.namesforfields(msmd.fieldsforscan(scn)[0])[0] not in self.cal_sources_string_array:
                         self.cal_sources_array.append(msmd.fieldsforscan(scn)[0])
                         self.cal_sources_string_array.append(msmd.namesforfields(msmd.fieldsforscan(scn)[0])[0])
                         self.cal_sources_string=self.cal_sources_string+msmd.namesforfields(msmd.fieldsforscan(scn)[0])[0]+','

          self.cal_sources_string=self.cal_sources_string[:-1]

          self.cal_sources_scans=[]
          for source in self.cal_sources_array:
               for scan in msmd.scansforfield(source):
                    self.cal_sources_scans.append(scan)

          self.all_antennas=np.array(msmd.antennanames())
          inner_antennas=np.array(['LA','PT','OV','KP','FD'])
          outer_antennas=np.array(['BR','HN','NL','MK','SC'])
          self.missing_inner_antennas=np.setdiff1d(inner_antennas,self.all_antennas)
          self.missing_outer_antennas=np.setdiff1d(outer_antennas,self.all_antennas)
          self.inner_antennas_available=np.setdiff1d(inner_antennas,self.missing_inner_antennas)
          self.outer_antennas_available=np.setdiff1d(outer_antennas,self.missing_outer_antennas)
          for antenna in self.all_antennas:
               if (antenna not in self.inner_antennas_available) & (antenna not in self.outer_antennas_available):
                    self.outer_antennas_available=np.append(self.outer_antennas_available,antenna)

          
          self.bad_scans=[]
          for scan in self.possible_cal_scans:
               if len(self.possible_cal_scans[scan])<4:
                    self.bad_scans.append(scan)
               else:
                    vlba_antenna=False
                    for antenna in self.possible_cal_scans[scan]:
                         if (self.all_antennas[antenna] in self.inner_antennas_available) | (self.all_antennas[antenna] in self.outer_antennas_available):
                              vlba_antenna=True
                              break
                    if not vlba_antenna:
                         self.bad_scans.append(scan)

          self.initial_fringefit_scans=[]
          self.intiial_fringefit_scans_string=''
          self.final_fringefit_scans=[]
          self.final_fringefit_scans_string=''
          for scan in self.scans:
               if scan not in self.bad_scans:
                    self.final_fringefit_scans.append(scan)
                    self.final_fringefit_scans_string=self.final_fringefit_scans_string+str(scan)+','
                    if scan in self.cal_sources_scans:
                         self.initial_fringefit_scans.append(scan)
                         self.initial_fringefit_scans_string=self.initial_fringefit_scans_string+str(scan)+','

          self.initial_fringefit_scans_string=self.initial_fringefit_scans_string[:-1]
          self.final_fringefit_scans_string=self.final_fringefit_scans_string[:-1]


          self.num_antennas=len(msmd.antennaids())
          spws=''
          for sp in range(0,msmd.nspw()):
               spws=spws+'{specwin}:0~{low_flag_chan};{high_flag_chan}~{num_chans},'.format(specwin=sp,
                                                                                            low_flag_chan=int(msmd.nchan(sp)*0.12),
                                                                                            high_flag_chan=msmd.nchan(sp)-1-int(msmd.nchan(sp)*0.12),
                                                                                            num_chans=msmd.nchan(sp)-1)
                                                                                            
          self.nspw=msmd.nspw()
          self.nchan=msmd.nchan(0)                                                          
          self.flag_spw_string=spws[:-1]
          
          self.fields=msmd.fieldnames()
          
          msmd.close()
          self.need_split=(len(list(self.spw_freq_dict))>1)
#***********************************************************************************************************************************************          
     def find_strong_cal_source(self,cal_scan_list=[],num_scans=0,refants='LA,PT,FD,KP',slint='inf'):
          '''This function takes a measurement set (ms_file), a list of possible calibration scans (cal_scan_list), 
     	     the number of scans you would like it to process (default = 50) and a list of reference antennas. Then
      	     it will attempt to do a fringefit on each scan, and put the total signal to noise and the percentage 
             of scans that have been flagged from the calibration table into two dictionaries (pflags_for_scan 
             (flags) and snr_sum_for_scan (total snr).
          ''' 
          self.refant_selection_dict[refants]={}
          self.__refants=refants.split(',')
          self.__pflags_for_scan={}
          self.__snr_sum_for_scan={}
          scns_cmplt=0
          if not cal_scan_list:
               cal_scns=self.cal_scan_list
          else:
               cal_scns=cal_scan_list
          if not num_scans:
               nm_scns=len(cal_scns)
          else:
               nm_scns=num_scans
          for scn in cal_scns:
               self.__valid_refant=False
               for ant in self.possible_cal_scans[scn]:
                    if self.ant_number_name_dict[ant] in self.__refants:
                         tb.open(self.ms_file)
                         antenna1=tb.getcol('ANTENNA1')
                         antenna2=tb.getcol('ANTENNA2')
                         flags=tb.getcol('FLAG')
                         if flags.ndim==1:
                              dim1=len(flags)/len(antenna1)/self.nchan
                              dim2=self.nchan
                              flags=flags.reshape((dim1,dim2,-1))
                         scans=tb.getcol('SCAN_NUMBER')
                         total_flags=0
                         total_datapoints=0
                         for i in range(len(flags)):
                              for j in range(len(flags[i])):
                                   total_flags+=float(sum(flags[i][j][np.where(((antenna1==ant) | (antenna2==ant))& (scans==scn))]))
                                   total_datapoints+=float(len(flags[i][j][np.where(((antenna1==ant) | (antenna2==ant))& (scans==scn))]))
                         percent_flagged=total_flags/total_datapoints
                         if percent_flagged<0.9:
                              self.__valid_refant=True
                              break
               if scns_cmplt<nm_scns:
                    if self.__valid_refant:
                         ctsks.fringefit(vis=self.ms_file,
                                   caltable=self.cal_file_root+'_test.sbd',
                                   scan=str(scn),
                                   solint=slint,
                                   zerorates=True,
                                   refant=refants,
                                   minsnr=10.0,
                                   parang=True)
                         try:
                              tb.open(self.cal_file_root+'_test.sbd')
                              self.__pflags_for_scan[scn]=float(sum(tb.getcol('FLAG')[0][0]))/float(len(tb.getcol('FLAG')[0][0]))
                              self.__snr_sum_for_scan[scn]=sum(tb.getcol('SNR')[0][0])
                              tb.close()
                         except:
                              self.__pflags_for_scan[scn]=1.0
                              self.__snr_sum_for_scan[scn]=0.0
                         ctsks.rmtables(self.cal_file_root+'_test.sbd')
                         scns_cmplt+=1
          self.refant_selection_dict[refants]['pflags']=self.__pflags_for_scan
          self.refant_selection_dict[refants]['snr']=self.__snr_sum_for_scan       
#***********************************************************************************************************************************************          
     def find_bandpass_scan(self,cal_scan_list=[],refants='LA,PT,FD,KP',gaincal_corr=[]):
          self.bandpass_snr_dict={}
          if not cal_scan_list:
               cal_scns=self.cal_scan_list
          else:
               cal_scns=cal_scan_list
          for scn in cal_scns:
               if gaincal_corr:
                    ctsks.bandpass(vis=self.ms_file,
                                   caltable=self.cal_file_root+'_test.bpass',
                                   gaintable=gaincal_corr,
                                   scan=str(scn),
                                   solnorm=True,
                                   solint='inf',
                                   refant=refants,
                                   bandtype='B',
                                   parang=True)
               else:
                    ctsks.bandpass(vis=self.ms_file,
                                   caltable=self.cal_file_root+'_test.bpass',
                                   scan=str(scn),
                                   solnorm=True,
                                   solint='inf',
                                   refant=refants,
                                   bandtype='B',
                                   parang=True)
               self.bandpass_snr_dict[scn]=ctsks.calstat(caltable=self.cal_file_root+'_test.bpass',axis='snr')['SNR']['mean']
          top_bp_snr_scan=list(self.bandpass_snr_dict)[0]
          for scn in self.bandpass_snr_dict:
               if self.bandpass_snr_dict[scn]>self.bandpass_snr_dict[top_bp_snr_scan]:
                    top_bp_snr_scan=scn
          self.bandpass_scan=top_bp_snr_scan         
#***********************************************************************************************************************************************          
     def split_cal_files(self,band=''):
          if not band:
               band=self.ms_file.split('/')[-1].split('_')[-1].split('.')[0]
          self.split_cal_ms_fils=[]
          if len(band)==1:
               bnd_strng='_'+band
          else:
               bnd_strng=''
          for fld in self.cal_sources_string_array:
               ctsks.split(vis=self.ms_file,
                           datacolumn='corrected',
                           field=fld, 
                           outputvis=os.path.join(self.cal_ms_dir,fld+'_gc{}.ms'.format(bnd_strng)))
               self.split_cal_ms_fils.append(os.path.join(self.cal_ms_dir,fld+'_gc{}.ms'.format(bnd_strng)))
#***********************************************************************************************************************************************          
     def split_ms_files(self,band=''):
          if not band:
               band=self.ms_file.split('/')[-1].split('_')[-1].split('.')[0]
          print(band)
          self.split_ms_fils=[]
          if len(band)==1:
               bnd_strng='_'+band
          else:
               bnd_strng=''
          print(bnd_strng)
          for fld in self.fields:
               ctsks.split(vis=str(self.ms_file),
                           datacolumn='corrected',
                           field=str(fld), 
                           outputvis=str(os.path.join(self.split_ms_dir,fld+'{}.ms'.format(bnd_strng))))
               self.split_ms_fils.append(os.path.join(self.split_ms_dir,fld+'{}.ms'.format(bnd_strng)))
#***********************************************************************************************************************************************          
     def phasecal_split_ms_files(self,band=''):
          if not band:
               band=self.ms_file.split('/')[-1].split('_')[-1].split('.')[0]
          self.split_ms_fils=[]
          for fld in self.source_fields:
               if os.path.exists(os.path.join(self.split_ms_dir,'{source}_{cal}_{bnd}.ms'.format(source=fld,cal=self.source_cal_dict[fld],bnd=band))):
                    ctsks.rmtables(tablenames=os.path.join(self.split_ms_dir,'{source}_{cal}_{bnd}.ms'.format(source=fld,cal=self.source_cal_dict[fld],bnd=band)))
               if os.path.exists(os.path.join(self.split_ms_dir,'{source}_{cal}_{bnd}.ms'.format(source=fld,cal=self.source_cal_dict[fld],bnd=band))+'.flagversions'):
                    shutil.rmtree(os.path.join(self.split_ms_dir,'{source}_{cal}_{bnd}.ms'.format(source=fld,cal=self.source_cal_dict[fld],bnd=band))+'.flagversions')
               ctsks.split(vis=self.ms_file,
                           datacolumn='corrected',
                           outputvis=os.path.join(self.split_ms_dir,'{source}_{cal}_{bnd}.ms'.format(source=fld,cal=self.source_cal_dict[fld],bnd=band)),
                           field='{source},{cal}'.format(source=fld,cal=self.source_cal_dict[fld]))
               self.split_ms_fils.append(os.path.join(self.split_ms_dir,'{source}_{cal}_{bnd}.ms'.format(source=fld,cal=self.source_cal_dict[fld],bnd=band)))

#***********************************************************************************************************************************************          
#***********************************************************************************************************************************************          
#***********************************************************************************************************************************************          

#***************************************
# Uses numpy to determine a 4th order
# polynomial from the values for
# gaincurve response
#***************************************

def transform_poly(coeff, min_elev=0, max_elev=90):
    f = np.poly1d(coeff[::-1])
    g = lambda x: np.sqrt(f(90 - x))
    x = np.linspace(min_elev, max_elev, 64, endpoint=True)
    y = g(x)
    return np.poly1d(np.polyfit(x, y, 3))


#***************************************
# VLBA frequency ranges to find what
# band the observations were taken in
# needs to be expanded to include bands
# other than s/x
#***************************************

def find_band_letter(freq):
    if freq>1.0e9 and freq<2.0e9:
        return 'L'
    elif freq>2.2e9 and freq<2.5e9:
        return 'S'
    elif freq>7.95e9 and freq<12.0e9:
        return 'X'
    elif freq>3.5e9 and freq<7.95e9:
        return 'C'
    elif freq>12.0e9 and freq<18.0e9:
        return 'Ku'
    elif freq>18.0e9 and freq<26.5e9:
        return 'K'
    elif freq>26.5e9 and freq<40.0e9:
        return 'Ka'
    elif freq>40.0e9 and freq<75.0e9:
        return 'Q'
    elif freq>75.0e9 and freq<110.0e9:
        return 'W'

#***************************************
# Makes dictionary for all antennas incl
# in the observation
#***************************************

def make_ant_dict(hdu):
    ant_dict={}
    for i in range(0,len(hdu['antenna'].data['anname'])):
        ant_dict[hdu['antenna'].data['antenna_no'][i]]=hdu['antenna'].data['anname'][i]
    return ant_dict


#Outline:
#*******************************************************************   
# 1. Define useful functions. These will be used at various points 
# in the pipeline. Might as well define them now!
#*******************************************************************

def find_start_end_times(listobs_file):
     mo_dict={'jan':'1',
              'feb':'2',
              'mar':'3',
              'apr':'4',
              'may':'5',
              'jun':'6',
              'jul':'7',
              'aug':'8',
              'sep':'9',
              'oct':'10',
              'nov':'11',
              'dec':'12'}
     start_time_dict={}
     end_time_dict={}
     start_date=''
     with open(listobs_file , 'r') as fl: 
          start_recording_scan_time=10000
          stop_recording_scan_time=10000
          for cnt,line in enumerate(fl):
               split_line=line.split()
               if (split_line):
                    if (split_line[0]=='Date'): 
                         start_recording_scan_time=cnt+1 
                    if (split_line[0]=='(nRows'): 
                         stop_recording_scan_time=cnt
                    if (split_line[0]=='Observed'):
                         line_start_date=split_line[2].split('/')[0]
                         start_date_list=line_start_date.split('-')
                         start_date_list[1]=mo_dict[start_date_list[1].lower()]
                         start_date='{year}-{month}-{day}'.format(year=start_date_list[2],month=start_date_list[1],day=start_date_list[0])
               if (cnt>start_recording_scan_time) & (cnt<stop_recording_scan_time): 
                    start_time_dict[int(split_line[3])]=split_line[0]
                    end_time_dict[int(split_line[3])]=split_line[2]
     return start_date,start_time_dict,end_time_dict

def find_bandpass_scan(ms_file,cal_file_root,cal_time_dict,refants='LA,PT,FD,KP',gaincal_corr=''):
     bandpass_snr_dict={}
     for scn in cal_time_dict:
          if gaincal_corr:
               ctsks.bandpass(vis=ms_file,
                              caltable='bandpass_test.bpass',
                              gaintable=[gaincal_corr],
                              scan=str(scn),
                              solnorm=True,
                              solint='inf',
                              refant=refants,
                              bandtype='B',
                              parang=True)
          else:
               ctsks.bandpass(vis=ms_file,
                              caltable='bandpass_test.bpass',
                              scan=str(scn),
                              solnorm=True,
                              solint='inf',
                              refant=refants,
                              bandtype='B',
                              parang=True)
          bandpass_snr_dict[scn]=ctsks.calstat(caltable='bandpass_test.bpass',axis='snr')['SNR']['mean']
     return bandpass_snr_dict

def calibrator_dict(ms_file,cal_scan_list,num_scans=50,refants='LA,PT,FD,KP'):
     '''This function takes a measurement set (ms_file), a list of possible calibration scans (cal_scan_list), 
        the number of scans you would like it to process (default = 50) and a list of reference antennas. Then
        it will attempt to do a fringefit on each scan, and put the total signal to noise and the percentage 
        of scans that have been flagged from the calibration table into two dictionaries (pflags_for_scan 
        (flags) and snr_sum_for_scan (total snr).
     ''' 
     pflags_for_scan={}
     snr_sum_for_scan={}
     scns_cmplt=0
     for scn in cal_scan_list:
          if scns_cmplt<num_scans:
               ctsks.fringefit(vis=ms_file,
                               caltable=cal_file_root+'_test.sbd',
                               scan=str(scn),
                               solint='inf',
                               zerorates=True,
                               refant=refants,
                               minsnr=10.0,
                               parang=True)
               try:
                    tb.open(cal_file_root+'_test.sbd')
                    total_flags=0
                    total_datapoints=0
                    total_snr=0
                    flgs=tb.getcol('FLAG')
                    for i in len(flgs):
                         for j in len(flgs[j]):
                              total_flags+=float(sum(flgs[i][j]))
                              total_datapoints+=float(len(flgs[i][j]))
                              total_snr+=float(sum(tb.getcol('SNR')[i][j]))
                    pflags_for_scan[scn]=total_flags/total_datapoints
                    snr_sum_for_scan[scn]=total_snr
                    tb.close()
               except:
                    pflags_for_scan[scn]=1.0
                    snr_sum_for_scan[scn]=0.0
               ctsks.rmtables(cal_file_root+'_test.sbd')
               scns_cmplt+=1
     return pflags_for_scan,snr_sum_for_scan

def write_cal_dict(pflag_dict={},snr_dict={},outfile='Printed_Dict.csv'):
     ''' Takes the dictionaries from above and writes them as text files. This can be used later for diagnostic
         purposes
     '''
     if pflag_dict:
          with open(outfile,'w') as pfl:
               pfl.write('Scan,SNR,Flags\n')
               for scn in snr_dict:
                    pfl.write('{scan},{snr},{flg}\n'.format(scan=scn,snr=snr_dict[scn],flg=pflag_dict[scn]))
     else:
          with open(outfile,'w') as pfl:
               pfl.write('Scan,SNR\n')
               for scn in snr_dict:
                    pfl.write('{scan},{snr}\n'.format(scan=scn,snr=snr_dict[scn]))

def make_refant_list(inner_ants,outer_ants,initial_ant):
     ''' This will take all of the VLBA antennas (since we know which antennas are in the inner part and which are
         in the outer part), and returns a list of reference anetnnas. The order is random (except for the first 
         antenna from the list, which should be predetermined (initial_ant)) with the five inner antennas always 
         being the first five, and the five outer antennas being the last five
     '''
     refant=initial_ant
     inner_ants=np.setdiff1d(inner_ants,np.array(initial_ant.split(',')))
     for i in reversed(range(len(inner_ants))):
          rand_int=random.randint(0,i)
          refant=refant+','+inner_ants[rand_int]
          inner_ants=inner_ants[inner_ants!=inner_ants[rand_int]]
     for i in reversed(range(len(outer_ants))):
          rand_int=random.randint(0,i)
          refant=refant+','+outer_ants[rand_int]
          outer_ants=outer_ants[outer_ants!=outer_ants[rand_int]]
     print(refant)
     return refant

def find_bandpass_scan(ms_fl,cl_fle_rt,cl_scns,refants='LA,PT,FD,KP'):
     bandpass_snr_dict={}
     for scn in cl_scns:
          ctsks.bandpass(vis=ms_fl,
                         caltable=cl_fle_rt+'_bandpass_test.bpass',
                         scan=str(scn),
                         solnorm=True,
                         solint='inf',
                         refant=refants,
                         bandtype='B',
                         parang=True)
          bandpass_snr_dict[scn]=ctsks.calstat(cl_fle_rt+'_bandpass_test.bpass',axis='snr')['SNR']['mean']
          ctsks.rmtables(os.path.join(cal_file_root,'bandpass_test.bpass'))
     return bandpass_snr_dict
     
def get_unimaged(img_dir='',splt=''):
     imgs=[f.split('_X')[0] for f in os.listdir(img_dir) if '_X' in f]
     srcs=[f.split('.ms')[0] for f in os.listdir(splt) if '.flag' not in f]
     imgs=np.array(imgs)
     srcs=np.array(srcs)
     unimaged=np.setdiff1d(srcs,imgs)
     unimaged_ms=[f+'.ms' for f in unimaged]
     return unimaged_ms
     
def calc_pixel_size(msin):
#     baseline_length_dict={'SC':{'HN':2853000,'NL':3645000,'FD':4143000,'LA':4458000,'PT':4579000,'KP':4839000,'OV':5460000,'BR':5767000,'MK':8611000},
#                           'HN':{'SC':2853000,'NL':1611000,'FD':3105000,'LA':3006000,'PT':3226000,'KP':3623000,'OV':3885000,'BR':3657000,'MK':7502000},
#                           'NL':{'SC':3645000,'HN':1611000,'FD':1654000,'LA':1432000,'PT':1663000,'KP':2075000,'OV':2328000,'BR':2300000,'MK':6156000},
#                           'FD':{'SC':4143000,'HN':3105000,'NL':1654000,'LA':608000,'PT':564000,'KP':744000,'OV':1508000,'BR':2345000,'MK':5134000},
#                           'LA':{'SC':4458000,'HN':3006000,'NL':1432000,'FD':608000,'PT':236000,'KP':652000,'OV':1088000,'BR':1757000,'MK':4970000},
#                           'PT':{'SC':4579000,'HN':3226000,'NL':1663000,'FD':564000,'LA':236000,'KP':417000,'OV':973000,'BR':1806000,'MK':4795000},
#                           'KP':{'SC':4839000,'HN':3623000,'NL':2075000,'FD':744000,'LA':652000,'PT':417000,'OV':845000,'BR':1913000,'MK':4466000},
#                           'OV':{'SC':5460000,'HN':3885000,'NL':2328000,'FD':1508000,'LA':1088000,'PT':973000,'KP':845000,'BR':1214000,'MK':4015000},
#                           'BR':{'SC':5767000,'HN':3657000,'NL':2300000,'FD':2345000,'LA':1757000,'PT':1806000,'KP':1913000,'OV':1214000,'MK':4398000},
#                           'MK':{'SC':8611000,'HN':7502000,'NL':6156000,'FD':5134000,'LA':4970000,'PT':4795000,'KP':4466000,'OV':4015000,'BR':4398000}}
     tb.open(msin+'/ANTENNA')
     positions=tb.getcol('POSITION')
     stations=tb.getcol('STATION')
     tb.close()
     antenna_geocentric_coord_dict={}
     baseline_length_dict={}
     for ant in range(len(stations)):                                                                                           
          antenna_geocentric_coord_dict[stations[ant]]={}
          antenna_geocentric_coord_dict[stations[ant]]['X']=positions[0][ant]
          antenna_geocentric_coord_dict[stations[ant]]['Y']=positions[1][ant]
          antenna_geocentric_coord_dict[stations[ant]]['Z']=positions[2][ant]     
     for pair in itertools.permutations(list(antenna_geocentric_coord_dict),2):
          if pair[0] not in list(baseline_length_dict): 
               baseline_length_dict[pair[0]]={}      
          baseline_length_dict[pair[0]][pair[1]]=np.sqrt(np.abs(antenna_geocentric_coord_dict[pair[0]]['X']-antenna_geocentric_coord_dict[pair[1]]['X'])**2+np.abs(antenna_geocentric_coord_dict[pair[0]]['Y']-antenna_geocentric_coord_dict[pair[1]]['Y'])**2+np.abs(antenna_geocentric_coord_dict[pair[0]]['Z']-antenna_geocentric_coord_dict[pair[1]]['Z'])**2)  
     msmd.open(msin)
     spws=msmd.spwsforscan(msmd.scannumbers()[0])
     mean_freq_spw_selector=int(len(spws)/2)
     frequency=msmd.meanfreq(spws[mean_freq_spw_selector])
     wavelength_m=3e8/frequency
     max_baseline=0
     for pair in itertools.permutations(msmd.antennanames(),2):
        if baseline_length_dict[pair[0]][pair[1]]>max_baseline:
             max_baseline=baseline_length_dict[pair[0]][pair[1]]
     hpbw_rad=wavelength_m/max_baseline
     hpbw_deg=hpbw_rad*180/np.pi
     hpbw_mas=hpbw_deg*3600*1000
     pixel_size_mas=hpbw_mas/5.
     msmd.close()
     return '{pxl_size}mas'.format(pxl_size=np.round(pixel_size_mas,3))

def image_loop(ms_file, im_name, im_size=512, robust_weight=2.0, clean_to=4.0, start_high_mask_threshold=8.0, end_high_mask_threshold=5.0, low_mask_threshold=10.0, prllel=False,n_iter=4000,mask_type='auto-multithresh',sdlobthresh=0.7,fieldname=''):
     no_mask_created=False
     cell_size=calc_pixel_size(ms_file)
     model_rms=0
     tmp_noisethreshold=start_high_mask_threshold
     counter=0
     while model_rms == 0:
          os.system('rm -rf {all_im_prods}*'.format(all_im_prods=im_name))
          print('*********Removed files {} times*******'.format(counter+1))
          #casalog.post('*********Removed files {} times*******'.format(counter+1))
          ctsks.tclean(vis=ms_file,
                 field=fieldname,
                 datacolumn='corrected',
                 imagename=im_name,
                 imsize=im_size,
                 cell=cell_size,
                 deconvolver='hogbom',
                 weighting='briggs',
                 robust=robust_weight,
                 niter=n_iter,
                 gain=0.1,
                 nsigma=clean_to,
                 usemask=mask_type,
                 sidelobethreshold=sdlobthresh,
                 noisethreshold=tmp_noisethreshold,
                 lownoisethreshold=low_mask_threshold,
                 restart=True,
                 savemodel='modelcolumn',
                 parallel=prllel)
          print('**********Finished tclean {} times*********'.format(counter+1))
          #casalog.post('**********Finished tclean {} times*********'.format(counter+1))
          tmp_imstat_dict=ctsks.imstat(imagename=im_name+'.model')
          model_rms=tmp_imstat_dict['rms'][0]
          print('***********Finished imstat and model rms. Found rms = {}*********'.format(model_rms))
          #casalog.post('***********Finished imstat and model rms. Found rms = {}*********'.format(model_rms))
          if model_rms == 0:
               print('***********The masking threshold is {msk}, the cleaning limit is {nth} sigma***********'.format(msk=tmp_noisethreshold,nth=clean_to))
               #casalog.post('The masking threshold is {msk}, the cleaning limit is {nth} sigma'.format(msk=tmp_noisethreshold,nth=clean_to))
               if tmp_noisethreshold>end_high_mask_threshold:
                    tmp_noisethreshold=tmp_noisethreshold-1.0
                    low_mask_threshold=tmp_noisethreshold*0.6666
               if tmp_noisethreshold<(clean_to+2.0):
                    clean_to=clean_to-1.0
               if tmp_noisethreshold<=end_high_mask_threshold:
                    model_rms=1
               counter+=1
     #if (tmp_noisethreshold>=end_high_mask_threshold) & (model_rms>0):
          #casalog.post('Image made for source {src} using noisethreshold {nth}'.format(src=im_name.split('/')[-1],nth=tmp_noisethreshold))
     #if (tmp_noisethreshold==end_high_mask_threshold) & (model_rms==0):
          #casalog.post('Possible no image made for source {src} because masking threshold is at limit of  {nth} sigma'.format(src=im_name.split('/')[-1],nth=tmp_noisethreshold))
     #     no_mask_created=True
     return tmp_noisethreshold,low_mask_threshold,no_mask_created

def are_lists_the_same(lista,listb):
     all_elements_the_same=True
     if len(lista)==len(listb):
          for elem in lista:import os


def test_tec_maps():
     fls=os.listdir('.')
     gnss_tec_file_exists=False
     for fl in fls:
         if gnss_tec_file_exists:
             return True
         if ('igsg' in fl) or ('igrg' in fl):
             gnss_tec_file_exists=True
     return gnss_tec_file_exists

@contextlib.contextmanager
def redirect_stdout(target):
    original = sys.stdout
    try:
        sys.stdout = target
        yield
    finally:
        sys.stdout = original

def add_fake_tsys_rdv(idifitsfile=''):
     if not idifitsfile:
          print('Need to define idifits file!')
     vlba_ants = ['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC'] # Need list of VLBA antennas to find non-vlba antennas. These will be 
     hdu=fits.open(idifitsfile)
     f=StringIO()
     with redirect_stdout(f):
         hdu.info()
     s=f.getvalue()
     hdu_index_name_dict={}
     for line in s.split('\n'):
          if line:
               try:
                    int(line.split()[0])
                    hdu_index_name_dict[int(line.split()[0])]=line.split()[1]
               except:
                    continue
     ANT_NUMBER_NAME_DICT={}
     for i in range(len(hdu['ANTENNA'].data)):
          ANT_NUMBER_NAME_DICT[hdu['ANTENNA'].data[i]['ANTENNA_NO']]=hdu['ANTENNA'].data[i]['ANNAME']
     tsys_ants=np.unique(hdu['SYSTEM_TEMPERATURE'].data['ANTENNA_NO'])
     non_tsys_ants=np.setdiff1d(np.array(list(ANT_NUMBER_NAME_DICT)),tsys_ants)
     default_tsys_ant=0
     for nm in ANT_NUMBER_NAME_DICT:
          if nm in tsys_ants:
               default_tsys_ant=nm
               break
     idx=np.where(hdu['SYSTEM_TEMPERATURE'].data['ANTENNA_NO']==default_tsys_ant)[0]
     tmp_tsys_table_dict={}
     i=0
     for cols in hdu['SYSTEM_TEMPERATURE'].columns:
          tmp_tsys_table_dict[cols.name]={}
          tmp_tsys_table_dict[cols.name]['DATA']=hdu['SYSTEM_TEMPERATURE'].data[cols.name]
          tmp_tsys_table_dict[cols.name]['FORMAT']=cols.format
          tmp_tsys_table_dict[cols.name]['UNIT']=cols.unit
          i+=1
     for antnum in non_tsys_ants:
          for i in range(len(hdu['SYSTEM_TEMPERATURE'].data[idx])):
              for cols in hdu['SYSTEM_TEMPERATURE'].columns:
                    if cols.name!='ANTENNA_NO':
                         tmp=np.append(tmp_tsys_table_dict[cols.name]['DATA'],hdu['SYSTEM_TEMPERATURE'].data[idx][i][cols.name])
                         tmp_tsys_table_dict[cols.name]['DATA'] = tmp
                    else:
                         tmp=np.append(tmp_tsys_table_dict[cols.name]['DATA'],antnum)
                         tmp_tsys_table_dict[cols.name]['DATA']=tmp
     tmp_tsys_table_dict['TSYS_1']['DATA']=np.reshape(tmp_tsys_table_dict['TSYS_1']['DATA'],(-1,len(hdu['SYSTEM_TEMPERATURE'].data['TSYS_1'][0])))
     tmp_tsys_table_dict['TANT_1']['DATA']=np.reshape(tmp_tsys_table_dict['TANT_1']['DATA'],(-1,len(hdu['SYSTEM_TEMPERATURE'].data['TANT_1'][0])))
     new_tsys_col_dict={}
     for col_name in tmp_tsys_table_dict:
          new_tsys_col_dict[col_name]=fits.Column(name=col_name,
                                                  array=tmp_tsys_table_dict[col_name]['DATA'],
                                                  format=tmp_tsys_table_dict[col_name]['FORMAT'],
                                                  unit=tmp_tsys_table_dict[col_name]['UNIT'])
     tsys_col_array=[new_tsys_col_dict[f] for f in list(new_tsys_col_dict)]
     organized_tsys_col_array=[]
     for col in hdu['SYSTEM_TEMPERATURE'].header['TTYPE*']:
          match_val=hdu['SYSTEM_TEMPERATURE'].header[col]
          for row in tsys_col_array:
               if match_val == row.name:
                    organized_tsys_col_array.append(row)
     tsys_hdu=fits.BinTableHDU.from_columns(organized_tsys_col_array,name='SYSTEM_TEMPERATURE')
     new_idifits_hdu=fits.HDUList([])
     for indx in hdu_index_name_dict:
          if hdu_index_name_dict[indx]!='SYSTEM_TEMPERATURE':
               new_idifits_hdu.append(hdu[hdu_index_name_dict[indx]])
          else:
               print('Appending System Temperature HDU')
               new_idifits_hdu.append(tsys_hdu)
     for val in list(hdu['SYSTEM_TEMPERATURE'].header):
          if val not in list(new_idifits_hdu['SYSTEM_TEMPERATURE'].header):
               new_idifits_hdu['SYSTEM_TEMPERATURE'].header[val]=hdu['SYSTEM_TEMPERATURE'].header[val]
     new_idifits_hdu.writeto(idifitsfile,clobber=True)
     hdu.close()
     new_idifits_hdu.close()

