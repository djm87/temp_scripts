# -*- coding: utf-8 -*-
"""
Extract time from ge.par and meta data
"""
import os
import glob
import pandas as pd
import json
import numpy as np

# folder_list = ['V5-7_insitu']
folder_list = ['V5-7_insitu','W1-2_insitu', 'V5-8_insitu', 'v7-1_insitu', 'V7-2_insitu', 'V7-3_insitu','V7-4_insitu', 
              'V7-7_insitu','V3-3_insitu', 'V3-6_insitu', 'V7-8_insitu', 'V5-2_insitu']

def extract_time(folder: str) -> pd.core.frame.DataFrame:

    ''' 
    par and json file are created during each scan 
    the json file contains the header file for the par file
    headers: date, time, epoch, SCAN_N, GE filenumber, omega, temperature, BF_DF_flag
    time from par file is likely the creation time 
    temperature reading is likely from the pyrometer 
    
    extract the modification time of the GE files
    the modification time is likely when the last scan is written to the GE file
    '''
    # import in par and json file
    par_file = glob.glob("/../../1_Raw_Data/CHESS_lados-3271-b/"+folder+"/*.par")
    json_file = glob.glob("/../../1_Raw_Data/CHESS_lados-3271-b/"+folder+"/*.json")
    header_names = json.loads(open(json_file[0], "r").readlines()[0])
    par_df = pd.read_csv(par_file[0], sep = " ", names=list(header_names.values()))
    par_df=par_df.rename(columns={"epoch": "Time_creat", "GE filenumber": "GE_filenumber", "temperature": "Temperature_pyro"})
    
    # extract modification time of the ge2 files
    ge2_ims = glob.glob("/../../1_Raw_Data/CHESS_lados-3271-b/"+folder+"_copy/*.ge2")+glob.glob("/../../1_Raw_Data/CHESS_lados-3271-b/"+folder+"_copy/*/*.ge2")
    mod_time_df = pd.DataFrame({'GE_filenumber':[],'Time_mod':[]}) 
    
    for im in ge2_ims:
        im_name = os.path.basename(im)
        im_index = int(im_name.split('_')[1].split('.')[0]) 
        mod_time = os.path.getmtime(im) 
        mod_time_df_dict = {'GE_filenumber':(im_index), 'Time_mod':mod_time} 
        mod_time_df = mod_time_df.append(mod_time_df_dict, ignore_index=True)
    
    mod_time_df.sort_values(by='GE_filenumber')
    
    # merge two dataframes based on the GE filenumber
    merged_df = pd.merge(par_df, mod_time_df, how='outer', on='GE_filenumber').sort_values(by='GE_filenumber').reset_index(drop=True)
    merged_df['Time_diff'] = merged_df['Time_mod']-merged_df['Time_creat']
    
    # print out the average difference between file creation and modification
    # the average difference is around 1/exposure time 
    avg_time_diff = np.average(merged_df['Time_diff'][merged_df['Time_diff'].notnull()])
    print("average difference between scan creation time and modified time is for "+folder+" is "
          +str(round(avg_time_diff, 3))+' s')
    
    merged_df.to_csv('./'+folder+'/'+folder+'_time.csv')
    return(merged_df) 
    
    
for folder in folder_list:
    globals()[f'{folder}_time_df']=extract_time(folder)