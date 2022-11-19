#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perform sequential Rietveld refinement for Ti64 + Steel
"""
import os
import sys
#sys.path.insert(0, '/home/bonnie/anaconda3/envs/GSASII/gsasii/GSASII')
sys.path.insert(0, 'C://Users/sanoj737/anaconda3/envs/G2/GSASII')
import GSASIIscriptable as G2sc
import glob
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [atoi(c) for c in re.split('(\d+)', text)]

def DataWrap(fil): return os.path.join("./", fil)

def PathWrap(fil): return os.path.join("./", fil)


# folder and gpx output names
folder_name = 'V5-7_insitu'
base_fname = '_int1_dc'
gpx_name_in = folder_name+'_single'+base_fname+'.gpx'
gpx_name_out = gpx_name_in.replace('_single', '_seq')

# import in desired histograms
search_string = folder_name+base_fname+'/*.fxye'
hists = []

bad_scans = [4764]
good_scans = [4330, 4769]

for hist in glob.glob(search_string):
    hist_name = os.path.basename(hist)
    scan_num = int(hist_name.split('_')[1].split('.')[0])
    if scan_num not in bad_scans and scan_num >= good_scans[0] and scan_num <= good_scans[1]:
        hists.append(hist)

hists.sort(key=natural_keys)
hists.pop(0)

diff = {"data": hists,
        "prm": "../../1_Calibration/Instrument_july13_1int.instprm",
        "fmthint": "GSAS powder",
        "databank": 1}

# define the key points
start_num = 4329 # start scan number
beta_app_num = 4422 # beta appear number
beta_str_num = 4432 # beta strong number
alpha_disapp_num = 4445 # alpha disappear number
alpha_reapp_num = 4483 # alpha reappear number
iso_hold_end_num = 4661 # end of isothermal hold number
end_num = 4769 # end scan number 


def main():

    gpx = G2sc.G2Project(gpxfile=PathWrap(gpx_name_in))  # load a project
    gpx.save(PathWrap(gpx_name_out))
    
    # extract scan number based on histogram name
    def num_from_hist_name(hist_name):
        return int(hist_name.split('ff_')[1].split('_')[0])
    
    # create a list of G2PwdrData histogram objects based on histogram names
    def gpx_hist(hist_names): 
        return [gpx.histogram(hist_name) for hist_name in hist_names]
    
    # create a list of G2PwdrData histogram objects based on scan numbers
    def gpx_hist_from_num(num1, num2):
        hist_names = []
        for j, hist in enumerate(gpx.histograms()):
            hist_name = hist['data'][-1]
            scan_num = int(hist_name.split('ff_')[1].split('_')[0])
            if scan_num >= num1 and scan_num < num2:
                hist_names.append(hist_name)
        return gpx_hist(hist_names)
    
    # set histogram/phase entry for scan numbers within a range
    def HAP_entry_for_num(entry_list, num1, num2):
        for i, p in enumerate(gpx.phases()):
            for j, (sclEnt, hstrainEnt, preforiEnt, mustrainEnt) in enumerate(
                    zip(p.getHAPentryList(0, 'Scale'),
                        p.getHAPentryList(0, 'HStrain'),
                        p.getHAPentryList(0, 'Pref.Ori.'),
                        p.getHAPentryList(0, 'Mustrain'))):
                hist_name = sclEnt[0][0]
                scan_num = num_from_hist_name(hist_name)
                if scan_num >= num1 and scan_num < num2:
                    p.setHAPentryValue(sclEnt[0], entry_list[i]['Scale'])
                    p.setHAPentryValue(hstrainEnt[0], entry_list[i]['HStrain'])
                    p.setHAPentryValue(preforiEnt[0], entry_list[i]['Pref.Ori.'])
                    p.setHAPentryValue(mustrainEnt[0], entry_list[i]['Mustrain'])
    
    # calculate alpha and beta phase fraction relative to steel phase fraction
    # use when there's NO total phase fraction sum to 1 constraint
    def phafrac_from_wgtfrac(steel_wgtfrac, alpha_wgtfrac, steel_phafrac = 1):
        # calculate phase fraction of alpha and beta given that steel phase fraction is set to 1
        steel_molar_mass = gpx.phase(0)['General']['Mass'] 
        alpha_molar_mass = gpx.phase(1)['General']['Mass']
        beta_molar_mass = gpx.phase(2)['General']['Mass']
        
        beta_wgtfrac = 1-steel_wgtfrac-alpha_wgtfrac
        
        beta_phafrac = beta_wgtfrac*steel_phafrac*steel_molar_mass/(steel_wgtfrac*beta_molar_mass)
        alpha_phafrac = alpha_wgtfrac*steel_phafrac*steel_molar_mass/(steel_wgtfrac*alpha_molar_mass)
        
        return(round(alpha_phafrac, 3), round(beta_phafrac, 3))
        
    # set sequential refinement parameters 
    def gpx_refine(num1, num2, cycle_num=5, seq_copy_flag=True, rev_seq_flag=False):
        gpx.set_Controls('sequential', gpx_hist_from_num(num1, num2))
        gpx.set_Controls('cycles', cycle_num)
        gpx.set_Controls('seqCopy', seq_copy_flag)
        gpx.set_Controls('Reverse Seq', rev_seq_flag)
        gpx.refine()
    
    # define the HAP entry initialization at beta_str_num, alpha_disapp_num, end_num
    beta_str_phafrac = phafrac_from_wgtfrac(steel_wgtfrac=0.15, alpha_wgtfrac=0.7)
    beta_str_ent = [{'HStrain': [[-0.0025, 0.0], [True, False]],
                     'Scale': [1, False],
                     'Mustrain': ['isotropic', [15000, 1000.0, 0.0], 
                                  [True, False, False],[0, 0, 1], [0.01, 0.01],[False, False]],
                     'Pref.Ori.':['SH', 1.0, True,[0, 0, 1], 8, 
                                  {'C(4,1)': -0.042, 'C(6,1)': 0.085,'C(8,1)': 0.14}, [''], 0.1],} ,
                    {'HStrain': [[-0.0025, -0.001], [True, True]],
                     'Scale': [beta_str_phafrac[0], True],
                     'Mustrain': ['isotropic', [15000, 1000.0, 0.0], 
                                  [True, False, False],[0, 0, 1], [0.01, 0.01, 0.01], [False, False, False]],
                     'Pref.Ori.':['SH', 1.0, True, [0, 0, 1], 4,
                                  {'C(2,0)': 0.07, 'C(4,0)': 0.06},[''],0.1],} ,
                    {'HStrain': [[-0.004, 0.0], [True, False]],
                     'Scale': [beta_str_phafrac[1], True], 
                     'Mustrain': ['isotropic', [15000.0, 1000.0, 0.0], 
                                  [True, False, False],[0, 0, 1], [0.01, 0.01], [False, False]],
                     'Pref.Ori.': ['SH', 1.0, True, [0, 0, 1], 4, 
                                   {'C(4,1)': -0.427}, [''], 0.1]} ,]

    
    no_alpha_phafrac = phafrac_from_wgtfrac(steel_wgtfrac=0.2, alpha_wgtfrac=0.0)
    no_alpha_ent = [{'HStrain': [[-0.003, 0.0], [True, False]],
                     'Scale': [1, False],
                     'Mustrain': ['isotropic', [15000, 1000.0, 0.0], 
                                  [True, False, False],[0, 0, 1], [0.01, 0.01],[False, False]],
                     'Pref.Ori.':['SH', 1.0, True,[0, 0, 1], 8, 
                                  {'C(4,1)': -0.444, 'C(6,1)': 0.077,'C(8,1)': 0.257}, [''], 0.1],} ,
                    {'HStrain': [[0.0, 0.0], [False, False]],
                     'Scale': [0, False],
                     'Mustrain': ['isotropic', [15000, 1000.0, 0.0], 
                                  [False, False, False],[0, 0, 1], [0.01, 0.01, 0.01], [False, False, False]],
                     'Pref.Ori.':['SH', 1.0, False, [0, 0, 1], 4,
                                  {'C(2,0)': 0.0, 'C(4,0)': 0.0},[''],0.1],} ,
                    {'HStrain': [[-0.006, 0.0], [True, False]],
                     'Scale': [no_alpha_phafrac[1], True], 
                     'Mustrain': ['isotropic', [15000.0, 1000.0, 0.0], 
                                  [True, False, False],[0, 0, 1], [0.01, 0.01], [False, False]],
                     'Pref.Ori.': ['SH', 1.0, True, [0, 0, 1], 4, 
                                   {'C(4,1)': -1.484}, [''], 0.1]} ,]
    
    
    ending_phafrac = phafrac_from_wgtfrac(steel_wgtfrac=0.2, alpha_wgtfrac=0.7)
    ending_ent =   [{'HStrain': [[-0.0001, 0.0], [True, False]],
                     'Scale': [1, False],
                     'Mustrain': ['isotropic', [15000, 1000.0, 0.0], 
                                  [True, False, False],[0, 0, 1], [0.01, 0.01],[False, False]],
                     'Pref.Ori.':['SH', 1.0, True,[0, 0, 1], 8, 
                                  {'C(4,1)': 0.256, 'C(6,1)': -0.172,'C(8,1)': -0.581}, [''], 0.1],} ,
                    {'HStrain': [[-7e-5, -0.0001], [True, True]],
                     'Scale': [ending_phafrac[0], True],
                     'Mustrain': ['isotropic', [15000, 1000.0, 0.0], 
                                  [False, False, False],[0, 0, 1], [0.01, 0.01, 0.01], [False, False, False]],
                     'Pref.Ori.':['SH', 1.0, False, [0, 0, 1], 4,
                                  {'C(2,0)': 0.189, 'C(4,0)': -0.1},[''],0.1],} ,
                    {'HStrain': [[-0.002, 0.0], [True, False]],
                     'Scale': [ending_phafrac[1], True], 
                     'Mustrain': ['isotropic', [15000.0, 1000.0, 0.0], 
                                  [True, False, False],[0, 0, 1], [0.01, 0.01], [False, False]],
                     'Pref.Ori.': ['SH', 1.0, True, [0, 0, 1], 4, 
                                   {'C(4,1)': -1.484}, [''], 0.1]} ,]
            
            
    # turn off or off unit cell, pref. ori., and mustrain refinement initially
    # turn on histrain Dijs, and not eA for cubic structures
        
    for i, p in enumerate(gpx.phases()):
        p.set_refinements({'Cell': False})
        p.set_HAP_refinements({'Pref.Ori.': True})
        p.set_HAP_refinements({'Mustrain': {"type": 'isotropic', 'refine': True}}) 
        if i == 0:
            # for steel turn on D11 refinement and turn off eA refinment
            for j, hstrainEnt in enumerate(gpx.phase(i).getHAPentryList(0, 'HStrain')):
                p.setHAPentryValue(hstrainEnt[0], [[0.0, 0.0], [True, False]])
            
        elif i==1:
            # for alpha turn on D11 and D33 refinement
            p.set_HAP_refinements({"HStrain": True})
        else:
            # for beta set an initial Mustrain value
            # should already been done in the previous single import script
            p.set_HAP_refinements({"Mustrain": {'value': 15000}})

    # gpx.get_Constraints('HAP')[:] = []
    
    # load in rest of the powder data
    for data in diff["data"]:
        hist = gpx.add_powder_histogram(DataWrap(data), DataWrap(diff["prm"]),
                                        fmthint=diff["fmthint"],
                                        databank=diff["databank"],
                                        instbank=diff["databank"],
                                        phases='all')
    
        

    # copy HAP values, background, instrument params. limits, and sample params.
    gpx.copyHistParms(0, 'all', ['b', 'i', 'l', 's'])
    for p in gpx.phases():
        p.copyHAPvalues(0, 'all')
    
    
    # gpx.add_EqnConstr(1.0,
    #                  [f'{i}:*:Scale' for i, _ in enumerate(gpx.phases())],
    #                  [1]*len(gpx.phases()))
    
    
    
    # when beta appears and till the end: 
    # turn on phase fraction of beta, and set it to around beta weight fraction of 0.1
    # turn on HStrain D11 flag of beta and set D11 to -0.003 
    beta_app_phafrac = phafrac_from_wgtfrac(steel_wgtfrac=0.2, alpha_wgtfrac=0.7)
    for j, (sclEnt, hstrainEnt) in enumerate(zip(gpx.phase(1).getHAPentryList(0, 'Scale'),
                                                 gpx.phase(1).getHAPentryList(0, 'HStrain'))):
        hist_name = sclEnt[0][0]
        scan_num = num_from_hist_name(hist_name)
        if scan_num >= (beta_app_num-1) and scan_num <= end_num:
            gpx.phase(2).setHAPentryValue(sclEnt[0], [beta_app_phafrac[1], True])
            gpx.phase(2).setHAPentryValue(hstrainEnt[0], [[-0.003, 0.0], [True, False]])
    
    
    # when alpha disappears and before alpha disappears: 
    # turn off Scale of alpha and set to 0, 
    # set the Mustrain to be around 15000
    # set and the HStrain and Pref. Ori. to 0 
    # no_alpha_hist_names = []
    for j, (sclEnt, hstrainEnt, preforiEnt, mustrainEnt) in enumerate(
            zip(gpx.phase(1).getHAPentryList(0, 'Scale'),
                gpx.phase(1).getHAPentryList(0, 'HStrain'),
                gpx.phase(1).getHAPentryList(0, 'Pref.Ori.'),
                gpx.phase(1).getHAPentryList(0, 'Mustrain'))):
        hist_name = sclEnt[0][0]
        scan_num = num_from_hist_name(hist_name)
        if scan_num >= (alpha_disapp_num-1) and scan_num < (alpha_reapp_num):
            gpx.phase(1).setHAPentryValue(sclEnt[0], [0.0, False])
            gpx.phase(1).setHAPentryValue(hstrainEnt[0], [[0.0, 0.0], [False, False]])
            gpx.phase(1).setHAPentryValue(preforiEnt[0], ['SH', 1.0, False, [0, 0, 1], 4, 
                                                    {'C(2,0)': 0.0, 'C(4,0)': 0.0}, [''], 0.1])
            gpx.phase(1).setHAPentryValue(mustrainEnt[0], ['isotropic',[15000, 1000.0, 0.0],
                                                    [False, False, False], [0, 0, 1], 
                                                    [0.01, 0.01, 0.01],[False, False, False]])
            #no_alpha_hist_names.append(hist_name)
    # supposedly another way to turn off the HAP entry by specifying the histograms explicitly
    #gpx.phase(1).set_HAP_refinements({'Mustrain': {"value": 18500}}, histograms=gpx_hist(no_alpha_hist_names)) 
    
    
    # sequential refinement from start to before beta reappears: forward, copy results to next
    gpx.set_Controls('sequential', gpx_hist_from_num(start_num, beta_app_num))
    gpx.set_Controls('cycles', 5)
    gpx.set_Controls('seqCopy', True)
    gpx.set_Controls('Reverse Seq', False)
    gpx.refine()
    
    # when beta is strong initializing entry values from a previous rietveld result
    HAP_entry_for_num(beta_str_ent, beta_str_num, beta_str_num+1)
    
    # sequential refinement from beta strong to beta appears: reverse, copy results to next
    gpx_refine(beta_app_num, beta_str_num+1, cycle_num=5, seq_copy_flag=True, rev_seq_flag=True)
    
    # sequential refinement from beta strong to alpha disappears: forward, copy results to next
    gpx_refine(beta_str_num, alpha_disapp_num, cycle_num=5, seq_copy_flag=True, rev_seq_flag=False)
    
    # sequential refinement from alpha disappear to before reappears: forward, copy results to next
    HAP_entry_for_num(no_alpha_ent, alpha_disapp_num, alpha_disapp_num+1)
    gpx_refine(alpha_disapp_num, alpha_reapp_num, cycle_num=10, seq_copy_flag=True, rev_seq_flag=False)
    
    # sequential refinement from end to before alpha reappears: reverse, copy results to next
    HAP_entry_for_num(ending_ent, end_num, end_num+1)
    gpx_refine(alpha_reapp_num, end_num, cycle_num=5, seq_copy_flag=True, rev_seq_flag=True)
    
    # sequential refinement from alpha reappear to end of isothermal hold: forward, copy results to next
    gpx_refine(alpha_reapp_num, iso_hold_end_num, cycle_num=5, seq_copy_flag=True, rev_seq_flag=False)
    
    # sequential refinment without copying to next
    gpx_refine(start_num, end_num, cycle_num=5, seq_copy_flag=False, rev_seq_flag=False)
    
    # sequential refinement from end to 
    gpx.save()
    



if __name__ == '__main__':
    main()
