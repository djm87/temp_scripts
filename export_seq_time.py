#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export Rietveld refinement results 
"""

import os
import sys
# sys.path.insert(0, '/home/bonnie/anaconda3/envs/GSASII/gsasii/GSASII')
sys.path.insert(0, 'C://Users/sanoj737/anaconda3/envs/G2/GSASII')
import GSASIIscriptable as G2sc
import math
import csv
import numpy as np
from scipy.optimize import root_scalar
from datetime import datetime
from typing import List
import pandas as pd
import glob

folder = 'V5-7_insitu'
gpx_name = folder+'_seq_int1_sum_norm_dc_noDispXY.gpx'
pwd_data_type = 'tif' #'fxye'
lamda_val = 0.1839 # angstrom
T_ref_val = 25 # deg C
exposure_time = 0.16 # s
ref_idx_val = 0
beta_CTE_val = 8.8 # *10^-6/K
beta_CTE_err_val = 0.5 # *10^-6/K
bankid = 1

def PathWrap(fil) -> str:
    """Make path relative to current directory."""
    return os.path.join("./", fil)


def DataWrap(fil) -> str:
    """Make path relative to current directory."""
    return os.path.join("./", fil)


def error_divmult(Q: float, par: List[float], pare: List[float]) -> float:
    """Propagates the error for addition."""
    # see https://www.statisticshowto.com/error-propagation/
    sum_div = 0
    for p, pe in zip(par, pare):
        if p != 0:
            sum_div += (pe/p)**2
    return abs(Q)*math.sqrt(sum_div)


def error_sine(theta, thetae) -> float:
    """Propagates the error for sine."""
    # see https://math.stackexchange.com/questions/2653896/error-propagation-for-sine-function-in-terms-of-degrees-and-radians
    return math.cos(theta)*thetae


def error_cos(theta, thetae) -> float:
    """Propagates the error for cos."""
    # see https: // math.stackexchange.com/questions/2653896/error-propagation-for-sine-function-in-terms-of-degrees-and-radians
    return math.sin(theta)*thetae


def error_add(pare: List[float]) -> float:
    """Propagates the error for addition."""
    # see https: // www.statisticshowto.com/error-propagation /
    sum_sq = 0
    for pe in pare:
        sum_sq += (pe)**2
    return math.sqrt(sum_sq)


def error_powers(Q, x, xe, n) -> float:
    """Propagates the error for powers with exact numbers."""
    # see https: // math.stackexchange.com/questions/2653896/error-propagation-for-sine-function-in-terms-of-degrees-and-radians
    if x == 0.0:
        return 0.0
    return abs(Q)*abs(n)*xe/abs(x)
    


def error_roots(Q, x, xe, n) -> float:
    """Propagates the error for roots."""
    # see https: // math.stackexchange.com/questions/2653896/error-propagation-for-sine-function-in-terms-of-degrees-and-radians
    
    return abs(Q)*xe/(abs(x)*abs(n))
 
   
def error_exact(A, xe):
    """Propagates the error for exact Measurement."""
    # see https: // math.stackexchange.com/questions/2653896/error-propagation-for-sine-function-in-terms-of-degrees-and-radians
    return abs(A)*xe


def d_from_2theta(twoth: float, twoth_err: float, lamda: float = lamda_val ) -> List[float]:
    """Converts two theta to dspacing"""
    par_sin = np.sin(np.deg2rad(twoth/2))
    par_sin_err = error_sine(np.deg2rad(twoth/2), np.deg2rad(twoth_err/2))
    d = lamda / 2 / par_sin
    d_err = lamda/2 * error_divmult(1/par_sin, [1, par_sin], [0, par_sin_err])
    
    return d, d_err    

def strain_cal(par: float, par_ref: float, par_esd: float, par_ref_esd: float, scale: float = 1e6) -> List[float]:
    """Converts lattice parameter to strain (microstrain default)."""
    strain = (par / par_ref - 1)*scale
    strain_err = scale*error_divmult(par/par_ref, [par, par_ref], [par_esd, par_ref_esd])
    return strain, strain_err



def SS316_temperature_cal(strain: float, strain_err: float, T_ref: float = T_ref_val) -> List[float]:
    """
    Convert strain to temperature for SS316

    Parameters
    ----------
    strain : float
        Strain in microstrain.
    strain_err: float
        Strain error in microstrain.
    T_ref : float
        Temperature in Centigrade.

    Returns
    -------
    float
        Temperature and error in Centigrade.

    """
    
    temperature = T_ref + 63032*strain/1e6 - 1003300*(strain/1e6)**2 + 20171000*(strain/1e6)**3
    temperature_err1 = 63032*strain_err/1e6   
    temperature_err2 = 1003300*error_powers((strain/1e6)**2, strain/1e6, strain_err/1e6, 2)
    temperature_err3 = 20171000*error_powers((strain/1e6)**3, strain/1e6, strain_err/1e6, 3)
    temperature_err = error_add([temperature_err1, temperature_err2, temperature_err3])
    
    return temperature, temperature_err


def temp_cal_from_pos(pos:List[float], pos_err:List[float], ref_idx: int =ref_idx_val, 
                      T_ref:float = T_ref_val, lamda:float= lamda_val)-> List[float]:
    '''
    Convert peak position from SPF to temperature for SS316
    Note that the peak position and error are lists
    '''
    par, par_err = np.vectorize(d_from_2theta)(pos, pos_err)
    strain, strain_err = np.vectorize(strain_cal)(par, par[ref_idx_val], par_err, par_err[ref_idx_val])
    temp, temp_err = np.zeros(len(pos)), np.zeros(len(pos))
    
    # the strain of at the lattice parameter reference will be 0
    # which fails during error calculation for temperature 
    temp[ref_idx_val] =  T_ref + 63032*strain[ref_idx_val]/1e6 - 1003300*(strain[ref_idx_val]/1e6)**2 + 20171000*(strain[ref_idx_val]/1e6)**3
    ref_mask = np.ones(len(temp), bool)
    ref_mask[ref_idx_val] = False
    temp[ref_mask], temp_err[ref_mask] = np.vectorize(SS316_temperature_cal)(strain[ref_mask], strain_err[ref_mask], T_ref)
    return temp, temp_err


def v_conc_cal(par: float, par_ref: float, par_esd: float, par_ref_esd: float, temp: float, temp_err: float, scale: float = 1e6,
                T_ref: float = T_ref_val, beta_CTE: float = beta_CTE_val, beta_CTE_err: float = beta_CTE_err_val)-> List[float]:
    '''
    Calculate vanadium concentration in beta phase from lattice paramater
    '''
    # total strain
    beta_str_total, beta_str_total_err = strain_cal(par, par_ref, par_esd, par_ref_esd)
    # non-thermal expansion strain
    beta_str_nonTE = beta_str_total-beta_CTE*(temp-T_ref)
    beta_str_nonTE_err1 =  error_divmult(Q=beta_CTE*(temp-T_ref), par=[beta_CTE, (temp-T_ref)], 
                                         pare=[beta_CTE_err, temp_err])
    beta_str_nonTE_err2 = error_add(pare=[beta_str_total_err, beta_str_nonTE_err1])
    # V concentration in atomic and weight fraction
    beta_V_at = 0.20771-13.463*beta_str_nonTE/scale-0.004
    beta_V_at_err = 13.463*beta_str_nonTE_err2/scale
    beta_V_wt = 0.20771-13.375*beta_str_nonTE/scale
    beta_V_wt_err = 13.375*beta_str_nonTE_err2/scale
    
    return(beta_V_at, beta_V_at_err, beta_V_wt, beta_V_wt_err)

def phafrac_alpha_cal(scale_alpha: float, scale_beta: float, scale_alpha_esd: float, scale_beta_esd: float)-> List[float]:
    '''
    Calculate alpha phase fraction
    '''
    phafrac_alpha = scale_alpha/(scale_alpha+scale_beta)
    if scale_alpha_esd is None:
        scale_alpha_esd = 0
    if scale_beta_esd is None:
        scale_beta_esd = 0
    phafrac_alpha_esd1 = error_add([scale_alpha_esd, scale_beta_esd])
    phafrac_alpha_esd2 = error_divmult(scale_alpha/(scale_alpha+scale_beta), 
                                       [scale_alpha, (scale_alpha+scale_beta)], 
                                       [scale_alpha_esd, phafrac_alpha_esd1])
    
    return (phafrac_alpha, phafrac_alpha_esd2)

def read_reflist(reflist: list) -> dict:
    """Move the reflist obj to dictionary of peaks."""
    dic = {}
    for ref in reflist:
        dic2 = {}
        dic2["multiplicity"] = ref[3]
        dic2["d"] = ref[4]
        dic2["pos"] = ref[5]
        dic2["sig"] = ref[6]
        dic2["gam"] = ref[7]
        dic2["Fosq"] = ref[8]
        dic2["Fcsq"] = ref[9]
        dic2["phase"] = ref[10]
        dic2["Icorr"] = ref[11]
        dic2["Prfo"] = ref[12]
        dic2["Trans"] = ref[13]
        dic2["ExtP"] = ref[14]
        # print(ref)
        dic[f'{[int(ref[i]) for i in range(0,3)]}'] = dic2
    return dic


def rietveld_output_dict():
    """Perform single peak refinements for each pkfits dictionary and histogram."""
    gpx = G2sc.G2Project(gpxfile=PathWrap(gpx_name))  # load a project
    seq = gpx.seqref()
    keys = ['GE_filenumber', 'Hist_name', 'Frame_number', 'Azm', 'Temp_ss_riet', 'Temp_ss_riet_esd',
            'V_conc_at', 'V_conc_at_esd', 'V_conc_wt', 'V_conc_wt_esd', 'Phase_frac_alpha', 'Phase_frac_alpha_esd',
            'Rwp', 'chisq', 'DelChi2', 'Scale', 'Scale_esd', 'WgtFrac', 'WgtFrac_esd', 'Mustrain',
            'Mustrain_esd', 'Mustrainmx', 'Mustrainmx_esd', 'Size', 'Size_esd', 'Sizemx', 'Sizemx_esd', 'a',
            'a_esd', 'a_str', 'a_str_esd', 'b', 'b_esd', 'b_str', 'b_str_esd', 'c', 'c_esd', 'c_str', 'c_str_esd',
            'alpha', 'alpha_esd', 'alpha_str', 'alpha_str_esd', 'beta', 'beta_esd', 'beta_str', 'beta_str_esd',
            'gamma', 'gamma_esd', 'gamma_str', 'gamma_str_esd', 'vol', 'vol_esd', 'vol_str', 'vol_str_esd']
    # 
    d = {}
    for key in keys:
        if key not in ['GE_filenumber', 'Hist_name', 'Frame_number', 'Azm', 'Temp_ss_riet', 'Temp_ss_riet_esd',
                       'V_conc_at', 'V_conc_at_esd', 'V_conc_wt', 'V_conc_wt_esd', 'Phase_frac_alpha', 'Phase_frac_alpha_esd', 
                       'Rwp', 'chisq', 'DelChi2']:
            for j, _ in enumerate(gpx.phases()):
                d[f'{key}_{j}'] = [None]*len(seq.histograms())
        else:
            d[key] = [None]*len(seq.histograms())

    for i, hist in enumerate(seq.histograms()):
        data = seq[hist]
        
        # histogram names
        d['Hist_name'][i] = data['title']
        
        # scan number and Azm depending on powder data type
        if pwd_data_type == 'tif':
            
            d['GE_filenumber'][i] = int(d['Hist_name'][i].split('.tif')[0].split('_')[1])
            d['Frame_number'][i] = int(d['Hist_name'][i].split('.tif')[0].split('_')[-1])
            d['Azm'][i] = int(d['Hist_name'][i].split('.tif')[1].split('Azm=')[1].split('.')[0])
        else:
            d['GE_filenumber'][i] = int(d['Hist_name'][i].split('.fxye')[0].split('_')[1])  
            if 'ge2' in d['Hist_name'][i]:
                d['Frame_number'][i] = int(d['Hist_name'][i].split('.fxye')[0].split('_#')[1].split('_')[0])
            else:
                d['Frame_number'][i] = 0
            d['Azm'][i] = int(d['Hist_name'][i].split('.fxye')[0].split('_A')[1].split('_')[0])
        
        # Rvals
        d['chisq'][i] = data['Rvals']['chisq']
        d['Rwp'][i] = data['Rvals']['Rwp']
        d['DelChi2'][i] = data['Rvals']['DelChi2']

        # From hist get the Rietveld reflection information
        reflist = []
        for phase, val in gpx.histograms(hist)[0].reflections().items():
            reflist.append(read_reflist(val['RefList']))

        # Add reflections if they don't exist and initialize all others to None
        for j, ref in enumerate(reflist):
            for hkl, val in ref.items():
                for key in ['d', 'pos', 'sig', 'gam', 'Icorr', 'Prfo']:
                    if f'{key}_{j}:{hkl}' not in d:
                        d[f'{key}_{j}:{hkl}'] = [None]*len(seq.histograms())
                    d[f'{key}_{j}:{hkl}'][i] = val[key]

        for j, phase in enumerate(gpx.phases()):

            # Cell parameters and esd
            cell = seq.get_cell_and_esd(j, i)
            d[f'a_{j}'][i] = cell[0][0]
            d[f'b_{j}'][i] = cell[0][1]
            d[f'c_{j}'][i] = cell[0][2]
            d[f'alpha_{j}'][i] = cell[0][3]
            d[f'beta_{j}'][i] = cell[0][4]
            d[f'gamma_{j}'][i] = cell[0][5]
            d[f'vol_{j}'][i] = cell[0][6]
            d[f'a_esd_{j}'][i] = cell[1][0]
            d[f'b_esd_{j}'][i] = cell[1][1]
            d[f'c_esd_{j}'][i] = cell[1][2]
            d[f'alpha_esd_{j}'][i] = cell[1][3]
            d[f'beta_esd_{j}'][i] = cell[1][4]
            d[f'gamma_esd_{j}'][i] = cell[1][5]
            d[f'vol_esd_{j}'][i] = cell[1][6]

            # Cell parameterstrains and esd
            d[f'a_str_{j}'][i], d[f'a_str_esd_{j}'][i] = strain_cal(
                d[f'a_{j}'][i], d[f'a_{j}'][0], d[f'a_esd_{j}'][i], d[f'a_esd_{j}'][0])
            d[f'b_str_{j}'][i], d[f'b_str_esd_{j}'][i] = strain_cal(
                d[f'b_{j}'][i], d[f'b_{j}'][0], d[f'b_esd_{j}'][i], d[f'b_esd_{j}'][0])
            d[f'c_str_{j}'][i], d[f'c_str_esd_{j}'][i] = strain_cal(
                d[f'c_{j}'][i], d[f'c_{j}'][0], d[f'c_esd_{j}'][i], d[f'c_esd_{j}'][0])
            d[f'alpha_str_{j}'][i], d[f'alpha_str_esd_{j}'][i] = strain_cal(
                d[f'alpha_{j}'][i], d[f'alpha_{j}'][0], d[f'alpha_esd_{j}'][i], d[f'alpha_esd_{j}'][0])
            d[f'beta_str_{j}'][i], d[f'beta_str_esd_{j}'][i] = strain_cal(
                d[f'beta_{j}'][i], d[f'beta_{j}'][0], d[f'beta_esd_{j}'][i], d[f'beta_esd_{j}'][0])
            d[f'gamma_str_{j}'][i], d[f'gamma_str_esd_{j}'][i] = strain_cal(
                d[f'gamma_{j}'][i], d[f'gamma_{j}'][0], d[f'gamma_esd_{j}'][i], d[f'gamma_esd_{j}'][0])
            d[f'vol_str_{j}'][i], d[f'vol_str_esd_{j}'][i] = strain_cal(
                d[f'vol_{j}'][i], d[f'vol_{j}'][0], d[f'vol_esd_{j}'][i], d[f'vol_esd_{j}'][0])

            # Extract the phase parameters and esd
            if seq.get_Variable(i, f'{j}:{i}:Scale') is not None:
                d[f'Scale_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Scale')[0]
                d[f'WgtFrac_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:WgtFrac')[0]
                d[f'Mustrain_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Mustrain;i')[0]
                d[f'Mustrainmx_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Mustrain;mx')[0]
                d[f'Size_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Size;i')[0]
                d[f'Sizemx_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Size;mx')[0]
                d[f'Scale_esd_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Scale')[1]
                d[f'WgtFrac_esd_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:WgtFrac')[1]
                d[f'Mustrain_esd_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Mustrain;i')[1]
                d[f'Mustrainmx_esd_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Mustrain;mx')[1]
                d[f'Size_esd_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Size;i')[1]
                d[f'Sizemx_esd_{j}'][i] = seq.get_Variable(i, f'{j}:{i}:Size;mx')[1]
        
        # Calculate temperature from steel lattice parameter
        d['Temp_ss_riet'][i], d['Temp_ss_riet_esd'][i]  = SS316_temperature_cal(d['a_str_0'][i], d['a_str_esd_0'][i], T_ref=T_ref_val)
        # Calculate alpha phase fraction
        d['Phase_frac_alpha'][i], d['Phase_frac_alpha_esd'][i] = phafrac_alpha_cal(d['Scale_1'][i], d['Scale_2'][i], d['Scale_esd_1'][i], d['Scale_esd_2'][i])
        
        
    # Calculate vanadium concentration    
    for i, hist in enumerate(seq.histograms()):
        if d['a_str_2'][i] > 0:
            d['V_conc_at'][i], d['V_conc_at_esd'][i], d['V_conc_wt'][i], d['V_conc_wt_esd'][i] \
                = v_conc_cal(par = d['a_2'][i], par_ref = d['a_2'][-1], par_esd = d['a_esd_2'][i], par_ref_esd = d['a_esd_2'][-1],
                             temp = d['Temp_ss_riet'][i], temp_err = d['Temp_ss_riet_esd'][i])
        else:
            d['V_conc_at'][i], d['V_conc_at_esd'][i], d['V_conc_wt'][i], d['V_conc_wt_esd'][i] = 0, 0, 0, 0
        
        
    return d


if __name__ == '__main__':
    d = rietveld_output_dict()
    df = pd.DataFrame(d)
    
    # import in time csv
    time_df = pd.read_csv(PathWrap(folder+'_time.csv'), index_col=[0])
    
    # merge two dataframe
    riet_df = time_df[['Time_creat', 'Time_mod', 'Time_diff', 'Temperature_pyro', 'GE_filenumber']].merge(df,
                       how='inner', on='GE_filenumber')
    
    
    
    # fill in missing value 
    # pyrometer temperature: forward or backward fill
    riet_df['Temperature_pyro'] = riet_df['Temperature_pyro'].fillna(method='bfill').fillna(method='ffill')
    # time difference: average
    riet_df['Time_diff'] = riet_df['Time_diff'].fillna(round(riet_df['Time_diff'].mean(),3))
    # creation time: modification time - average time difference
    riet_df['Time_creat'] = riet_df['Time_creat'].fillna(riet_df['Time_mod']-riet_df['Time_diff'])
    riet_df = riet_df.drop(columns=['Time_mod', 'Time_diff'])
    # riet steel temperature esd: 0
    riet_df['Temp_ss_riet_esd'] = riet_df['Temp_ss_riet_esd'].fillna(0)
    
    # calculate time of each frame based on exposure time
    riet_df['Time_creat'] = riet_df['Time_creat']+riet_df['Frame_number']*exposure_time
    
    # calculate elapsed time
    riet_df.insert(loc=riet_df.columns.get_loc("Time_creat")+1, column='Time_elap_creat', 
                   value=riet_df['Time_creat']-riet_df['Time_creat'][0])
    
    # import in spf csv
    if os.path.exists(PathWrap((gpx_name.split('.gpx')[0])+'_spf.csv')):
        spf_df = pd.read_csv(PathWrap((gpx_name.split('.gpx')[0])+'_spf.csv'), index_col=[0])
        spf_df = spf_df.dropna(how='all')
        
        # calculate temperature based on spf of steel peaks from spf
        steel_peaks = ['111', '200', '220', '311', '222']
        
        for num, peak in enumerate(steel_peaks):
            temp_ss, temp_ss_err = temp_cal_from_pos(spf_df['pos0_'+peak[0]+'_'+peak[1]+'_'+peak[2]], 
                                                     spf_df['pos_esd0_'+peak[0]+'_'+peak[1]+'_'+peak[2]], T_ref=T_ref_val)
            riet_df.insert(loc=riet_df.columns.get_loc("Temp_ss_riet_esd")+num*2+1, column='Temp_ss_spf_'+peak[0]+peak[1]+peak[2], value=temp_ss)
            riet_df.insert(loc=riet_df.columns.get_loc("Temp_ss_riet_esd")+num*2+2, column='Temp_ss_spf_esd_'+peak[0]+peak[1]+peak[2], value=temp_ss_err)

        
    # write to csv
    riet_df.to_csv(PathWrap((gpx_name.split('.gpx')[0])+'_riet.csv'))
    
