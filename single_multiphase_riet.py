#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perform single Rietveld refinement for Ti64 + Steel
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
base_fname = '_int1_sum_norm_dc'
gpx_name_out = folder_name+'_single'+base_fname+'.gpx'

# import in desired histograms
search_string = folder_name+base_fname+'/*.fxye'
hists = []

bad_scans = []
good_scans = [4329, 4769]

for hist in glob.glob(search_string):
    hist_name = os.path.basename(hist)
    scan_num = int(hist_name.split('_')[1].split('.')[0])
    if scan_num not in bad_scans and scan_num >= good_scans[0] and scan_num <= good_scans[1]:
        hists.append(hist)

hists.sort(key=natural_keys)

# phases: 0 for steel, 1 for alpha, and 2 for beta
phases = [{'cif_name': "../../1_Calibration/steel.cif", 'phasename': 'steel', 'fmthint': 'CIF', 'initial_frac': 1},
          {'cif_name': "../../1_Calibration/alphap.cif", 'phasename': 'alphap', 'fmthint': 'CIF', 'initial_frac': 7},
          {'cif_name': "../../1_Calibration/beta.cif", 'phasename': 'beta', 'fmthint': 'CIF', 'initial_frac': 0.0}]

# background
back = {'type': 'chebyschev-1', 'nocoeffs': 7}

# limit
limits = [2.8, 11.2]

# sample parameter
sample_param = {'gonio_rad': 976.1, 'disp_x': 181.484, 'disp_y': 4027.77}

# size and microstrain
samps = [{"Size": {'type': 'isotropic', 'refine': False, 'value': 1, 'LGmix': 0.0},
         "Mustrain": {'type': 'isotropic', 'refine': False, 'value': 15000, 'LGmix': 0.0}, },
         {"Size": {'type': 'isotropic', 'refine': False, 'value': 1, 'LGmix': 0.0},
         "Mustrain": {'type': 'isotropic', 'refine': False, 'value': 15000, 'LGmix': 0.0}, },
         {"Size": {'type': 'isotropic', 'refine': False, 'value': 1, 'LGmix': 0.0},
         "Mustrain": {'type': 'isotropic', 'refine': False, 'value': 15000, 'LGmix': 0.0}, }
         ]

diff = {"data": [hists[0]],
        "prm": "../../1_Calibration/Instrument_july11.instprm",
        "fmthint": "GSAS powder",
        "databank": 1}

preferred_orientation = [['SH', 1.0, False, [0, 0, 1], 8, {
    'C(4,1)': 0.0, 'C(6,1)': 0.0, 'C(8,1)': 0.0}, [''], 0.1],
    ['SH', 1.0, False, [0, 0, 1], 4, {
        'C(2,0)': 0.0, 'C(4,0)': 0.0}, [''], 0.1],
    ['SH', 1.0, False, [0, 0, 1], 4, {'C(4,1)': 0.0}, [''], 0.1]]


def main():

    gpx = G2sc.G2Project(filename=gpx_name_out)
    for data in diff["data"]:
        hist = gpx.add_powder_histogram(DataWrap(data), DataWrap(diff["prm"]),
                                        fmthint=diff["fmthint"], databank=diff["databank"], instbank=diff["databank"])
        # setting the sample parameter for each histogram
        hist['Sample Parameters']['Gonio. radius'] = sample_param['gonio_rad']
        hist['Sample Parameters']['DisplaceX'] = [sample_param['disp_x'],False]
        hist['Sample Parameters']['DisplaceY'] = [sample_param['disp_y'],False]
        
    # add a phases to the project
    for i, (phase, samp) in enumerate(zip(phases, samps)):
        gpx.add_phase(PathWrap(phase['cif_name']),
                      phasename=phase['phasename'], fmthint=phase['fmthint'], histograms=gpx.histograms())

        # Loop over the histograms and set the phase fractions
        for j, sclEnt in enumerate(gpx.phase(i).getHAPentryList(0, 'Scale')):
            gpx.phase(i).setHAPentryValue(sclEnt[0], [phase['initial_frac'], False])
            gpx.phase(i).setSampleProfile(j, 'size',
                                          samp['Size']['type'],
                                          samp['Size']['value'],
                                          LGmix=samp['Size']['LGmix'])
            gpx.phase(i).setSampleProfile(j, 'microstrain',
                                          samp['Mustrain']['type'],
                                          samp['Mustrain']['value'],
                                          LGmix=samp['Mustrain']['LGmix'])
        
        # Loop over histograms and set the preferred orientation
        for j, sclEnt in enumerate(gpx.phase(i).getHAPentryList(0, 'Pref.Ori.')):
            gpx.phase(i).setHAPentryValue(sclEnt[0], preferred_orientation[i])

    # constraint phases to sum to 1
    # gpx.add_EqnConstr(1.0,
    #                   [f'{i}:0:Scale' for i, _ in enumerate(phases)],
    #                   [1]*len(phases))
    
    # set controls  
    gpx['Controls']['data']['max cyc'] = 5
    # gpx.set_Controls('cycles', 5)
    
    # refine background and histogram scale
    pardict = [{'set': {'Limits': limits,
                        'Background': {'type': back['type'], "no. coeffs": back['nocoeffs'], 'refine': True}
                        }
                }]
               
    # gpx.set_refinement(pardict[0])
    gpx.do_refinements(pardict)   # calculate pattern
    
    # refine previous + phase fraction for alpha only 
    pardict = {'Scale': True}
    gpx.phase(1).set_HAP_refinements(pardict)
    gpx.do_refinements()
    gpx.save()

    # refine pevious + lattice parameters for steel and alpha
    pardict = {'Cell': True}
    gpx.phase(0).set_refinements(pardict)
    gpx.phase(1).set_refinements(pardict)
    gpx.do_refinements()
    gpx.save()

    
    # #p.getHAPentryValue(mstrEnt[0][0])[2][0] = False
    # #gpx.add_HoldConstr(['2::A0', '2::A2'])
    # gpx.do_refinements()
    
    # refine previous + texture for steel and alpha
    pardict = {'Pref.Ori.': True}
    gpx.phase(0).set_HAP_refinements(pardict)
    gpx.phase(1).set_HAP_refinements(pardict)
    gpx.do_refinements()
    gpx.save()
    
    
    # refine previous + strain for steel and alpha
    pardict = {'Mustrain': {"type": 'isotropic', 'refine': True}}
    gpx.phase(0).set_HAP_refinements(pardict)
    gpx.phase(1).set_HAP_refinements(pardict)
    gpx.do_refinements()
    gpx.save()



if __name__ == '__main__':
    main()
