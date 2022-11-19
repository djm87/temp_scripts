import shutil
import time
from pathlib import Path
from typing import List
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cv2
import fabio
import numpy as np
import pyFAI
from functools import partial
from multiprocessing import Pool, freeze_support
import pandas as pd

NCPUS = 30

watch_dir = Path.cwd() / "tempdata/weld_17_156_RT_P15_center_BS250"
watch_dir.mkdir(exist_ok=True)

poni_file = "cal_nov17_451pm.poni"
mask_file = "mask.edf"
ai = pyFAI.load(poni_file)
mask = fabio.open(mask_file)

image_log_dir = 'image_log'
image_dir = 'image'
spectra_dir = 'txt_spectra'
image_spectra_dir = 'image_spectra'

meta_keys = ['sample','focus_method','distance','exposure','attenuation','ID']
#function to return files in a directory
def filesInDirectory(my_dir: Path, ext: str = '*.tif')->List[str]:
    return list(my_dir.glob(ext))

#function comparing two lists
def listComparison(OriginalList: List[str], NewList: List[str]):
    return list(set(NewList) - set(OriginalList))

def setup_out_dir(out_dir: Path)->None:
    out_dir.mkdir(exist_ok=True)
    
def copy_file(files: Path ,out_dir: Path):
    setup_out_dir(out_dir)
    out_file = out_dir / file.name
    if not out_file.is_file():
        shutil.copy(file,out_file)

def log_scale_image(image: np.array,data_size:int=np.iinfo(np.int32).max):
    # Apply log transformation 
    c = 10*data_size / np.log(1 + np.max(image)*0.2)
    offset=np.abs(np.min(image))+ 1
    #image c * (np.log(image+offset))
    
    return c * (np.log(image+offset))

def compress_image(image: np.array):
    #Returns 8bit normalized image
    return cv2.normalize(image,None,0,255,cv2.NORM_MINMAX,dtype=cv2.CV_8U)

def equalize_image(image: np.array):
    return cv2.equalizeHist(image)

def write_image(image: np.array,file:Path)->None:
    cv2.imwrite(str(file),image)

def image_integrate(ai,img_array,mask,npoints):
    ttheta,intensity,sigma = ai.integrate1d_ng(img_array,
                            npoints,
                            mask=mask,
                            radial_range=[2.5,12],
                            method=("bbox", "csr", "cython"),
                            unit="2th_deg",
                            error_model="Poisson")

    return ttheta,intensity,sigma

def save_spectra_image(ttheta,intensity,wavelength,file: Path):
    plt.plot(calc_dspacing(ttheta,wavelength), np.sqrt(intensity),linestyle='-', marker='+',markersize=5, color='k')
    plt.xlim([0.5,3])
    plt.ylim([0,60])
    plt.xlabel('d-space (A)')
    plt.ylabel('sqrt(Intensity) arb units')
    plt.show()
    plt.savefig(file,format='png',dpi=100)

def write_spectra(ttheta,intensity,file,sigma=None):
	if sigma is None:
		np.savetxt(file,
			np.c_[ttheta,intensity],
			delimiter='\t')
			#header='2theta,Intensity')
	else:
		np.savetxt(file,
			np.c_[ttheta,intensity,sigma],
			delimiter='\t',
            header=f"{file}\nBANK 1 {len(ttheta)} {len(ttheta)} CONS {ttheta[0]} {ttheta[1]-ttheta[0]} 0 0 FXYE",
            comments='')
def calc_dspacing(ttheta: np.array,wavelength: float):
    return 1e10*wavelength/(2*np.sin(np.radians(ttheta)/2))    

#function that performs tasks when new files are detected
def process_image(file: str):
    #load data
    img = fabio.open(str(file))
    img_array = img.data
    
    file = Path(file)
    parent = file.parent
    #Compressed and scaled image
    filename = file.parent / image_log_dir / (Path(file).stem + '.png')
    if not filename.parent.is_dir():
        filename.parent.mkdir(exist_ok=True)
    write_image(
        compress_image(
            log_scale_image(img_array)),str(filename))
    filename = file.parent / image_dir / (Path(file).stem + '.png')
    if not filename.parent.is_dir():
        filename.parent.mkdir(exist_ok=True)
    write_image(
        equalize_image(
            compress_image(img_array)),str(filename))

    #Integration
    ttheta,intensity,sigma = image_integrate(ai,img_array,mask.data,npoints=2500)

    #filename = file.parent / image_spectra_dir / (Path(file).stem + '.png')
    #if not filename.parent.is_dir():
    #    filename.parent.mkdir(exist_ok=True)
    #save_spectra_image(ttheta,intensity,ai.wavelength,str(filename))

    filename = file.parent / spectra_dir / (Path(file).stem + '.txt')
    if not filename.parent.is_dir():
        filename.parent.mkdir(exist_ok=True)
    write_spectra(ttheta,intensity,str(filename))

    filename = file.parent / spectra_dir / (Path(file).stem + '.fxye')
    write_spectra(ttheta,intensity,str(filename),sigma)

def rebuild_cinema():
    cinema_d = []
    files = watch_dir.rglob("*.tif")
    for file in files:
        rel_path = file.parent.relative_to(watch_dir)
        file_d={}
        base_name = file.stem
        #for key,val in zip(meta_keys,base_name.split('_')):
        file_d['id']=base_name.split('_')[-1]
        file_d["FILE_pattern"]=str( rel_path / Path(image_dir) / base_name) + '.png'
        file_d["FILE_pattern_log"]=str( rel_path / Path(image_log_dir) / base_name) + '.png'
        #file_d["FILE_spectra"]=str( rel_path /Path(image_spectra_dir) /  base_name) + '.png'
        file_d["FILE_spectra_txt"]=str( rel_path / Path(spectra_dir)  / base_name) + '.txt'
        #file_d["FILE_pattern_raw"]=str(  base_name) + '.tif'
        cinema_d.append(file_d)
    cinema_df = pd.DataFrame.from_dict(cinema_d)             
    cinema_df.to_csv(str(watch_dir /'data.csv'), index=False, na_rep='NaN')

#function to watch directory

def file_watcher(my_dir: Path, pollTime: int = 5, process_all: bool = False):
    pool = Pool(NCPUS)
    while True:
        time.sleep(pollTime)

        if Path("process_log.csv").is_file() and not process_all:
            df = pd.read_csv("process_log.csv")
            log  = list(df["0"])
        else:
            log=[]

        fileDiff = listComparison(log, [str(file) for file in watch_dir.rglob('*.tif')])   
        if len(fileDiff) != 0:
            #process_image(fileDiff[0])
            pool.map(process_image, fileDiff)
            rebuild_cinema()
            d = {"0": log+fileDiff}
            log_df = pd.DataFrame.from_dict(d)             
            log_df.to_csv('process_log.csv', index=False)
        print("Done for now")

if __name__ == "__main__":
    #file_watcher(my_dir=watch_dir,process_all=False)
    freeze_support()
    #rebuild_cinema()
    file_watcher(watch_dir,process_all=False)
