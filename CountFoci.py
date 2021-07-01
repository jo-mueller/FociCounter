# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 15:13:21 2021

@author: johan
"""

import os
import pandas as pd
import numpy as np
import napari
from tqdm import tqdm
from aicsimageio import AICSImage
import tifffile as tf
import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
from skimage.morphology import binary_erosion
from skimage.segmentation import expand_labels
from skimage import segmentation as skseg
from skimage import restoration
from skimage.feature import blob_dog
from scipy import ndimage as ndi
import copy
import matplotlib.pyplot as plt


def findFoci(FociMap, LabelMap, **kwargs):
    "Find local maxima on foci image"
    
    size = kwargs.get('size', 20)
    coordinates = peak_local_max(FociMap, min_distance=size)
    
    
def readFile(filename):
    
    if filename.endswith('tiff') or filename.endswith('tif'):
        img = tf.imread(filename)
        
        if np.argmin(img.shape) == 2:
            img = img.transpose((2,0,1))
        
        DAPI = img[2]
        Foci = img[1]
    else:
        img = AICSImage(filename)
        img = img.get_image_data("CYX", Z=0, T=0, S=0, B=0)
        DAPI = img[1]
        Foci = img[0]
        
    return DAPI, Foci

def runStarDIst(img):
    
    from stardist.models import StarDist2D 
    from csbdeep.utils import normalize

    # prints a list of available models 
    StarDist2D.from_pretrained() 
    
    # creates a pretrained model
    model = StarDist2D.from_pretrained('2D_versatile_fluo')
    
    labels, _ = model.predict_instances(normalize(img), prob_thresh=0.2)
    
    return labels

def filter_cells(labels, focimap, **kwargs):
    
    minsize = kwargs.get('Minsize', 15)
    pixsize = kwargs.get('PixSize', 0.1135)
    
    skseg.clear_border(labels)
    Celllist = np.unique(labels)
    Sizes = [np.sum(labels == x)*pixsize**2 for x in Celllist]
    

    
    df = pd.DataFrame(columns=['Cell ID', 'Size', 'Intensity'])
    df['Cell ID'] = Celllist
    df['Size'] = Sizes
    
    # collect intensity values of nuclei in foci-channel and add to df
    intensity = []
    for i, sample in df.iterrows():
        mean = np.mean(focimap[labels == sample['Cell ID']])
        intensity.append(mean)
        
    df['Intensity'] = intensity
    threshold = np.quantile(df['Intensity'], 0.9)

    # Iterate over nuclei and check criteria
    for i, sample in df.iterrows():
        
        if sample.Intensity > threshold:
            labels[labels == sample['Cell ID']] = 0
        
        if sample.Size < minsize:            
            labels[labels == sample['Cell ID']] = 0
    
    return labels

def process_image(DAPI, Foci):
    
    # Preprocess
    Foci = Foci - restoration.rolling_ball(Foci, radius=20)
    
    # segment nuclei
    labelmap = runStarDIst(DAPI)
    labelmap = filter_cells(labelmap, Foci, Minsize=15)
    
    # FInd foci
    nucleus_map = labelmap != 0
    df = pd.DataFrame(columns=['Nucleus ID'])
    df['Nucleus ID'] = np.unique(labelmap)
    df = df[df['Nucleus ID'] != 0]

    # First, find all foci
    foci_map = copy.deepcopy(Foci)
    foci = peak_local_max(foci_map, min_distance=5, threshold_abs=np.quantile(Foci, 0.9))
    # foci = peak_local_max(foci_map,
    #                       labels=labelmap,
    #                       min_distance=10,
    #                       threshold_abs=np.quantile(Foci, 0.8),
    #                       footprint=np.ones((10, 10)))
    
    df_foci = pd.DataFrame(columns=['x', 'y', 'label', 'intensity', 'Size'])
    
    df_foci.x = [x[0] for x in foci]
    df_foci.y = [x[1] for x in foci]
    
    point_labels = []
    point_intensity = []
    for i, pt in enumerate(foci):
        point_labels.append(labelmap[int(pt[0]), int(pt[1])])
        point_intensity.append(Foci[int(pt[0]), int(pt[1])])
        
    df_foci.label = point_labels
    df_foci.intensity = point_intensity
    df_foci = df_foci[df_foci.label != 0].sort_values(by='label')
    
    # iterate over cells
    sizes = []
    n_foci = []
    for label in df['Nucleus ID'].unique():
        sizes.append(np.sum(labelmap == label))
        n_foci.append((df_foci.label == label).sum())
    
    df['Nucleus Size'] = sizes
    df['N foci'] = n_foci
    
    # iterate over cells
    proximity_map = np.zeros_like(Foci, dtype=int)
    perimeter_map = np.zeros_like(Foci, dtype=int)
    for i, sample in tqdm(df.iterrows()):
        
        label = sample['Nucleus ID']
        
        # Make a copy of this nucleus
        x1 = np.argwhere(labelmap == label)[:,0].min()+1
        x2 = np.argwhere(labelmap == label)[:,0].max()+1
        y1 = np.argwhere(labelmap == label)[:,1].min()+1
        y2 = np.argwhere(labelmap == label)[:,1].max()+1
        
        Nucleus_DAPI = DAPI[x1:x2, y1:y2]
        Nucleus_Foci = Foci[x1:x2, y1:y2]
        Nucleus_labelmap = labelmap[x1:x2, y1:y2]
        
        # get local background
        BG = np.quantile(Nucleus_Foci, 0.1)
        
        # get those foci that correspond to this cell a.k.a. label
        _df_foci = df_foci[df_foci.label == label]
        
        # put foci locations on empty image as seeds
        _proximity_map = np.zeros_like(Nucleus_DAPI, dtype='uint16')
        _proximity_map[_df_foci.x-x1, _df_foci.y-y1] = _df_foci.index
        
        # expand seeds and crop with cell label region
        _proximity_map = expand_labels(_proximity_map, distance=40)
        _proximity_map[Nucleus_labelmap != label] = 0
        
        # calc. FAHM for every foci's proximity zone (full area half maximum)
        _perimeter_map = np.zeros_like(_proximity_map, dtype='uint16')
        for j, fc in _df_foci.iterrows():
            HM = BG + 0.5*(fc.intensity - BG)  # Half-maximum above local background
            local_foci_map = skseg.flood(Nucleus_Foci, (fc.x-x1, fc.y-y1),
                                         tolerance=(fc.intensity - HM),
                                         connectivity=8) * fc.name
            local_foci_map[_proximity_map != fc.name] = 0
            df_foci.loc[j, 'Size'] = np.sum(local_foci_map)
            _perimeter_map += local_foci_map.astype('uint16')
        
        proximity_map[x1:x2, y1:y2] += _proximity_map
        perimeter_map[x1:x2, y1:y2] += _perimeter_map
        
    return df, df_foci
    
def plot_image_data(df_c, df_f):
    fig, axes = plt.subplots(nrows=1, ncols=3)
    
    axes[0].hist(df_c['Nucleus Size'], bins=30)
    axes[0].set_xlabel('Nucleus Size')
    axes[0].set_ylabel('N')
    
    axes[1].hist(df_c['N foci'], bins=30)
    axes[1].set_xlabel('Foci per cell')
    axes[1].set_ylabel('N')
    
    axes[2].hist(df_f['Size'], bins=30)
    axes[2].set_xlabel('Foci size')
    axes[2].set_ylabel('N')
    
    

if __name__ == "__main__":
    directory = r'D:\Documents\Promotion\Projects\FociCounter\Data'
    
    files = os.listdir(directory)
    files = [os.path.join(directory, x) for x in files if (x.endswith('tiff') or x.endswith('czi'))]
    df = pd.DataFrame(data=files, columns=['Filename'])
    
    DAPI, Foci = readFile(df.Filename.loc[1])    
    _df_cells, _df_foci = process_image(DAPI, Foci)
    
    plot_image_data(_df_cells, _df_foci)