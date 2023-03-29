#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:44:16 2020
Updated on Mon Feb 20 16:19:33 2023

@author: jgra
"""

# In[1]
# Import Libraries and Tools
import numpy as np

# In[2]
# Define Save function
def saveSeg(filename,seg,h,correct_limb_brightening,resize_param,
            foreground_weight,background_weight,init_mask,init_mask_method,
            fill_init_holes,init_alpha,alpha,narrowband,N,
            image_preprocess=None):
    '''
    Save Function for use in collaboration with ACWE output generated using 
    acweFunctions_v4 or greater.
    
    Parameters
    ----------
    filename : str
        Full File Path where final Segmentation will be stored
    seg : [bool] OR [float]
        [MxN] Boolean array of single ACWE segmentation OR [IxMxN] 
        floating point array of multiple segmentations in the case of a 
        confidence map
        
        I : number of segmentations
        
        M, N : segmentation Dimensions
    h : dict
        Full header of original .fits file used to generate this
        segmentation. It should be noted that the header can be provided 
        as generated using the fits.open protocol from astropy.io. This
        function will, however, convert the output to a dictionary before
        saving it within the metadata.
    correct_limb_brightening : bool
        Set to True to indicate that the limb brightening correction of 
        [2] was performed on the .fits file
    resize_param : int
        The factor by which the image was downsampled prior to performing
        ACWE. 
    foreground_weight : float
        Weight term for the foreground (CH) homogeneity within the ACWE 
        energy functional.
    background_weight : float
        Weight term for the background (quiet Sun and all remaining on 
        disk features) homogeneity within the ACWE energy functional.
    init_mask : [bool]
        initial mask prior to performing ACWE
    init_mask_method : str
        Text description of the method used to generate initial mask to 
        facilitate future methods of generating this mask.
        
        For consistency, please use the pusdo-code string "alpha*mean(qs)" 
        for the traditional method (described in [1]).
        
        If another method is used to generate the initial mask, please
        provide either the name or a written description of the process
        here.
    fill_init_holes : bool
        State if holes in the initial mask are filled prior to performing ACWE
        
        if True holes were filled, if False holes were not
    init_alpha : float
        An alpha parameter is used as a threshold to generate the initial 
        mask. This is the initial choice for the alpha parameter.
        
        Use NaN when the initial mask generation method does not use a 
        thresholding method.
    alpha : float
        For the traditional method of mask generation, there remains the 
        possibility that the initial mask will fail. In acweFunctions_v4
        there is the option to mitigate for this issue by incrementally
        updating alpha until a successful segmentation is generated, as 
        such place here the actual alpha (threshold) parameter used to 
        generate the final mask.
        
        When init_alpha and alpha are the same, indicate so by populating
        this space with the same value as init_alpha.
        
        Use NaN when the initial mask generation method does not use a 
        thresholding method.
    narrowband : int 
        Constraint on ACWE evolution to ensure iterative optimization 
        process does not result in overcorrection of contour boundary.
    N : int
        Number of iterations of ACWE between checks for convergence.
    image_preprocess : str, optional
        String describing any additional processing done to the original
        solar EUV image prior to performing ACWE.
        
        Default Value: None
    Outputs
    -------
    At the specified file path there will be a compressed .npz file
    containing, in order the original header h, as a dictionary, a header
    AH outlining the ACWE Process, and the ACWE segmentation(s) seg. 
    Since two of the outputs are headers, reopening the files requires
    that the allow_pickle setting within the numpy load function to be True.
    
    References
    ----------
    [1] 
        L. E. Boucheron, M. Valluri, and R. T. J. McAteer, "Segmentation 
        of Coronal Holes Using Active Contours Without Edges," Solar 
        Physics, vol. 291, pp. 2353-2372, 2016.
    [2] 
        C. Verbeeck, V. Delouille, B. Mampaey, & R. De Visscher, "The 
        SPoCA-suite: Software for extraction, characterization and 
        tracking of active regions and coronal holes on EUV images," 
        Astronomy & Astrophysics, vol. 561, pp. A29, 2014.
    '''
    
    # convert original header into dictionary
    H = dict(h)
    
    # Record key elements of ACWE process in second header
    segHeader = {
                'CORRECT_LIMB_BRIGHTENING' : correct_limb_brightening,
                'IMAGE_PREPROCESS'         : image_preprocess,
                'RESIZE_PARAM'             : resize_param,
                'FOREGROUND_WEIGHT'        : foreground_weight,
                'BACKGROUND_WEIGHT'        : background_weight,
                'INIT_MASK'                : init_mask,
                'INIT_MASK_METHOD'         : init_mask_method,
                'Fill_INIT_HOLES'          : fill_init_holes,
                'INIT_ALPHA'               : init_alpha,
                'ALPHA'                    : alpha,
                'NARROWBAND'               : narrowband,
                'ITTER_BETWEEN_CHK'        : N
                }
    
    # save as .npz file
    np.savez_compressed(filename,H,segHeader,seg)

# In[3]:
# Define open function 
def openSeg(filename):
    '''
    Open final acwe segmentation file function. This function is fully
    compatible with segmentations saved using acweSaveSeg_v2.py and 
    later.

    Parameters
    ----------
    filename : str
            Full File Path where final Segmentation will be stored

    Returns
    -------
    FITSHEADER : dict
        Full header of original .fits file used to generate this
        segmentation. 
    ACWEHEADER : dict
        Header outlining ACWE process. Refer to the saveSeg function to see 
        information about data within header.
    SEG : [bool] OR [float]
        ACWE Segmentation
    '''
    
    # Open .npz file and get list of "arrays"
    data = np.load(filename, allow_pickle=True)
    lst = data.files
    
    # Separate into Header and Image
    FITSHEADER = data[lst[0]].tolist() # Header of Original .fits File
    ACWEHEADER = data[lst[1]].tolist() # Record key elements of ACWE process
    SEG = data[lst[2]] # Segment
    
    # Return results
    return FITSHEADER, ACWEHEADER, SEG
