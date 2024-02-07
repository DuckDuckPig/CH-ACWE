#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    test_acwe.py from "ACWE_python_v3" converted into callable functions.

Created on Mon Jun  8 14:28:32 2020
Edited on  Thu Apr  1 14:25:54 2021 - Adjusted for STEREO A and B
           Wed Jul  7 12:14:39 2021 - Rolling Alpha Correction,
                                      Optimized Confidence Map Implementation
           Wed Apr 27 09:59:53 2022 - Return Init Mask, Documentation
           Thu May 26 16:52:18 2022 - Return Init Mask Without Hole Filling
           Mon Feb 20 14:42:47 2023 - Condenced to tools used in paper, with 
                                      aditinal documentation
           Wed Feb  7 12:11:59 2024 - Corrected bug in stopping critera

@author: jgra
"""

# In[1]
# Import Libraries and tools
import numpy as np
from matplotlib import pyplot as plt
import skimage.transform
import copy
from .ACWE_python_v3 import correct_limb_brightening
import scipy as sp
from .ACWE_python_v3 import acwe

# In[2]:
# Resizing Function
def resize_EUV(J,h,resize_param=8,interpolation='Bi-cubic'):
    '''
    Resizes solar EUV image based on user-specified resize parameter, and 
    return useful metadata about resized image.
    
    Parameters
    ----------
    J : [float]
        Solar EUV image
    h : dict
        original .fits header
    resize_param : int
        Resize Parameter which defines the degree to which the image should be
        spatially downsampled. Note that the parameter can be interpreted as 
        the denominator of a fraction such that a value of '2' would indicate 
        that the output image should be 1/2 the original spatial resolution 
        in each dimension, 4 indicates 1/4 scale, etc. By convention, this 
        parameter is 8 for SDO-AIA images and 4 for STEREO A and STEREO B 
        images, this produces an output that is 512x512 pixels in size.
    interpolation : str, optional
        Interpolation method for the resizing process. Valid options are 
        'Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
        'Bi-quartic', and 'Bi-quintic'.
        
        Default Value: 'Bi-cubic'
    Returns
    -------
        I : [float]
            Solar EUV Image, resized to user-specified dimensions.
        im_size : [int]
            Array with provides the dimensions of the image.
        sun_radius : float
            Radius of the Sun in the resized image.
        sun_center : [int]
            Coordinates of the center of the sun in the resized image.
    '''
    
    # Determine Downsample Method
    downsample = ['Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
                  'Bi-quartic','Bi-quintic']
    for order in range(len(downsample)):
        if downsample[order] == interpolation:
            break
    
    # Resize image
    if resize_param > 1:
        I = skimage.transform.resize(J,np.asarray(J.shape)/resize_param,
                                     order=order,preserve_range=True,
                                     anti_aliasing=True)
    else:
        I = copy.deepcopy(J)
        
    # Determine characteristics of image
    im_size = np.asarray(I.shape) # size of the image
    
    # Get Solar Radius from Metadata
    try: # SDO/AIA
        sun_radius = h['R_SUN']/resize_param
    except: # STERO A or B
        sun_radius = (h['RSUN']/h['CDELT1'])/resize_param
    
    # Get Solar Center from Metadata
    sun_center = np.asarray([int(round(h['CRPIX1']))-1,
                             int(round(h['CRPIX2']))-1])
    sun_center = sun_center/resize_param
    
    # Return Resized Image, Image Dimensions and Solar Radius & Center
    return I,im_size,sun_radius,sun_center


# In[3]
# Masking Functions

# Circle Mask
def make_circle_mask(c,im_dims,r):
    '''
    Defines a binary image with image dimensions im_dims of a circle with 
    center c and radius r
    
    Parameters
    ----------
    c : [int,int]
        x, and y coordinates of the center of a circle, obeying right-hand 
        rule where the upper left corner is the origin.
    im_dims : [int,int]
        Image dimensions in format [x,y] where x is the height and y is the 
        length
    r : float
        radius of circle
    Returns
    -------
    c_mask : [bool]
        Mask of size im_dims where circle of radius r, centered at c, is given
        the value of 1 and all other regions are assigned a value of 0.
    '''
    cx = c[0]
    cy = c[1]
    ix = im_dims[0]
    iy = im_dims[1]
    x,y = np.meshgrid(np.arange(-(cx),(ix-cx),1),np.arange(-(cy),(iy-cy),1))
    c_mask = (x**2+y**2)<=r**2
    return c_mask

# Initial Masks
def inital_masks(I,im_size,sun_radius,sun_center,alpha=0.3,rollingAlpha=0):
    '''
    Function returns circle mask that separates on-disk and off disk areas
    and initial mask for performing ACWE.
    
    Parameters
    ----------
    I : [float]
            Solar EUV Image, resized to user-specified dimensions.
    im_size : [int]
            Array with provides the dimensions of the image.
    sun_radius : float
        Radius of the Sun in the resized image.
    sun_center : [int]
        Coordinates of the center of the sun in the resized image.
    alpha : float, optional
        Threshold parameter as expressed in [1], alpha will be multiplied by 
        mean quiet Sun intensity (QS) to generate the threshold for the 
        initial mask. It is recommended that an alpha of 0.3 be use for AIA.
        An optimal value for STEREO is still being determined.
        
        Default Value: 0.3
    rollingAlpha : float, optional
        If no mask is produced using threshold alpha * QS, the threshold will
        be replaced with (alpha + i*rollingAlpha) * QS where 'i' is the number
        of times an empty mask was produced. This will result in a mask being 
        produced, even if no CH is present on disk. Set to 0 to disable this
        process.
        
        Default Value: 0 (Mask will always be alpha * QS)
    Returns
    -------
    sd_mask : [bool]
        Mask that separates on-disk and off-disk areas
    m : [bool]
        Initial mask, note holes are not filled prior to returning mask.
    alphar : float, optional
        The alpha parameter that was actually used to generate the mask, only 
        returned if rollingAlpha != 0
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
    
    # Define solar disk mask
    sd_mask = make_circle_mask(sun_center,im_size,sun_radius)
    
    # Determine threshold value for initialization of AC as percentage of QS;  
    # estimate QS as maximum bin of histogram
    Ihist,ihist_edges = np.histogram(I[sd_mask],100) # 100 bin histogram of SD
    ihist_centers = (ihist_edges[:-1]+ihist_edges[1:])/2. # bin centers
    ihistmax = Ihist.argmax() # index of maximum bin
    initial_thresh = alpha*ihist_centers[ihistmax] # bin center of maximum bin
    m = (I<=initial_thresh)*sd_mask # initial mask
    
    if rollingAlpha != 0:
        # Save old alpha parameter for exception check
        alphar = alpha * 1
    
        # Adjust Alpha Parameter, if needed and instructed 
        # to do so to ensure valid initial mask
        while np.sum(m.astype(int)) == 0:
            alphar += 0.01 # increase alpha by 1% and regenerate initial mask
            initial_thresh = alphar*ihist_centers[ihistmax] # bin center of maximum bin
            m = (I<=initial_thresh)*sd_mask # new initial mask
        
    # Return Results
        return sd_mask,m,alphar
    else:
        return sd_mask,m

# In[4]
# ACWE Iterations

# Single ACWE Segmentation
def itterate_acwe(I,im_size,sd_mask,m,foreground_weight=1,
                  background_weight=1/50.,narrowband=2,N=10,
                  fillInitHoles=True,verbose=False):
    '''
    Runs coronal hole (CH) segmentation using active contours without edges 
    (ACWE) as described in [1].
    
    Parameters
    ----------
    I : [float]
        Solar EUV Image, resized to user-specified dimensions, with correction
        for limb brightening, if needed.
    im_size : [int]
        Array which provides the dimensions of the image.
    sd_mask : [bool]
        Mask that separates on-disk and off-disk areas.
    m : [bool]
        Inital mask.
    foreground_weight : float, optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of 
        the backround_weight compared to forground_weight.
        
        Default Value: 1
    background_weight : float, optional
        Weight term for the background (quiet Sun and all remaining on disk
        features) homogeneity within the energy functional.
        
        Default Value: 1/50.0
    narrowband : int, optional
        Constraint on ACWE evolution to ensure iterative optimization process
        does not result in overcorrection of contour boundary.
        
        Default Value: 2
    N : int, optional
        Number of iterations of ACWE between checks for convergence
        
        Default Value: 10
    fillInitHoles : bool, optional
        Fill holes in initial mask using morphology prior to performing ACWE
        
        Default Value: True
    verbose : bool, optional
        Display ACWE evolution in real time.
        
        Default Value: False
    Returns
    -------
    seg : [bool]
        final segmentation mask, in dimensions 
        np.asarray(J.shape)/resize_param
    References
    ----------
    [1] 
        L. E. Boucheron, M. Valluri, and R. T. J. McAteer, "Segmentation 
        of Coronal Holes Using Active Contours Without Edges," Solar
        Physics, vol. 291, pp. 2353-2372, 2016.
    '''
    
    if verbose:
        plt.ion() # interactive plotting on 
    
    # Set up variables for ACWE iterations
    if fillInitHoles:
        m_seg = sp.ndimage.morphology.binary_fill_holes(m) # fill holes
    else:
        m_seg = m # image to keep track of current initialization of ACWE
        
    # Valid Mask - Perform ACWE
    if np.sum(m.astype(int))!=0:
        
        I_seg = I # image of current segmentation
        I_seg[~sd_mask] = I[~m_seg&sd_mask].mean() # set pixels outside SD to
                                                   # mean of background to 
                                                   # force ACWE to ignore
        counter = 0 # to keep track of proxy of iterations
        seg_diff_cum = np.zeros(im_size) # to keep track of how many times
                                         # pixels change classes over
                                         # iterations
        iterate = 1 # flag to continue iterating
        
        # Continue iterating with N iterations of the ACWE evolution, 
        # followed by check for convergence. Running ACWE for default of N=10
        # iterations is a good trade-off between checking too often and not
        # often enough for convergence.
        
        if verbose:
            print('% Diff            % New Diff')
        while iterate:
            seg = acwe.acwe(I_seg,m_seg,N,(0,0,foreground_weight,
                            background_weight),narrowband,verbose) # Evolve ACWE 
                                                                   # for N
                                                                   # Iterations
            
            # update current segmentation
            I_seg = I 
            I_seg[~sd_mask] = I[~seg&sd_mask].mean() 
            
            # difference in seg from previous iteration to now
            seg_diff = seg.astype(int) - m_seg.astype(int) 
            
            m_seg = seg # update m_seg image
            counter = counter + 1 # iterate counter
            
            # compute percentage of pixels that changed between previous 
            # iteration, and now
            percent_diff = float(abs(seg_diff).sum())/float(seg.sum())*100.0
            # keep track of how many times pixels have changed classes
            seg_diff_cum = seg_diff_cum + abs(seg_diff)
            # percentage of currently new pixels that have never changed
            # classes before
            percent_new_diff = float(((seg_diff_cum==1)*abs(seg_diff)).sum())/\
                               float(((seg_diff_cum>=1)*abs(seg_diff)).sum()+\
                               np.finfo(float).eps)*100 
            if verbose:
                print(str(percent_diff) + ' ' + str(percent_new_diff))
            if percent_new_diff==0 | ~(seg.sum()>0):
                iterate = 0
        
        # Return Segmentation
        return seg
    else:
        return m * 1

# ACWE Confidence Maps
def itterate_acwe_confidence_map(I,im_size,sd_mask,m,foreground_weight=1,
                                 background_weights=[1/50.],narrowband=2,
                                 N=10,fillInitHoles=True,verbose=False):
    '''
    Runs coronal hole (CH) segmentation using active contours without edges 
    (ACWE) as described in [1].
    
    Parameters
    ----------
    I : [float]
        Solar EUV Image, resized to user-specified dimensions, with correction 
        for limb brightening, if needed.
    im_size : [int]
        Array with provides the dimensions of the image.
    sd_mask : [bool]
        Mask that separates on-disk and off-disk areas
    m : [bool]
        Initial mask.
    foreground_weight : float, optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of the 
        backround_weight compared to forground_weight. This function only
        accpets one foreground weight, which is kept constant through all
        cases.
        
        Default Value: 1
    background_weights : [float], optional
        Weight terms for the background (quiet Sun and all remaining on disk 
        features) homogeneity within the energy functional. These values do 
        not need to be ordered.
        
        Default Value: [1/50.0]
    narrowband : int, optional
        Constraint on ACWE evolution to ensure iterative optimization process 
        does not result in overcorrection of contour boundary.
        
        Default Value: 2
    N : int, optional
        Number of iterations of ACWE between checks for convergence
        
        Default Value: 10
    fillInitHoles : bool, optional
        Fill holes in initial mask using morphology prior to performing ACWE
        
        Default Value: True
    verbose : bool, optional
        Display ACWE evolution in real time.
        
        Default Value: False
    Returns
    -------
    Segs : [float]
        final segmentation mask, in dimensions
        np.hstack([len(background_weight),np.asarray(J.shape)/resize_param])
    References
    ----------
    [1] 
        L. E. Boucheron, M. Valluri, and R. T. J. McAteer, "Segmentation 
        of Coronal Holes Using Active Contours Without Edges," Solar 
        Physics, vol. 291, pp. 2353-2372, 2016.
    '''
    
    if verbose:
        plt.ion() # interactive plotting on 
    
    # Generate ordered list of background weights
    background_weight_ordered = np.unique(np.sort(background_weights))
    # Doubles are neither expected nor recommended, however accounting for this
    # now should further optimize this implementation of this function for the
    # case where a user may wish to use this method to weigh different 
    # background_rates
    
    # Set up variables for ACWE iterations
    if fillInitHoles:
        m_seg = sp.ndimage.morphology.binary_fill_holes(m) # fill holes
    else:
        m_seg = m # image to keep track of current initialization of ACWE
    
    # Prepare for ACWE
    outputShape = np.asarray(m_seg.shape)
    outputShape = np.hstack([len(background_weights),outputShape]).astype(int)
    Segs = np.empty(outputShape); Segs[:] = np.nan
    
    # Valid Mask - Perform ACWE
    if np.sum(m.astype(int))!=0:
    
        # evolve from smallest to largest background weight
        # Under the assumption that iterating outward from the initial threshold
        # will result in the less iterations that expanding outward to the largest
        # segmentation then working back in again, this method was chosen to 
        # further optimize runtime.
        for background_weight in background_weight_ordered:
            
            # Finish set up of variables for ACWE iterations
            I_seg = I # image of current segmentation
            I_seg[~sd_mask] = I[~m_seg&sd_mask].mean() # set pixels outside SD to
                                                       # mean of background to
                                                       # force ACWE to ignore
            counter = 0 # to keep track of proxy of iterations
            seg_diff_cum = np.zeros(im_size) # to keep track of how many times
                                             # pixels change classes over
                                             # iterations
            iterate = 1 # flag to continue iterating
            
            # continue iterating with N iterations of the ACWE evolution, followed
            # by check for convergence
            # running ACWE for default of N=10 iterations is a good trade-off
            # between checking too often and not often enough for convergence
            if verbose:
                print('% Diff            % New Diff')
            while iterate:
                seg = acwe.acwe(I_seg,m_seg,N,(0,0,foreground_weight,
                                background_weight),narrowband,verbose) # Evolve
                                                                       # ACWE for N
                                                                       # Iterations
                
                # update current segmentation
                I_seg = I 
                I_seg[~sd_mask] = I[~seg&sd_mask].mean() 
                
                # difference in seg from previous iteration to now
                seg_diff = seg.astype(int) - m_seg.astype(int) 
                
                m_seg = seg # update m_seg image
                counter = counter + 1 # iterate counter
                
                # compute percentage of pixels that changed between previous
                # iteration and now
                percent_diff = float(abs(seg_diff).sum())/float(seg.sum())*100.0
                # keep track of how many times pixels have changed classes
                seg_diff_cum = seg_diff_cum + abs(seg_diff)
                # percentage of currently new pixels that have never changed
                # classes before
                percent_new_diff = float(((seg_diff_cum==1)*abs(seg_diff)).sum())/\
                                   float(((seg_diff_cum>=1)*abs(seg_diff)).sum()+\
                                   np.finfo(float).eps)*100 
                if verbose:
                    print(str(percent_diff) + ' ' + str(percent_new_diff))
                if percent_new_diff==0 | ~(seg.sum()>0):
                    iterate = 0
            
            # Find and fill appropriate background weight index/indices
            index = np.where(background_weight==background_weights)[0]
            for i in index:
                Segs[i] = seg.astype(int) * 1
                
        # Return Segmentation
        return Segs
    else:
        Segs[:]=0
        return Segs

# In[5]
# Running ACWE
def run_acwe(J,h,resize_param=8,foreground_weight=1,background_weight=1/50.,
             alpha=0.3,narrowband=2,N=10,verbose=False,
             correctLimbBrightening=True,rollingAlpha=0,fillInitHoles=True):
    
    '''
    Primary function for running coronal hole (CH) segmentation using active 
    contours without edges (ACWE) as described in [1], including option for 
    limb brightening correction of [2].
    
    Parameters
    ----------
    J : [float]
        Solar EUV image stored as a numpy array
    h : dict
        .fits header for Solar EUV image J
    resize_param : int, optional
        The factor by which the image will be downsampled. In general
        operation ACWE usually operates over a 512X512 pixel image,
        thus it is recommended that 
        resize_param = np.min(np.asarray(J.shape)/np.array([512,512])).astype(int)
        
        Default Value: 8
    foreground_weight : float, optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of 
        the backround_weight to forground_weight.
        
        Default Value: 1
    background_weight : float, optional
        Weight term for the background (quiet Sun and all remaining on disk
        features) homogeneity within the energy functional.
        
        Default Value: 1/50.0
    alpha : float, optional
        Threshold parameter alpha will be multiplied by mean quiet sun 
        intensity to generate the threshold for the initial mask.

        Default Value: 0.3
    narrowband : int, optional
        Constraint on ACWE evolution to ensure iterative optimization process 
        does not result in overcorrection of contour boundary.
        
        Default Value: 2
    N : int, optional
        Number of iterations of ACWE between checks for convergence.
        
        Default Value: 10
    verbose : bool, optional
        Display ACWE evolution in real time.
        
        Default Value: False
    correctLimbBrightening : bool, optional
        Perform limb brightening correction of [2]. Note, process is not 
        recommended on STEREO images.
        
        Default Value: True
    rollingAlpha : float, optional
        If no mask is produced using threshold alpha * QS, the threshold will
        be replaced with (alpha + i*rollingAlpha) * QS where 'i' is the number
        of times an empty mask was produced. This will result in a mask being 
        produced, even if no CH is present on disk. Set to 0 to disable this
        process.
        
        Default Value: 0 (Mask will always be alpha * QS)
    fillInitHoles : bool, optional
        Fill holes in initial mask
        
        Default Value: True
        
    Returns
    -------
    seg : [bool]
        final segmentation mask, in dimensions 
        np.asarray(J.shape)/resize_param
    alphar : float, optional
        the final alpha parameter, returned if (and only if)
        oldThreshold == True and rollingAlpha == True
    m : [bool]
        Initial mask without any holes filled
    
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
    
    # Resize image
    I,im_size,sun_radius,sun_center = resize_EUV(J,h,resize_param)
    
    # Correct limb brightening per Verbeeck et al. 2014
    if correctLimbBrightening:
        I = correct_limb_brightening.correct_limb_brightening(I,sun_center,
                                                              sun_radius)

    #  Define solar disk mask and initial mask
    if rollingAlpha != 0:
        sd_mask,m,alphar = inital_masks(I,im_size,sun_radius,sun_center,
                                        alpha,rollingAlpha)
    else:
        sd_mask,m = inital_masks(I,im_size,sun_radius,sun_center,alpha,
                                 rollingAlpha)
    
    # Perform ACWE
    seg = itterate_acwe(I,im_size,sd_mask,m,foreground_weight,
                        background_weight,narrowband,N,fillInitHoles,verbose)
    
    # Return Results
    if rollingAlpha != 0:
        return seg,alphar,m
    
    else:
        return seg,m

# ACWE Confidence Map
def run_acwe_confidenceMap(J,h,resize_param=8,foreground_weight=1,
                           background_weights=[1/50.],alpha=0.3,narrowband=2,
                           N=10,verbose=False,correctLimbBrightening=True,
                           rollingAlpha=0,fillInitHoles=True):
    
    '''
    Function for generating confidence map based segmentation of coronal hole 
    (CH) using active contours without edges (ACWE) as described in [1], 
    including option for limb brightening correction of [2] by iterating over
    the background_weight parameter.
    
    Parameters
    ----------
    J : [float]
        Solar EUV image stored as a numpy array
    h : dict
        .fits header for Solar EUV image J
    resize_param : int, optional
        The factor by which the image will be downsampled. In general
        operation ACWE usually operates over a 512X512 pixel image, thus it is
        recommended that 
        resize_param = np.min(np.asarray(J.shape)/np.array([512,512])).astype(int)
        
        Default Value: 8
    foreground_weight : float, optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of the
        backround_weight to forground_weight. This function only accepts one 
        foreground weight, which is kept constant through all cases.
        
        Default Value: 1
    background_weights : [float], optional
        Weight terms for the background (quiet Sun and all remaining on disk
        features) homogeneity within the energy functional. These values do
        not need to be ordered.
        
        Default Value: [1/50.0]
    alpha : float, optional
        Threshold parameter alpha will be multiplied by mean quiet Sun 
        intensity to generate the threshold for the initial mask.

        Default Value: 0.3
    narrowband : int, optional
        Constraint on ACWE evolution to ensure iterative optimization 
        process does not result in overcorrection of contour boundary.
        
        Default Value: 2
    N : int, optional
        Number of iterations of ACWE between checks for convergence
        
        Default Value: 10
    verbose : bool, optional
        Display ACWE evolution in real time.
        
        Default Value: False
    correctLimbBrightening : bool, optional
        Perform limb brightening correction of [2]. Note, process is not 
        recommended on STEREO images.
        
        Default Value: True
    rollingAlpha : float, optional
        If no mask is produced using threshold alpha * QS, the threshold will
        be replaced with (alpha + i*rollingAlpha) * QS where 'i' is the number
        of times an empty mask was produced. This will result in a mask being 
        produced, even if no CH is present on disk. Set to 0 to disable this
        process.
        
        Default Value: 0 (Mask will always be alpha * QS)
    fillInitHoles : bool, optional
        Fill holes in initial mask
        
        Default Value: True
    
    Returns
    -------
    Segs : [float]
        final segmentation mask, in dimensions 
        np.hstack([len(background_weight),np.asarray(J.shape)/resize_param])
        
        Note:
            The relationship between background_weight and this final 
            segmentation is such that Segs[i] is the segmentation generated
            using background_weight[i]
    alphar : float, optional
        the final alpha parameter, returned if (and only if)
        oldThreshold == True and rollingAlpha == True 
    m : [bool]
        Initial mask without any holes filled
    
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
    
    # Resize image
    I,im_size,sun_radius,sun_center = resize_EUV(J,h,resize_param)
    
    # Correct limb brightening per Verbeeck et al. 2014
    if correctLimbBrightening:
        I = correct_limb_brightening.correct_limb_brightening(I,sun_center,
                                                              sun_radius)

    # Define solar disk mask and initial mask
    if rollingAlpha != 0:
        sd_mask,m,alphar = inital_masks(I,im_size,sun_radius,sun_center,
                                        alpha,rollingAlpha)
    else:
        sd_mask,m = inital_masks(I,im_size,sun_radius,sun_center,alpha,
                                 rollingAlpha)
    
    # Return ACWE
    Segs = itterate_acwe_confidence_map(I,im_size,sd_mask,m,foreground_weight,
                                        background_weights,narrowband,N,
                                        fillInitHoles,verbose)
    
    # Return Results
    if rollingAlpha != 0:
        return Segs,alphar,m
    
    else:
        return Segs,m
