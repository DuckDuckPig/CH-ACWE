#-------------------------------------------------------------------------------
# Test Script for Active Contours Without Edges
#
# test_acwe
#
# Script to run segmentation of coronal holes using active contours without
# edges (ACWE) as described in [1], including limb brightening correction of
# [2].  
#
# User can define a base directory under which .fits files reside.  The code
# as written will loop over all .fits files in that base directory and compute
# the ACWE coronal hole segmentation and display the results.  As written, the
# final segmentation is not saved.  The user can add a line to save the result
# in the desired format immediately after exiting the while loop.
#
# References:
# [1] L. E. Boucheron, M. Valluri, and R. T. J. McAteer, "Segementation of 
#     Coronal Holes Using Active Contours Without Edges," Solar Physics, vol.
#     291, pp. 2353-2372, 2016.
# [2] C. Verbeeck, V. Delouille, B. Mampaey, & R. De Visscher, "The SPoCA-
#     suite: Software for extraction, characterization and tracking of active
#     regions and coronal holes on EUV images," Astronomy & Astrophysics, vol.
#     561, pp. A29, 2014.
#
# Notes:
# Ideas taken from the kevin-keraudren translation of the MATLAB creaseg code.
# Requires the following python modules: glob, numpy, scipy, skimage, astropy,
# matplotlib
#
# Copyright 2019 Laura Boucheron
# This file is part of ACWE.
#
# ACWE is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# ACWE is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# ACWE.  If not, see <https://www.gnu.org/licenses/>.

import glob
import numpy as np
import scipy as sp
import skimage.transform
from astropy.io import fits
from matplotlib import pyplot as plt
import correct_limb_brightening
import acwe

def make_circle_mask(c,im_dims,r):
    # defomes a binary image with image dimensions im_dims of a circle with 
    # center c and radius r
    cx = c[0]
    cy = c[1]
    ix = im_dims[0]
    iy = im_dims[1]
    x,y = np.meshgrid(np.arange(-(cx),(ix-cx),1),np.arange(-(cy),(iy-cy),1))
    c_mask = (x**2+y**2)<=r**2
    return c_mask

plt.ion() # interactive plotting on 

# define the base directory under which fits files exist
# the code will loop over all fits files in the directory and segment CHs 
# from each
base_directory = './'
files = sorted(glob.glob(base_directory+'*.fits'))

# default parameters from [1]
resize_param = 8
foreground_weight = 1
background_weight = 1/50.
alpha = 0.3
narrowband = 2

for filename in files:
    print(filename)
    # Read in fits header and fits image
    hdulist = fits.open(filename)
    hdulist.verify('silentfix') #no clue why this is needed for successful data read
    h = hdulist[1].header
    J = hdulist[1].data
    if resize_param>1:
        #I = skimage.transform.resize(J,np.asarray(J.shape)/resize_param,\
        #                     order=3,preserve_range=True,anti_aliasing=True)
        I = skimage.transform.resize(J,np.asarray(J.shape)/resize_param,\
                             order=3,preserve_range=True,anti_aliasing=True)
    else:
        I = J

    # Determine characteristics of image
    im_size = np.asarray(I.shape) # size of the image
    sun_radius = h['R_SUN']/resize_param # solar radius from metadata

    # Correct limb brightening per Verbeeck et al. 2014
    I = correct_limb_brightening.correct_limb_brightening(I,sun_radius)

    # Define solar disk mask
    sd_mask = make_circle_mask(im_size/2,im_size,sun_radius)
 
    # Determine threshold value for initialization of AC as percentage of QS;  
    # estimage QS as maximum bin of histogram
    Ihist,ihist_edges = np.histogram(I[sd_mask],100) # 100 bin histogram of SD
    ihist_centers = (ihist_edges[:-1]+ihist_edges[1:])/2. # bin centers
    ihistmax = Ihist.argmax() # index of maximum bin
    initial_thresh = alpha*ihist_centers[ihistmax] # bin center of maximum bin
    m = (I<=initial_thresh)*sd_mask # initial mask
    m = sp.ndimage.morphology.binary_fill_holes(m) # fill holes
    
    # Set up variables for ACWE iterations
    m_seg = m # image to keep track of current initialization of ACWE
    I_seg = I # image of current segmentation
    I_seg[~sd_mask] = I[~m_seg&sd_mask].mean() # set pixels outside SD to mean 
                                               # of background to force ACWE to
                                               # ignore
    counter = 0 # to keep track of proxy of iterations
    seg_diff_cum = np.zeros(im_size) # to keep track of how many times pixels
                                     # change classes over iterations
    iterate = 1 # flag to continue iterating
    
    # continue iterating with 10 iterations of the ACWE evolution, followed by
    # check for convergence
    # running ACWE for 10 iterations is a good tradeoff between checking too
    # often and not often enough for convergence
    print('% Diff            % New Diff')
    while iterate:
        seg = acwe.acwe(I_seg,m_seg,10,(0,0,foreground_weight,\
                        background_weight),narrowband,1) # evolve ACWE for 10
                                                         # iterations
        
        # update current segmentation
        I_seg = I 
        I_seg[~sd_mask] = I[~seg&sd_mask].mean() 
        
        # difference in seg from previous iteration to now
        seg_diff = seg.astype(int) - m_seg.astype(int) 
        
        m_seg = seg # update m_seg image
        counter = counter + 1 # iterate counter
        
        # compute percentage of pixels that changed between previous iteration
        # and now
        percent_diff = float(abs(seg_diff).sum())/float(seg.sum())*100.0
        # keep track of how many times pixels have changed classes
        seg_diff_cum = seg_diff_cum + abs(seg_diff)
        # percentage of currently new pixels that have never changed classes
        # before
        percent_new_diff = float((seg_diff_cum*(abs(seg_diff)==1)).sum())/\
                           float((seg_diff_cum*(abs(seg_diff))>=1).sum()+\
                           np.finfo(float).eps)*100  
        print(str(percent_diff) + ' ' + str(percent_new_diff))
        if percent_new_diff==0 | ~(seg.sum()>0):
            iterate = 0
           
    hdulist.close()