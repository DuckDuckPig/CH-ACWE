#-------------------------------------------------------------------------------
# Correct Limb Brightening
#
# I_smooth = correct_limb_brightening(I,sun_radius)
#
# Implements the limb brighetning correction of Verbeeck et al. 2014 [1].
#
# Inputs:
#     I: input image
#     sun_radius: radius of the sun in image I in pixels
#
# Output:
#     I_smooth: corrected image
#
# References:
# [1] C. Verbeeck, V. Delouille, B. Mampaey, & R. De Visscher, "The SPoCA-
#     suite: Software for extraction, characterization and tracking of active
#     regions and coronal holes on EUV images," Astronomy & Astrophysics, vol.
#     561, pp. A29, 2014.
#
# Notes:
# Requires the following python modules: numpy
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
# You should have received a copy of the GNU General Public License along
# with ACWE.  If not, see <https://www.gnu.org/licenses/>.


# Edited By Jeremy A. Grajeda, Apirl 2, 2021

import numpy as np

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

def correct_limb_brightening(I,sun_center,sun_radius):
    im_size = np.asarray(I.shape)
    # make solar disk masks for the different regions of correction per [1]
    sd_mask = make_circle_mask(sun_center,im_size,sun_radius) 
    r1, r2, r3, r4 = 0.7, 0.95, 1.08, 1.12
    sd_mask1 = make_circle_mask(sun_center,im_size,sun_radius*r1)
    sd_mask2 = make_circle_mask(sun_center,im_size,sun_radius*r2)
    sd_mask3 = make_circle_mask(sun_center,im_size,sun_radius*r3)
    sd_mask4 = make_circle_mask(sun_center,im_size,sun_radius*r4)
    
    # compute average intensity within each annulus of 1 pixel wide
    F = np.zeros(im_size)
    for r in np.arange(r1*sun_radius,r4*sun_radius,1):
        annulus1 = make_circle_mask(sun_center,im_size,r)
        annulus2 = make_circle_mask(sun_center,im_size,r+1)
        annulus = (annulus2^annulus1)>0
        F[annulus] = (annulus*I).sum()/annulus.sum()
    # define corrected image per [1]
    I_corr = np.zeros(im_size)
    I_corr[F>0] = np.median(I[sd_mask])*I[F>0]/F[F>0]

    # define smoothed corrected image
    I_smooth = np.zeros(im_size)
    # no correction for r<r1 or r>r4
    I_smooth[sd_mask1] = I[sd_mask1]
    I_smooth[~sd_mask4] = I[~sd_mask4]
 
    # complete correction for r2<r<r3
    region = (sd_mask3^sd_mask2)>0
    I_smooth[region] = I_corr[region]

    # smoothed correction for r1<r<r2
    for r in np.arange(r1*sun_radius,r2*sun_radius,1):
        annulus1 = make_circle_mask(sun_center,im_size,r)
        annulus2 = make_circle_mask(sun_center,im_size,r+1)
        annulus = (annulus2^annulus1)>0
        f = 0.5*np.sin(np.pi/(r2-r1)*(r/sun_radius-(r1+r2)/2))+0.5
        I_smooth[annulus] = (1-f)*I[annulus] + f*I_corr[annulus]
  
    # smoothed correction for r3<r<r4
    for r in np.arange(r3*sun_radius,r4*sun_radius,1):
        annulus1 = make_circle_mask(sun_center,im_size,r)
        annulus2 = make_circle_mask(sun_center,im_size,r+1)
        annulus = (annulus2^annulus1)>0
        f = 0.5*np.sin(np.pi/(r4-r3)*(r/sun_radius+(r4-3*r3)/2))+0.5
        I_smooth[annulus] = (1-f)*I[annulus] + f*I_corr[annulus]
    return I_smooth
