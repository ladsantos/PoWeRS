#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from pwoogs import moog,estimate,utils
import matplotlib.pyplot as plt
from sys import argv

# WARNING: this code is not polished yet and is horrible to read

#script,choice,vsini,abund = argv
script,choice = argv
choice = int(choice)

# Synthesis parameters
line_file = 'sun_syn.lin'
lines = np.loadtxt(line_file,skiprows=1,usecols=(0,1))
interval = 1.0
res_power = 65000.

# Star parameters
star_info = np.genfromtxt('star.mod',skip_header=1,skip_footer=83,usecols=(0,1),delimiter='/  ')
T_star = star_info[0]

v_macro = 3.6+(T_star-5777.)/486.

# Other stuff
u = utils.arr_manage()

print 'Managing the data file.'
#choice = 0
spec_window = np.array([lines[choice,0]-interval/2,lines[choice,0]+interval/2])
u.cut(spec_window[0]-10.,spec_window[1]+10.,'spectrum_full.dat','spectrum.dat')
print 'Done.\n'

# The instrumental broadening
gauss = np.mean(spec_window)/res_power

print 'Finding the shift on the wavelength.'
radius_1 = 2    # number of points
data = np.loadtxt('spectrum.dat')
ind = u.find_index(lines[choice,0],data[:,0])
center = u.find_center(data[ind-radius_1+1:ind+radius_1+2])
wl_shift = lines[choice,0]-center
print 'Wavelength shift = %.4f\n' % wl_shift

print "Finding the correction factor for the continuum."
radius_2 = 3.0  # radius of the spectrum to analyze, angstroms
radius_3 = 3    # radius of points around the maximum of the spectrum
ind_min = u.find_index(lines[choice,0]-radius_2,data[:,0])
ind_max = u.find_index(lines[choice,0]+radius_2,data[:,0])
target_wls = np.loadtxt('continuum.dat')
corr = 1.0/np.mean(u.find_corr_from_ensemble(
    data[ind_min:ind_max,:],
    target_wls[choice,:],
    radius_3
    ))
print "Correction factor = %.4f" % corr

# Instatiating the function to write parameters for MOOG
r = estimate.rotation(
            spec_window,
            gauss,
            v_macro,
            'sun_syn.lin',
            choice,
            x_wl=wl_shift,
            y_mult=corr,
            bwing_w = 10.0,
            center_w = 20.0,
            perf_radius=8
            )

# Finding vsini and abundance
r.find(N=15,I=10)