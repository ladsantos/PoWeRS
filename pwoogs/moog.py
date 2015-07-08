#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
import numpy as np
import plotter

class driver(object):
    
    def __init__(self,params):
        
        # Output parameters
        self.standard_out = 'vm_long.out'
        self.summary_out = 'vm_li.out'
        self.smoothed_out = 'vm_smooth.out'
        
        # Input parameters
        self.model_in = 'star.mod'
        self.lines_in = 'sun_syn.lin'
        self.observed_in = 'spectrum.dat'
        
        # Modelling parameters
        self.syn_start = params[0]
        self.syn_end = params[1]
        self.wl_start = params[2]
        self.wl_end = params[3]
        self.xshift = params[4]
        self.abund = params[5]
        self.microv = params[6]
        self.rotv = params[7]
    
    # Writes the MOOG driver file batch.par
    def create_batch(self):
        
        with open('batch.par','w') as f:
            f.truncate()
            f.write('synth\n')
            f.write('standard_out  %s\n' % self.standard_out)
            f.write('summary_out  %s\n' % self.summary_out)
            f.write('smoothed_out  %s\n' % self.smoothed_out)
            f.write('model_in  %s\n' % self.model_in)
            f.write('lines_in  %s\n' % self.lines_in)
            f.write('observed_in  %s\n' % self.observed_in)
            f.write('atmosphere  1\n')
            f.write('molecules  1\n')
            f.write('lines    1\n')
            f.write('flux/int  0\n')
            f.write('damping    0\n')
            f.write('freeform  1\n')
            f.write('plot    3\n')
            f.write('abundances  14 1')
            f.write('''
    3 +0.0
    6 +0.0
    7 +0.0
   14 +0.0
   20 +0.0
   21 +0.0
   22 +0.0
   23 +0.0
   24 +0.0
   25 +0.0
   26 %.3f
   44 +0.0
   58 +0.0
   62 +0.0\n''' 
                % self.abund)
            f.write('isotopes   0  1\n')
            f.write('synlimits\n')
            f.write(' %.1f %.1f 0.01 2.0\n' % (self.syn_start,self.syn_end))
            f.write('obspectrum  5\n')
            f.write('plotpars  1\n')
            f.write(' %.2f %.2f 0.5 1.05\n' % (self.wl_start,self.wl_end))
            f.write(' 0.0  %.4f  0.001  1.000\n' % self.xshift)
            f.write(' r  0.052  %.1f  0.6  %.2f  0.0' % (self.rotv,self.microv))

class run(object):
    
    def __init__(self):
        self.params = np.genfromtxt('params.txt',usecols=(1))            
        self.d = driver(self.params).create_batch()
        os.system('MOOGSILENT > moog.log 2>&1')

        # Plotting
        self.p = plotter.line(self.params).plot()