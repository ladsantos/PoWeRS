#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from pwoogs import plotter,moog,utils
from crepe import normal

"""
This code is used to estimate the projected rotation speed of a star given a
spectral line, correction factors and an initial guess. It uses the 
cross-entropy (CE) method to do estimation. You need to have all the necessary
input files for MOOG (atmosphere model, observed spectrum, lines file and a
batch.par file, the latter can be empty). For now, it works only for single-line
analysis.
"""

class rotation(object):
    
    # spec_window is the spectral analysis window, a 1x2 numpy array
    # x_vel is the velocity shift to be applied to the x-axis
    # x_wl is the wavelengths shift to be applied to the x-axis
    # y_add is the additive shift to be applied to the spectrum
    # y_mult is the multiplicative shift to be applied to the spectrum
    # par is the a 2x2 numpy array containing:
    # [[ veloc guess min, veloc guess max ],
    #  [ abund guess min, abund guess max ]]
    # macro_v is the macroturbulence velocity
    # Z is the atomic number of the chemical element that produces the line
    def __init__(self,spec_window,x_vel,x_wl,y_add,y_mult,macro_v,Z):
        
        self.c = 2.998E18
        self.am = utils.arr_manage()
        self.n = normal.optimize()
        self.spec = spec_window
        self.vshift = x_vel
        self.wlshift = x_wl
        self.yadd = y_add
        self.ymult = y_mult
        self.m_v = macro_v
        self.Z = Z
        self.data = np.loadtxt('spectrum.dat')
        self.data[:,0] = self.data[:,0] + self.wlshift - self.data[:,0]*(
            self.c/(self.vshift*1E13+self.c)-1.0)
        self.data[:,1] = self.data[:,1] + self.yadd - \
            self.data[:,1]*self.ymult
        self.data_target = self.am.x_set_limits(self.spec[0],self.spec[1],self.data)
        
    # Function that writes to params.txt
    def write_params(self,log_rot_v,abund):
        
        self.rot_v = 10**log_rot_v
        #print self.rot_v
        with open('params.txt','w') as f:
            f.truncate()
            f.write(
                '''Parameter       Value 1         Value 2         Comment
l_limit         %.1f            %.1f            Lower limits: 1 = synthesis, 2 = plotting, angstrons
u_limit         %.1f            %.1f            Upper limits: 1 = synthesis, 2 = plotting, angstrons
synth_pars      0.01            2.0             1 = Step size of the synthesis, 2 = wavelength point to consider opacity contributions from neighboring transitions
atmosphere      1               None            See description on WRITEMOOG.ps
molecules       1               None            See description on WRITEMOOG.ps
trudamp         1               None            See description on WRITEMOOG.ps
lines           1               None            See description on WRITEMOOG.ps
flux/int        0               None            See description on WRITEMOOG.ps
damping         0               None            See description on WRITEMOOG.ps
xshift          %.1f            %.4f            Shift on the x-axis, 1 = velocity, 2 = wavelength, for the observed spectrum
yshift          %.4f            %.4f            Shift on the y-axis, 1 = additive, 2 = multiplicative, for the observed spectrum
smoothing       0.052           0.0             1 = Gaussian, 2 = Lorentzian
darkening       0.6             None            Limb darkening coefficient of a rotational broadening function
velocities      %.2f            %.3f            1 = macroturbulence, 2 = rotation
nabunds         1               None            # of elements
abund           %i              %.3f'''
                % (self.spec[0],self.spec[0],self.spec[1],self.spec[1],
                   self.vshift,self.wlshift,self.yadd,self.ymult,
                   self.m_v,self.rot_v,
                   self.Z,abund)
                )
                
    '''
    The performance function: first it creates the params.txt file, then runs
    moog in silent mode, interpolates the generated model to the points of
    the observed spectrum, and then simply calculates the sum of squared 
    differences, weighted by the inverse of the observed spectrum to the power
    of 0.1. In the future, I'll implement a way to change this power.
    '''
    def perf(self,p):
        
        foo = self.write_params(p[0],p[1])
        m = moog.run(silent=True)
        self.model = np.loadtxt('vm_smooth.out', skiprows=2)
        self.model_interp = np.interp(self.data_target[:,0],
                                      self.model[:,0],self.model[:,1])
        return np.sum((1.0/self.data_target[:,1])**0.1*(self.data_target[:,1]-\
            self.model_interp[:])**2)
    
    def find(self,guess_min,guess_max):
        
        self.p_min = guess_min
        self.p_max = guess_max
        self.p_mean = (self.p_min+self.p_max)/2.
        self.p_sigma = (self.p_max-self.p_min)/2.
        self.n = normal.optimize()
        
        # The estimation by CREPE is done in just one line:
        self.new_p_mean,self.new_p_sigma = self.n.estimate(self.perf,\
            self.p_mean,self.p_sigma,c_limit=1E-2,verbose=True,alpha=0.1)

        # Printing and plotting the results
        print 'log_v_rot = %.3f p/m %.3f' % (self.new_p_mean[0],\
            self.new_p_sigma[0])
        print 'abund = %.3f p/m %.3f' % (self.new_p_mean[1],\
            self.new_p_sigma[1])

        self.write_params(self.new_p_mean[0],self.new_p_mean[1])
        m = moog.run()