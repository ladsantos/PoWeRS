#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
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
    # v_macro is the macroturbulence velocity
    # Z is the atomic number of the chemical element that produces the line
    def __init__(self,spec_window,gauss,v_macro,line_file,line,**kwargs):
        
        # Default optional parameters:
        
        if ('x_vel' in kwargs):
            self.vshift = kwargs['x_vel']
        else:
            self.vshift = 0.0
            
        if ('x_wl' in kwargs):
            self.xshift = kwargs['x_wl']
        else:
            self.xshift = 0.0
        
        if ('y_add' in kwargs):
            self.yadd = kwargs['y_add']
        else:
            self.yadd = 0.0
            
        if ('y_mult' in kwargs):
            self.ymult = kwargs['y_mult']
        else:
            self.ymult = 0.0
        
        if ('gamma' in kwargs):
            self.gamma = kwargs['gamma']
        else:
            self.gamma = 0.0
        
        if ('perf_radius' in kwargs):
            self.radius = kwargs['perf_radius']
        else:
            self.radius = 5
        
        self.c = 2.998E18
        self.am = utils.arr_manage()
        self.n = normal.optimize()
        self.spec = spec_window
        self.gauss = gauss
        self.v_m = v_macro
        self.lines = np.loadtxt(line_file,skiprows=1,usecols=(0,1))
        self.Z = self.lines[line,1]
        self.line_center = self.lines[line,0]
        
    # Function that writes to params.txt
    def write_params(self,log_rot_v,abund):
        
        self.rot_v = 10**log_rot_v
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
smoothing       %.3f            0.0             1 = Gaussian, 2 = Lorentzian
darkening       0.6             None            Limb darkening coefficient of a rotational broadening function
velocities      %.2f            %.1f            1 = macroturbulence, 2 = rotation
nabunds         1               None            # of elements
abund           %i              %.3f'''
                % (self.spec[0],self.spec[0],self.spec[1],self.spec[1],
                   self.vshift,self.xshift,self.yadd,self.ymult,
                   self.gauss,self.v_m,self.rot_v,
                   self.Z,abund)
                )
    
    
    '''
    The performance function: first it creates the params.txt file, then runs
    moog in silent mode, interpolates the generated model to the points of
    the observed spectrum, and then simply calculates the sum of squared 
    differences, weighted by the inverse of the observed spectrum to the power
    of alpha.
    '''
    def perf(self,p):
        
        # Applying corrections to the observed spectrum
        self.data = np.loadtxt('spectrum.dat')
        #self.line = np.loadtxt(self.line_file)
        self.data[:,0] = self.data[:,0] + self.xshift - self.data[:,0]*(
            self.c/(self.vshift*1E13+self.c)-1.0)
        self.data[:,1] = self.data[:,1]*self.ymult + self.yadd
        self.data_target = self.am.x_set_limits(self.spec[0],self.spec[1],self.data)
        
        # Running MOOGSILENT
        self.write_params(np.log10(p[0]),p[1])
        m = moog.run(silent=True)
        
        # Evaluating the performance in a radius around the center of the line
        self.center_index = self.am.find_index(self.line_center,self.data_target[:,0])
        self.ci0 = self.center_index-self.radius
        self.ci1 = self.center_index+self.radius+1
        self.model = np.loadtxt('vm_smooth.out', skiprows=2)
        self.model_interp = np.interp(self.data_target[self.ci0:self.ci1,0],
                                      self.model[:,0],
                                      self.model[:,1])
        #plt.plot(self.data_target[self.ci0:self.ci1,0],self.data_target[self.ci0:self.ci1,1])
        #plt.plot(self.data_target[self.ci0:self.ci1,0],self.model_interp[:])
        #plt.show()
        
        return np.sum((1.0/self.data_target[self.ci0:self.ci1,1])**self.gamma\
            *(self.data_target[self.ci0:self.ci1,1]-\
                self.model_interp[:])**2)

"""    
    def find(self,guess_min,guess_max,**kwargs):
        
        if ('N' in kwargs):
            self.N = kwargs['N']
        else:
            self.N = 100
        
        if ('alpha' in kwargs):
            self.alpha = kwargs['alpha']
        else:
            self.alpha = 0.5

        if ('beta' in kwargs):
            self.beta = kwargs['beta']
        else:
            self.beta = 0.1
        
        if ('c_limit' in kwargs):
            self.c_limit = kwargs['c_limit']
        else:
            self.c_limit = 1E-2
            
        if ('s_limit' in kwargs):
            self.s_limit = kwargs['s_limit']
        else:
            self.s_limit = 1E-1
            
        if ('rho' in kwargs):
            self.rho = kwargs['rho']
        else:
            self.rho = 0.1
            
        if ('plot' in kwargs):
            self.plot = kwargs['plot']
        else:
            self.plot = True
        
        self.p_min = guess_min
        self.p_max = guess_max
        self.p_mean = (self.p_min+self.p_max)/2.
        self.p_sigma = (self.p_max-self.p_min)/2.
        
        # The estimation by CREPE is done in just one line:
        self.new_p_mean,self.new_p_sigma = self.n.estimate(self.perf,\
            self.p_mean,self.p_sigma,c_limit=self.c_limit,verbose=True,\
                alpha=self.alpha,beta=self.beta,N=self.N,rho=self.rho)

        # Printing and plotting the results
        print 'log_v_rot = %.3f p/m %.3f' % (self.new_p_mean[0],\
            self.new_p_sigma[0])
        #print 'xshift = %.4f p/m %.4f' % (self.new_p_mean[2],\
        #    self.new_p_sigma[2])
        print 'abund = %.3f p/m %.3f' % (self.new_p_mean[1],\
            self.new_p_sigma[1])

        self.write_params(self.new_p_mean[0],self.new_p_mean[1])
        
        if self.plot == True:
            m = moog.run()
        else:
            m = moog.run(silent=True)
"""