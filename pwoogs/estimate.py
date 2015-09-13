#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import time
from pwoogs import plotter,moog,utils

"""
This code is used to estimate the projected rotation speed of a star given a
spectral line, correction factors and an initial guess. You need to have all 
the necessary input files for MOOG (atmosphere model, observed spectrum, and 
a lines file). For now, it works only for single-line analysis.
"""


class vsini(object):
  
  
    # spec_window is the spectral analysis window, a 1x2 numpy array
    # gauss is the instrumental broadening parameter
    # v_macro is the macroturbulence velocity
    # line_file is the name of the file containing the chosen lines
    # line is which of the lines on the previous file to work on
    # SN is the signal-to-noise ratio of the spectrum
    def __init__(self, spec_window, gauss, v_macro, line_file, line, SN, 
                 **kwargs):
        
        # Default optional parameters:
        
        # x_vel is the velocity shift to be applied to the x-axis
        if ('x_vel' in kwargs):
            self.vshift = kwargs['x_vel']
        else:
            self.vshift = 0.0
        
        # x_wl is the wavelengths shift to be applied to the x-axis
        if ('x_wl' in kwargs):
            self.xshift = kwargs['x_wl']
        else:
            self.xshift = 0.0
        
        # y_add is the additive shift to be applied to the spectrum
        if ('y_add' in kwargs):
            self.yadd = kwargs['y_add']
        else:
            self.yadd = 0.0
        
        # y_mult is the multiplicative shift to be applied to the spectrum
        if ('y_mult' in kwargs):
            self.ymult = kwargs['y_mult']
        else:
            self.ymult = 0.0
        
        # perf_radius is the number of points around the line center where
        # to evaluate the performance of the synthetic spectrum
        if ('perf_radius' in kwargs):
            self.radius = kwargs['perf_radius']
        else:
            self.radius = 5
        
        # bwing_w is the weight to be applied to the blue side of the line
        # when evaluating the performance
        if ('bwing_w' in kwargs):
            self.bwing_w = kwargs['bwing_w']
        else:
            self.bwing_w = 5.0
        
        # bwing_w is the weight to be applied to the red side of the line
        # when evaluating the performance
        if ('rwing_w' in kwargs):
            self.rwing_w = kwargs['rwing_w']
        else:
            self.rwing_w = 1.0
        
        # center_w is the weight to be applied to the line center when 
        # evaluating the performance
        if ('center_w' in kwargs):
            self.center_w = kwargs['center_w']
        else:
            self.center_w = 10.0
        
        # Maximum number of points around the performance radius that are
        # allowed to be a bad fit (1 S/N sigma lower than observed signal)
        # If this limit is exceeded, the variable badfit_status will return
        # True after running find()
        # For high precision spectrum, set this to a very low number
        if ('badfit_tol' in kwargs):
            self.badfit_tol = kwargs['badfit_tol']
        else:
            self.badfit_tol = 10
        
        self.c = 2.998E18
        self.am = utils.arr_manage()
        self.spec = spec_window
        self.gauss = gauss
        self.v_m = v_macro
        self.lines = np.loadtxt(line_file, skiprows=1, usecols=(0,1))
        self.Z = self.lines[line,1]
        self.line_center = self.lines[line,0]
        self.spec_sigma = 1./SN
        
        
    # Function that writes to params.txt
    # WARNING: this needs to be written in a more readable fashion
    def write_params(self, vsini, abund):
        
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
velocities      %.2f            %.3f            1 = macroturbulence, 2 = v sin i
nabunds         1               None            # of elements
abund           %i              %.5f'''
                % (self.spec[0], self.spec[0], self.spec[1], self.spec[1],
                   self.vshift, self.xshift, self.yadd, self.ymult,
                   self.gauss, self.v_m,vsini, self.Z,abund)
                )
    
    
    '''
    The performance function: first it creates the params.txt file, then runs
    moog in silent mode, interpolates the generated model to the points of
    the observed spectrum, and then simply calculates the sum of squared 
    differences, weighted by the inverse of the observed spectrum to the power
    of alpha.
    '''
    def perf(self, p):
        
        # Applying corrections to the observed spectrum
        self.data = np.loadtxt('spectrum.dat')
        self.data[:,0] = self.data[:,0] + self.xshift - self.data[:,0] * \
            (self.c / (self.vshift*1E13 + self.c) - 1.0)
        self.data[:,1] = self.data[:,1] * self.ymult + self.yadd
        self.data_target = self.am.x_set_limits(self.spec[0], self.spec[1], 
                                                self.data)
        
        # Running MOOGSILENT
        #self.write_params(np.log10(p[0]),p[1])
        self.write_params(p[0], p[1])
        m = moog.run(silent=True)
        
        # Evaluating the performance in a radius around the center of the line
        self.center_index = self.am.find_index(self.line_center, 
                                               self.data_target[:,0])
        self.ci0 = self.center_index - self.radius
        self.ci1 = self.center_index + self.radius+1
        self.model = np.loadtxt('vm_smooth.out', skiprows=2)
        self.model_interp = np.interp(self.data_target[self.ci0:self.ci1,0],
                                      self.model[:,0],
                                      self.model[:,1])
        
        # Checking the fit on line wings
        self.check = self.data_target[self.ci0:self.ci1,1] - self.model_interp
        self.check = len(np.where(self.check > self.spec_sigma)[0])
        
        # Creating the weights vector
        self.w = np.zeros(2 * self.radius + 1, float)
        for i in range(self.radius):
            self.w[i] = self.bwing_w
            self.w[i + self.radius+1] = self.rwing_w
        self.w[self.radius] = self.center_w
        
        return np.sum(self.w[:] * (self.data_target[self.ci0:self.ci1,1] - \
                self.model_interp[:])**2) / np.sum(self.w)


    # WARNING: the following is not PEP8 compliant yet
    def find(self, **kwargs):
        
        # Number of points to try for each iteration
        if ('N' in kwargs):
            self.pts = kwargs['N']
        else:
            self.pts = 15
        
        # Narrowing factor when going to the next iteration
        if ('pace' in kwargs):
            self.pace = kwargs['pace']
        else:
            self.pace = 2.0
        
        # Initial guess range for abundance. It has to be a numpy array of 
        # length = 2
        if ('a_guess' in kwargs):
            self.a_guess = kwargs['a_guess']
        else:
            self.a_guess = np.array([-0.100,0.100])
        
        # Initial guess range for vsini. It has to be a numpy array of 
        # length = 2
        if ('v_guess' in kwargs):
            self.v_guess = kwargs['v_guess']
        else:
            self.v_guess = np.array([0.5,10.0])
        
        # Maximum number of iterations
        if ('max' in kwargs):
            self.max_i = kwargs['I']
        else:
            self.max_i = 20
        
        # Convergence limits: a numpy array with length 2, corresponding to the
        # limits of vsini and abundance, respectively
        if ('limits' in kwargs):
            self.limits = kwargs['limits']
        else:
            self.limits = np.array([0.1,0.001])
            
        # Plot the spectral line fit at the end?
        if ('plot' in kwargs):
            self.plot = kwargs['plot']
        else:
            self.plot = True
        
        # Set 'save' to a filename with an extension (e.g. png, eps) 
        # Overrides 'plot' to False
        if ('save' in kwargs):
            self.save = kwargs['save']
            self.plot = False
        else:
            self.save = None
        
        # Do you want the program to be silent? Not very useful at the moment
        if ('silent' in kwargs):
            self.silent = kwargs['silent']
        else:
            self.silent = False
        
        if self.silent == False:
            print "\nStarting estimation of vsini and abundance...\n"
        self.t0 = time.time()
        self.best_a = np.mean(self.a_guess)
        self.best_v = np.mean(self.v_guess)
        self.it = 1
        self.finish = False
        self.badfit_status = False
        
        while self.finish == False and self.it < self.max_i:
            
            # Evaluating vsini
            self.v_grid = np.linspace(self.v_guess[0],self.v_guess[1],self.pts)
            self.S = np.array([self.perf(np.array([self.v_grid[k],\
                self.best_a])) for k in range(self.pts)])
            self.best_v_ind = np.where(self.S == min(self.S))[0][0]
            self.best_v = self.v_grid[self.best_v_ind]
            
            # Evaluating abundance
            self.a_grid = np.linspace(self.a_guess[0],self.a_guess[1],self.pts)
            self.S= np.array([self.perf(np.array([self.best_v,self.a_grid[k]]))\
                for k in range(self.pts)])
            self.best_a_ind = np.where(self.S == min(self.S))[0][0]
            self.best_a = self.a_grid[self.best_a_ind]
            
            self.go_v = True
            self.go_a = True
            
            # Checking if the best values are too near the edges of the guess
            if self.best_v_ind == 0 or self.best_v_ind == self.pts-1:
                self.go_v = False
            elif self.best_a_ind == 0 or self.best_a_ind == self.pts-1:
                self.go_a = False

            # Calculating changes
            self.v_change = np.abs(self.best_v-np.mean(self.v_guess))
            self.a_change = np.abs(self.best_a-np.mean(self.a_guess))
            if self.silent == False:
                if self.v_change < self.limits[0] and self.a_change < self.limits[1]:
                    self.finish = True
                    print "final iteration = %i" % self.it
                else:
                    print "iteration = %i" % self.it
                print "best vsini = %.1f (changed %.3f)" % \
                    (self.best_v,self.v_change)
                print "best abund = %.3f (changed %.4f)\n" % \
                    (self.best_a,self.a_change)
            
            # Setting the new guess. If one of the best values are too near the
            # edges of the previous guess, it will not narrow its new guess range
            self.v_width = self.v_guess[1]-self.v_guess[0]
            self.a_width = self.a_guess[1]-self.a_guess[0]
            if self.go_v == True:
                self.v_guess = np.array([self.best_v-self.v_width/2/self.pace,\
                    self.best_v+self.v_width/2/self.pace])
            else:
                self.v_guess = np.array([self.best_v-self.v_width/2,\
                    self.best_v+self.v_width/2])
            if self.go_a == True:
                self.a_guess = np.array([self.best_a-self.a_width/2/self.pace,\
                    self.best_a+self.a_width/2/self.pace])
            
            # Checking if the v_guess contains vsini lower than 0.5. If True,
            # it will add a value to the array so that the lower limit is 0.5
            if self.v_guess[0] < 0.5 and self.silent == False:
                print "WARNING: vsini guess is less than 0.5. I will shift it."
                self.v_guess += 0.5-self.v_guess[0]
            
            else:
                self.a_guess = np.array([self.best_a-self.a_width/2,\
                    self.best_a+self.a_width/2])
            
            self.it += 1
        
        # Finalizing the routine
        if self.it == self.max_i and self.silent == False:
            print "Maximum number of iterations reached!\n"
        self.t1 = time.time()
        self.total_time = self.t1-self.t0
        if self.silent == False:
            print "Estimation took %.3f seconds." % self.total_time
        self.write_params(self.best_v,self.best_a)

        if self.plot == True:
            m = moog.run(silent=False)
        elif (self.plot == False or self.silent == True) and self.save != None:
            m = moog.run(silent=False,save=self.save)
        else:
            m = moog.run(silent=True)
            
        # Trigger bad fit warning
        self.S = self.perf(np.array([self.best_v,self.best_a]))
        if self.check > self.badfit_tol:
            self.badfit_status = True
            if self.silent == False:
                print """
    It is possible that the estimation is bad. This is usually caused
    by a slow-rotating star, line-blending, line center shift or bad continuum. 
    I suggest using a lower pace and more points, doing it manually or 
    discarding result for this particular line.
                """
        
        return self.best_v,self.best_a,self.badfit_status