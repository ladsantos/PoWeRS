#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

"""
This code is used to plot the results produced by moog.py
"""


class line(object):


    def __init__(self, params):

        # Modelling parameters
        self.c = 2.998E18
        self.wlshift = params[9,1]
        self.vshift = params[9,0]
        self.yadd = params[10,0]
        self.ymult = params[10,1]
        self.wl_start = params[0,1]
        self.wl_end = params[1,1]
        self.v_rot = params[13,1]
        self.v_macro = params[13,0]
        self.fwhm_gauss = params[11,0]
        self.abund = params[15,1]

        # Getting the numbers for the plot
        self.data = np.loadtxt('spectrum.dat')
        self.model = np.loadtxt('vm_smooth.out', skiprows=2)
        self.data[:,0] = self.data[:,0] + self.wlshift - self.data[:,0]* \
            (self.c / (self.vshift*1E13+self.c) - 1.0)
        self.data[:,1] = self.data[:,1] * self.ymult + self.yadd


    # Function that finds the index of a target value inside an array arr
    def find_index(self, target, arr):
        
        self.diff = np.abs(arr-target)
        self.index = np.where(self.diff == min(self.diff))[0][0]
        return self.index


    # Function that sets the x-axis limits of the analysis
    def x_set_limits(self, start, end, data2d):
        
        self.start_index = self.find_index(start, data2d[:,0])
        self.end_index = self.find_index(end, data2d[:,0])
        return data2d[self.start_index:self.end_index]


    # The plotting function
    def plot(self, **kwargs):
        
        # mode is set to plot on window by default. Any other value will
        # save a file with that value as name. For example: mode='plot.eps'
        # will save the plot as plot.eps 
        if ('mode' in kwargs):
            self.mode = kwargs['mode']
        else:
            self.mode = 'window'
            
        # Getting star name or using the default
        if ('star_name' in kwargs):
            self.name = kwargs['star_name']
        else:
            self.name = 'Star'
        
        # Setting the limits on the x-axis
        self.data_target = self.x_set_limits(self.wl_start, self.wl_end, 
                                             self.data)
        self.model_target = self.x_set_limits(self.wl_start, self.wl_end, 
                                              self.model)
        
        # And finally
        fig, ax = plt.subplots()
        ax.plot(self.data_target[:,0], self.data_target[:,1], '.')
        ax.plot(self.model_target[:,0], self.model_target[:,1])
        ax.ticklabel_format(useOffset=False)
        plt.title(
            r'%s, $v \sin{i} = %.2f$, $v_{macro} = %.2f$, ' \
                r'$FWHM_{gauss} = %.3f$, $abund = %.3f$' \
                % (self.name, self.v_rot, self.v_macro, self.fwhm_gauss, \
                   self.abund)
            )
        plt.xlabel(r'$\lambda$ ($\AA$)')
        plt.ylabel(r'$F_\lambda$ $d\lambda$')
        
        if self.mode == 'window':
            plt.show()
        else:
            print 'Saving plot as %s' % self.mode
            plt.savefig(self.mode)