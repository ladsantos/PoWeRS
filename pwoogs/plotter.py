#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

"""
This code is used to plot the results produced by moog.py
"""

class line(object):

    def __init__(self,params):

        # Modelling parameters
        self.c = 2.998E18
        self.wlshift = params[9,1]
        self.vshift = params[9,0]
        self.yadd = params[10,0]
        self.ymult = params[10,1]
        self.wl_start = params[0,1]
        self.wl_end = params[1,1]

        # Getting the numbers for the plot
        self.data = np.loadtxt('spectrum.dat')
        self.model = np.loadtxt('vm_smooth.out', skiprows=2)
        self.data[:,0] = self.data[:,0] + self.wlshift - self.data[:,0]*(
            self.c/(self.vshift*1E13+self.c)-1.0)
        self.data[:,1] = self.data[:,1] + self.yadd - \
            self.data[:,1]*self.ymult

    # Function that finds the index of a target value inside an array arr
    def find_index(self,target,arr):
        
        self.diff = np.abs(arr-target)
        self.index = np.where(self.diff == min(self.diff))[0][0]
        return self.index

    # Function that sets the x-axis limits of the analysis
    def x_set_limits(self,start,end,data2d):
        
        self.start_index = self.find_index(start,data2d[:,0])
        self.end_index = self.find_index(end,data2d[:,0])
        return data2d[self.start_index:self.end_index]

    # The plotting function
    def plot(self):
        
        # Setting the limits on the x-axis
        self.data_target = self.x_set_limits(self.wl_start,self.wl_end,self.data)
        self.model_target = self.x_set_limits(self.wl_start,self.wl_end,self.model)
        
        # And finally
        fig, ax = plt.subplots()
        ax.plot(self.data_target[:,0],self.data_target[:,1],'.')
        ax.plot(self.model_target[:,0],self.model_target[:,1])
        ax.ticklabel_format(useOffset=False)
        plt.xlabel(r'$\lambda$ ($\AA$)')
        plt.ylabel(r'$I_\lambda$')
        plt.show()