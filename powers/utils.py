#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


class arr_manage(object):

    
    # This routine finds the index of a value closest to the target in a 
    # numpy-array
    def find_index(self,target,array):
        
        self.diff = np.abs(array-target)
        self.index = np.where(self.diff == min(self.diff))[0][0]
        return self.index
    
    
    # This routine returns a section of an array given the start and end values
    def x_set_limits(self,start,end,data2d):
        
        assert start > data2d[0,0], 'Invalid start value'
        assert end < data2d[-1,0], 'Invalid end value'
        self.start_index = self.find_index(start,data2d[:,0])
        self.end_index = self.find_index(end,data2d[:,0])
        return data2d[self.start_index:self.end_index]
    
    
    # This routine writes the section of a datafile to a new one, for which 
    # the start and end values are chosen for a specific column of the 
    # datafile. By default, the column is 0
    def cut(self,start,end,file_full,file_target,**kwargs):
        
        # Setting default column to look for the start and end
        if ('col' in kwargs):
            self.col = kwargs['col']
        else:
            self.col = 0
        
        # ISSUE: the code works only with 2-column datafiles as of now
        self.data = np.loadtxt(file_full)
        self.start_index = self.find_index(start,self.data[:,self.col])
        self.end_index = self.find_index(end,self.data[:,self.col])
        self.N = self.end_index-self.start_index
        self.M = len(self.data[0,:])
        with open(file_target,'w') as f:
            f.truncate()
            for i in range(self.N):
                f.write('%.7f\t%.7f\n' % \
                    (
                    self.data[i+self.start_index,0],
                    self.data[i+self.start_index,1]
                    )
                )
        
        
    # This routine finds the center of a line and returns a wavelength linear 
    # shift in order to correct for the error in centralization
    def find_center(self,data):
        
        p = np.polyfit(x=data[:,0],y=data[:,1],deg=2)
        N = len(data[:,0])
        xn = np.linspace(data[0,0],data[N-1,0])
        yn = np.array([p[0]*xk**2+p[1]*xk+p[2] for xk in xn])
        return -p[1]/2./p[0]
    
    
    # This routine finds the multiplicative factor to correct the normalization
    # of a region of the spectrum, using the highest region (radius) of the 
    # spectrum inside a section contained on data
    def find_corr(self,data,radius,**kwargs):
        
        # Do you want the program to be silent?
        if ('silent' in kwargs):
            self.silent = kwargs['silent']
        else:
            self.silent = False
        
        self.m = max(data[:,1])
        self.ind = self.find_index(self.m,data[:,1])
        self.stdev = np.std(data[self.ind-radius:self.ind+radius+1,1])
        if self.silent == False:
            print 'Standard deviation of the spectrum on the ' \
                  'selected data points = %.4f' % self.stdev
        return np.mean(data[self.ind-radius:self.ind+radius+1,1])
    
    
    # This routine finds the multiplicative factor to correct the normalization
    # of a region of the spectrum, using an ensemble of user defined points, 
    # based on which the correction will be estimated around a region within
    # a radius of points
    def find_corr_from_ensemble(self,data,target_wls,radius):
        
        self.corrs = np.empty(len(target_wls),float)
        for i in range(len(target_wls)):
            self.ind = self.find_index(target_wls[i],data[:,1])
            self.corrs[i] = np.mean(data[self.ind-radius:self.ind+radius+1,1])
        return self.corrs