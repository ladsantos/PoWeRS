#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class arr_manage(object):
    
    # This routine finds the index of a specific value (target) in a numpy-array
    def find_index(self,target,arr):
        
        self.diff = np.abs(arr-target)
        self.index = np.where(self.diff == min(self.diff))[0][0]
        return self.index
    
    # This routine performs a "soft-pop" in a numpy array
    def x_set_limits(self,start,end,data2d):
        
        self.start_index = self.find_index(start,data2d[:,0])
        self.end_index = self.find_index(end,data2d[:,0])
        return data2d[self.start_index:self.end_index]
    
    # This routine performs a "hard-pop" from a data file
    def cut(self,start,end,file_full,file_target):
        
        self.data = np.loadtxt(file_full)
        self.start_index = self.find_index(start,self.data[:,0])
        self.end_index = self.find_index(end,self.data[:,0])
        self.N = self.end_index-self.start_index
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
        plt.plot(data[:,0],data[:,1],'.',label='Data')
        plt.plot(xn,yn,label='2nd order poly')
        plt.xlabel('Wavelength (angstroms)')
        plt.ylabel('Flux')
        plt.legend()
        plt.show()
        return -p[1]/2./p[0]
    
    # This routine finds the multiplicative factor to correct the normalization
    # of a region of the spectrum, using the highest region (radius) of the spectrum 
    # insde a section contained on data
    def find_corr(self,data,radius):
        
        self.m = max(data[:,1])
        self.ind = self.find_index(self.m,data[:,1])
        self.stdev = np.std(data[self.ind-radius:self.ind+radius+1,1])
        print 'Standard deviation of the spectrum on the selected data points = %.4f' % self.stdev
        return np.mean(data[self.ind-radius:self.ind+radius+1,1])
    
    # This routine finds the multiplicative factor to correct the normalization
    # of a region of the spectrum, using an ensemble of user defined points, based on which
    # the correction will be estimated around a region within a radius of points
    def find_corr_from_ensemble(self,data,target_wls,radius):
        
        self.corrs = np.empty(len(target_wls),float)
        for i in range(len(target_wls)):
            self.ind = self.find_index(target_wls[i],data[:,1])
            self.corrs[i] = np.mean(data[self.ind-radius:self.ind+radius+1,1])
        return self.corrs