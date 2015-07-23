#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

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