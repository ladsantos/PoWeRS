#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

class arr_manage(object):
    
    def find_index(self,target,arr):
        
        self.diff = np.abs(arr-target)
        self.index = np.where(self.diff == min(self.diff))[0][0]
        return self.index
    
    def x_set_limits(self,start,end,data2d):
        
        self.start_index = self.find_index(start,data2d[:,0])
        self.end_index = self.find_index(end,data2d[:,0])
        return data2d[self.start_index:self.end_index]