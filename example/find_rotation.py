#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from pwoogs import estimate,utils

spec_window = np.array([6151.1,6152.1])
#u = utils.arr_manage()
#u.cut(spec_window[0]-10.,spec_window[1]+10.,'spectrum_full.dat','spectrum.dat')
r = estimate.rotation(spec_window,3.47,26,y_mult=0.996)
r.find(
    np.array([0.0,-0.040,-0.0120]),
    np.array([1.0,-0.060,-0.0140]),
    gamma=10.0,
    N=200,
    alpha=0.2,
    beta=0.1,
    c_limit=1E-3
    )