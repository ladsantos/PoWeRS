#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from pwoogs import estimate

r = estimate.rotation(np.array([6026.5,6027.5]),3.47,26,y_mult=0.996)
r.find(np.array([0.2,-0.040,-0.0120]),np.array([0.7,-0.060,-0.0140]),gamma=2.0)