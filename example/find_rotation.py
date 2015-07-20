#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from pwoogs import rotation_CE

r = rotation_CE.rotation(np.array([6026.5,6027.5]),0.0,-0.0135,0.0,0.004,3.47,26)
r.find(np.array([0.2,-0.040]),np.array([0.7,-0.070]))