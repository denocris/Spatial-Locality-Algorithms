# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 09:14:57 2015

@author: nicola
"""

import numpy as np


class Vars:
    def __init__(self):
        self.x_pos = np.zeros((0,),dtype=float)
        self.y_pos = np.zeros(self.x_pos.shape)
        self.quad = np.zeros((2,4),dtype=float)
        self.corner = np.zeros((2,4),dtype=float)
        self.levgrid = 0
        self.levscan = int(4)
        self.side = int(16)
        self.key_list = np.empty((0,)) 
        return