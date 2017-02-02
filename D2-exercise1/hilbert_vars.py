# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 15:12:14 2015

@author: nicola
"""
import numpy as np

class Vars:
    def __init__(self):
        self.n_points = int(4096)
        self.rotation_table = np.array([3, 0, 0, 1], dtype=np.int32)
        self.sense_table = np.array([-1, 1, 1, -1], dtype=np.int32)
        self.quad_table = np.zeros((4,2,2), dtype=np.int32)#(0:3,0:1,0:1)  
        self.quad_pos = np.zeros((2,4), dtype = np.float)  
        self.corner = np.zeros((2,4), dtype = np.float)
        self.x_pos = np.zeros((self.n_points,),dtype=np.float)
        self.y_pos = np.zeros(self.x_pos.shape,dtype=np.float)
        self.levgrid = np.int32(0)
        self.levscan = np.int32(0)
        self.side = np.int32(0)
        self.key_list = np.empty((0,)) 