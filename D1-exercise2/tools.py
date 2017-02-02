# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 11:44:43 2015

@author: nicola
"""

def bit_set(v, index, x):
    """Set the index-th bit of v to x, and return the new value."""
    mask = 1 << index
    v &= ~mask
    if x:
        v |= mask
    return v
    
def bit_test(i, pos):
    """test if the pos-th bit of i is equal to 1, and return a bool."""
    return bool(i & (1 << pos))
    
def write_x_y_key(filename,x,y,key):
    f = open(filename,'w')
    for ix, iy, key_m in zip(x,y,key):
        f.write(str(ix)+', '+str(iy)+', '+str(int(key_m))+'\n')
    f.close()
    
    return