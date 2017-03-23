from node import Node
from qtree import QuadTree

import numpy as np

def read_file(filename):

    x = np.zeros((0,))
    y = np.zeros((0,))

    f = open(filename,'r')
    aux = f.readline()
    for line in f:
        line = line.split()
        xl = float(line[0])
        yl = float(line[1])
        x = np.hstack([x,xl])
        y = np.hstack([y,yl])
    f.close()

    return x/16., y/16.


if __name__ == '__main__':

    x, y = read_file('tree.dat')

    coords_list = []

    for ix, iy in zip(x,y):
        coords = (ix,iy)
    	coords_list.append(coords)

    qtree = QuadTree(5, coords_list)

    qtree.refine(qtree.get_root_node())

    qtree.plot()
