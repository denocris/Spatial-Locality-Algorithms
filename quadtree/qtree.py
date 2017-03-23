from node import Node
from shapely.geometry import Polygon

import matplotlib.pylab as plt

import matplotlib.colors as colors
import matplotlib.cm as cmx

from matplotlib.collections import PatchCollection
from descartes import PolygonPatch

import numpy as np

class QuadTree:
    node_list = []
    def __init__(self, maxdepth, particles_in):
        polygon = Polygon([(0,0),(1,0),(1,1),(0,1)])
    	node = Node(polygon)
        node.particles_in_node = particles_in
        node.type = "root"
        self.node_list.append(node)
        self.max_depth = maxdepth
        self.particle_list = particles_in
        return

    def refine(self, node):
        self.bloom(node)
        self.prune()
        return

    def prune(self):
        for node in self.children:
            if node.type == 'node':
                self.node_list.append(node)
                self.refine(node)
            if node.type == 'leaf':
                self.node_list.append(node)
        return


    def bloom(self, parent):

        self.children = parent.generate_4_child()
        for child in self.children:

            child.evaluete_part_in_node(parent.particles_in_node)
            child.level = parent.level + 1
            child.parent_level = parent.level

            num_particles = len(child.particles_in_node)

            if num_particles > 1:
				child.type = 'node'
            if (num_particles > 0) and (child.level == self.max_depth) == True:
				child.type = 'leef'
            if num_particles == 1:
				child.type = 'leef'
        return

    def get_root_node(self):
        return self.node_list[0]

    def plot(self):

		fig, ax = plt.subplots()

		jet = plt.get_cmap('jet')
		values = range(0, self.max_depth + 1)
		cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

		xp = np.zeros((0,))
		yp = np.zeros((0,))

		for node in self.node_list:
			c = node.polygon.centroid
			x = c.x
			y = c.y
			colorVal = scalarMap.to_rgba(values[node.level])
			patch = PolygonPatch(node.polygon,facecolor=colorVal,edgecolor='white')
			ax.add_patch(patch)
			plt.annotate(str(node.level), xy=(x, y))
			xp = np.hstack([xp,x])
			yp = np.hstack([yp,y])

		plt.plot(xp, yp, '-o')

		x = np.zeros((0,))
		y = np.zeros((0,))

		for p in self.particle_list:
			x = np.hstack([x,p[0]])
			y = np.hstack([y,p[1]])

		plt.plot(x,y,'o')


		plt.show()
		return
