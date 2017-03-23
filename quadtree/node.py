from shapely.geometry import Polygon
from shapely.geometry import LineString
from shapely.geometry import Point

import matplotlib.pylab as plt

def get_mid(p1,p2):
	c = (0.5 * (p1[0]+p2[0]), 0.5 * (p1[1]+p2[1]))
	return c

class Node:
    def __init__(self, poly):
        self.polygon = poly
        self.type = " "
        self.level = 0
        self.parent_level = -1
        self.particles_in_node = []
        return

    def generate_4_child(self):

        node_center = self.polygon.centroid.coords[0]
        point_list = self.polygon.exterior.coords

        lowl = point_list[0]
        lowr = point_list[1]
        uppr = point_list[2]
        uppl = point_list[3]

        mid_north = get_mid(uppl, uppr)
        mid_est   = get_mid(uppr, lowr)
        mid_south = get_mid(lowl, lowr)
        mid_west  = get_mid(lowl, uppl)

        NE = Node(Polygon([node_center,mid_north,uppr,mid_est]))
        SE = Node(Polygon([mid_south,node_center,mid_est,lowr]))
        SW = Node(Polygon([lowl,mid_west,node_center,mid_south]))
        NW = Node(Polygon([mid_west,uppl,mid_north,node_center]))

        node_list = [NE, SE, SW, NW]

        return node_list

    def evaluete_part_in_node(self, particles_in):
        self.particles_in_node = []
        for part in particles_in:
            if self.polygon.contains(Point(part)) == True:
                self.particles_in_node.append(part)
        return

    def plot(self):
		x = []
		y = []
		coords = self.polygon.exterior.coords
		for p in coords:
			x.append(p[0])
			y.append(p[1])

		plt.plot(x,y)
		return
