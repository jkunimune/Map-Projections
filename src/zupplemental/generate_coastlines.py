import math
import numpy as np
import shapefile

from helpers import plot, is_widdershins


SOURCE = 'ne_50m'
TRIM_ANTARCTICA = False
MAX_RANK = float('inf')


"""data from http://www.naturalearthdata.com/"""
sf = shapefile.Reader("data/{}_land".format(SOURCE))
for record, shape in zip(sf.records(), sf.shapes()):
	if float(record[2]) > MAX_RANK: 	continue

	if TRIM_ANTARCTICA:
		if shape.points[0][0] == -180 and shape.points[0][1] < -80: #if it is Antarctica
			coast = shape.points
			shape.points = [(coast[0][0],y) for y in range(-90,int(coast[0][1]))] + coast + [(coast[-1][0],y) for y in range(int(coast[-1][1]),-91,-1)]
	plot(shape.points, midx=shape.parts, close=False, fourmat='xd')
