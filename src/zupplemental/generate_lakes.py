import math
import numpy as np
import shapefile

from helpers import plot


SOURCE = 'ne_50m'
SKIP = 1
MAX_RANK = 2
INCLUDE_LAKES = False


"""data from http://www.naturalearthdata.com/"""
sf = shapefile.Reader("data/{}_lakes".format(SOURCE))
for record, lake in zip(sf.records(), sf.shapes()):
	if float(record[-2]) > MAX_RANK: #if this lake isn't important enough
		continue #skip it

	x1l, y1l, x2l, y2l = lake.bbox
	for i, continent in enumerate(final_coasts):
		xc, yc = zip(*continent)
		x1c, y1c, x2c, y2c = min(xc), min(yc), max(xc), max(yc)
		if x1c < x1l and y1c < y1l and x2c > x2l and y2c > y2l: #if this continent contains this lake
			curve = lake.points[0:len(lake.points):SKIP]
			if curve[-1] != lake.points[-1]: 	curve.append(lake.points[-1])
			split_points[continent[0]] = split_points.get(continent[0], []) + [len(continent)+p for p in lake.parts] #mark the movetos
			final_coasts[i] = continent + curve #add it
			break