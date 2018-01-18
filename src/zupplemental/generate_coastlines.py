import math
import numpy as np
import shapefile

from helpers import plot, is_widdershins


SOURCE = 'ne_50m'
SKIP = 1
EXPAND_ANTARCTICA = True
MAX_RANK = 2


"""data from http://www.naturalearthdata.com/"""
sf = shapefile.Reader("data/{}_coastline".format(SOURCE))
pending_coasts = {} #key is last endpoint, value is list of points; every point should appear exactly twice, in each direction
split_points = {} #key is first endpoint of master curve, value is list of indices where the curve must be split
final_coasts = []
for record, shape in zip(sf.records(), sf.shapes()):
	curve = shape.points[0:len(shape.points):SKIP]
	if curve[-1] != shape.points[-1]: 	curve.append(shape.points[-1])
	if len(curve) < 4: 	continue

	if curve[0] == curve[-1]: #this curve is complete; we're done here
		final_coasts.append(curve)
	else:
		if curve[0] in pending_coasts: #if this is an extension of another path
			curve = pending_coasts.pop(curve[0]) + curve[1:] #replace that
		elif curve[-1] in pending_coasts: #if this reversed is an extension of another path
			curve = pending_coasts.pop(curve[-1]) + list(reversed(curve[:-1])) #replace that
		pending_coasts[curve[-1]] = curve
		pending_coasts[curve[0]] = list(reversed(curve))

	# for key in pending_coasts:
	# 	assert pending_coasts[key][-1] == key, "The curve from {} to {} is keyed by {}".format(val[0], val[-1], key)
	# 	assert pending_coasts[key][0] in pending_coasts, "{} keys to a curve from {} to {}, but there is no key {}".format(key, val[0], val[-1], val[0])

for key, coast in pending_coasts.items():
	if pending_coasts[coast[0]] is not None: #if the mirror is not checked
		final_coasts.append(coast if is_widdershins(coast) else pending_coasts[coast[0]])
	pending_coasts[key] = None #mark this as checked

final_coasts = sorted(sorted(sorted(final_coasts, key=lambda c:c[0][0]), key=lambda c:c[0][1]), key=len, reverse=True)
if SOURCE == 'ne_10m':
	afroeurasia = final_coasts[1]
	split_points[afroeurasia[0]] = [len(afroeurasia)]
	final_coasts[1] += final_coasts.pop(12) #manually attach the Caspian to Afroeurasia
else:
	afroeurasia = final_coasts[0]
	for i in range(1,len(afroeurasia)):
		if math.hypot(afroeurasia[i][0]-afroeurasia[i-1][0], afroeurasia[i][1]-afroeurasia[i-1][1]) > 10:
			split_points[afroeurasia[0]] = [i] #manually detach the Caspian from Afroeurasia
			break
if EXPAND_ANTARCTICA:
	antarctica = final_coasts[2]
	final_coasts[2] = [(antarctica[0][0], y) for y in range(-90, math.ceil(antarctica[0][1]))] +\
			antarctica + [(antarctica[-1][0], y) for y in range(math.floor(antarctica[-1][1]), -91, -1)]

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

for coast in final_coasts:
	midx = [0] + split_points.get(coast[0], [])
	plot(coast, midx=midx, fourmat='xd', close=False)
