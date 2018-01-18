import math


def obliquify(lat1, lon1, lat0, lon0):
	""" go from relative to absolute coordinates """
	latf = math.asin(math.sin(lat0)*math.sin(lat1) - math.cos(lat0)*math.cos(lon1)*math.cos(lat1))
	innerFunc = math.sin(lat1)/math.cos(lat0)/math.cos(latf) - math.tan(lat0)*math.tan(latf)
	if lat0 == math.pi/2: # accounts for special case when lat0 = pi/2
		lonf = lon1+lon0
	elif lat0 == -math.pi/2: # accounts for special case when lat0 = -pi/2
		lonf = -lon1+lon0 + math.pi
	elif abs(innerFunc) > 1: # accounts for special case when cos(lat1) -> 0
		if (lon1 == 0 and lat1 < -lat0) or (lon1 != 0 and lat1 < lat0):
			lonf = lon0 + math.pi
		else:
			lonf = lon0
	elif math.sin(lon1) > 0:
		lonf = lon0 + math.acos(innerFunc)
	else:
		lonf = lon0 - math.acos(innerFunc)
		
	return latf, lonf


def plot(coords, midx=[0], close=True, fourmat='pr'):
	tag = '\t\t\t<path d="'
	last_move = None
	for i, coord in enumerate(coords):
		if coord == last_move:
			tag += 'Z'
			continue #save space by removing unnecessary repetitions

		if i in midx:
			letter = 'M'
			last_move = coord
		else:
			letter = 'L'

		if 'r' in fourmat:
			coord = (math.degrees(c) for c in coord)
		if 'p' in fourmat:
			y, x = coord
		elif 'x' in fourmat:
			x, y = coord

		tag += '{}{:.3f},{:.3f} '.format(letter, x, y)
	tag += 'Z" />' if close else '" />'
	print(tag.replace('.000',''))


def is_widdershins(coords): #this method is not exact, but it should work well enough for my purposes.
	x1, y1 = coords[0]
	x2, y2 = coords[len(coords)//4]
	x3, y3 = coords[len(coords)*3//4]

	return (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) > 0
