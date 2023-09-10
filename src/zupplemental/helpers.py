from __future__ import annotations

import math
import random as rng
from typing import Any, Iterable

import shapefile


def load_shaperecords(filename) -> Iterable[ShapeRecord]:
	""" load a shapefile, with a good error message, and transmuting from the annoying pyshp api to
	    this dictier one. also, optionally, import some small shapes from a similar higher-
	    fidelity shapefile, in case you want to have all the shaperecords present even if they're
	    too small to show.
	    https://www.naturalearthdata.com/
	"""
	try:
		sf = shapefile.Reader(f"shapefiles/{filename}")
	except shapefile.ShapefileException:
		raise FileNotFoundError(f"The shapefile {filename} is missing; please download it and put it in "
		                        f"src/zupplemental/shapefiles/")
	good_shaperecords = []
	for bad_shaperecord in sf.shapeRecords():
		good_shaperecords.append(normalize_shaperecord(bad_shaperecord))
	return good_shaperecords


def load_shapes_from_one_place_and_records_from_another(shape_filename, record_filename, identifier="adm0_a3") -> Iterable[ShapeRecord]:
	""" load two shapefiles and take the shapes from one but the records from the other.  any shapes
	    in the first without corresponding records in the other will be dropd.  any records in the
	    other without corresponding shapes in the first will be assigned a NULL geometry.  this is
	    useful if you want to load low-resolution country border shapes, but also want to know which
	    small countries are being excluded so you can mark them with circles.
	"""
	shape_regions = load_shaperecords(shape_filename)
	record_regions = load_shaperecords(record_filename)
	new_regions = []
	for region in record_regions:
		record = region.record
		shape = None
		for other_region in shape_regions:
			if region.record[identifier] == other_region.record[identifier]:
				shape = other_region.shape
				break
		if shape is None:
			shape = shapefile.Shape(shapefile.NULL)
		new_regions.append(ShapeRecord(shape, record))
	return new_regions


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

	while lonf > math.pi:
		lonf -= 2*math.pi
	while lonf < -math.pi:
		lonf += 2*math.pi
		
	return latf, lonf


def plot(coords, midx=[0], close=True, fourmat='pr', clazz=None, ident=None, tabs=3) -> str:
	class_attr = f'class="{clazz}" ' if clazz is not None else ''
	ident_attr = f'id="{ident}" ' if ident is not None else ''
	tag = '\t'*tabs + f'<path {class_attr}{ident_attr}d="'
	last_move = None
	for i, coord in enumerate(coords):
		if i > 0 and coords[i-1] == coords[i]:
			continue  # no point repeating commands, though that sometimes happens
		if coord == last_move:
			tag += 'Z'
			continue  # save space by removing unnecessary repetitions

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
	if close:
		tag += 'Z'
	tag += '" />\n'
	return tag.replace('.000', '')


def trim_edges(coast, coast_parts):
	"""remove the extra points placed along the edges of the Plate Carree map"""
	for i in range(len(coast)-1, -1, -1):
		x0, y0 = coast[i-1] if i > 0 else coast[i]
		x1, y1 = coast[i]
		x2, y2 = coast[i+1] if i < len(coast)-1 else coast[i]
		if (i not in coast_parts) and (abs(x0) > 179.99 or abs(y0) > 89.99) and (abs(x1) > 179.99 or abs(y1) > 89.99) and (abs(x2) > 179.99 or abs(y2) > 89.99):
			coast.pop(i)
			for j in range(len(coast_parts)):
				if coast_parts[j] > i:
					coast_parts[j] -= 1
	return coast


def lengthen_edges(coast):
	"""add more points to any long straight bits"""
	for i in range(len(coast)-2, -1, -1):
		x0, y0 = coast[i]
		x1, y1 = coast[i+1]
		if x0 == x1 and abs(y1-y0) > 1:
			step = 1 if y0 > y1 else -1
			for j in range(math.ceil(abs(y1-y0))-1):
				coast.insert(i+1, (x1, y1+step*(j+1)))
		elif y0 == y1 and abs(x1-x0) > 1:
			step = 1 if x0 > x1 else -1
			for j in range(math.ceil(abs(x1-x0))-1):
				coast.insert(i+1, (x1+step*(j+1), y1))
	return coast


def get_centroid(points, parts=None):
	"""Compute the centroid of the spherical shape"""
	minL = min([p[0] for p in points])
	minP = min([p[1] for p in points])
	maxL = max([p[0] for p in points])
	maxP = max([p[1] for p in points])

	if maxL-minL < 15 and maxP-minP < 60:
		return ((maxL+minL)/2, (maxP+minP)/2)
	elif parts: #if there are multiple parts, try guessing the centroid of just one; the biggest one by bounding box
		parts.append(len(points))
		max_area = 0
		for i in range(1, len(parts)):
			part = points[parts[i-1]:parts[i]]
			minL = min([p[0] for p in part])
			minP = min([p[1] for p in part])
			maxL = max([p[0] for p in part])
			maxP = max([p[1] for p in part])
			if (maxL-minL)*(maxP-minP) > max_area:
				max_area = (maxL-minL)*(maxP-minP)
				best_part = (parts[i-1], parts[i])
		return get_centroid(points[best_part[0]:best_part[1]])

	lines = []
	for i in range(1, len(points)):
		lines += [(*points[i-1], *points[i])]
	xc, yc, zc = 0, 0, 0
	j = 0
	num_in = 0
	min_latnum = math.sin(math.radians(minP))
	max_latnum = math.sin(math.radians(maxP))
	min_lonnum = math.radians(minL)
	max_lonnum = math.radians(maxL)
	rng.seed(0)
	while j < 4000 or num_in < 10: #monte carlo
		j += 1
		latr = math.asin(rng.random()*(max_latnum-min_latnum)+min_latnum)
		lonr = rng.random()*(max_lonnum-min_lonnum) + min_lonnum
		latd, lond = math.degrees(latr), math.degrees(lonr)
		num_crosses = 0
		for l1, p1, l2, p2 in lines: #count the lines a northward ray crosses
			if ((l1 <= lond and l2 > lond) or (l2 <= lond and l1 > lond)) and (p1*(l2-lond)/(l2-l1) + p2*(l1-lond)/(l1-l2) > latd):
				num_crosses += 1

		if num_crosses%2 == 1: #if odd,
			num_in += 1
			xc += math.cos(latr)*math.cos(lonr) #this point is in.
			yc += math.cos(latr)*math.sin(lonr) #update the centroid
			zc += math.sin(latr)

	loncr = math.atan2(yc, xc)
	latcr = math.atan2(zc, math.hypot(xc, yc))
	return (math.degrees(loncr), math.degrees(latcr))

def normalize_shaperecord(shaperecord: shapefile.ShapeRecord) -> ShapeRecord:
	""" the Natural Earth dataset has a really horrible problem with inconsistent casing.  some of their
	    records' field names are all-lowercase, some are all-uppercase, and some are camel-case.  which is
	    especially a problem because sometimes the same field across different shapefiles from their database
	    will have different case (a country has a NAME, but a mountain only has a name).  isn't the whole
	    point that the different shapefiles are supposed to go together well?  I mean, I get that they do
	    geometricly.  but this is such an easy thing to standardize.  so I'm going to do it for them now, at
	    runtime.
	"""
	record_fields = {}
	for attr in vars(shaperecord.record)["_Record__field_positions"].keys():
		if "__" not in attr:
			record_fields[attr.lower()] = getattr(shaperecord.record, attr)
	better_shaperecord = ShapeRecord(shaperecord.shape, record_fields)
	return better_shaperecord

class ShapeRecord:
	def __init__(self, shape: shapefile.Shape, record: dict[str, Any]):
		self.shape = shape
		self.record = record
