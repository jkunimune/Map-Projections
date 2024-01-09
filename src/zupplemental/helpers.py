from __future__ import annotations

import random as rng
from typing import Any, Iterable, Optional

import numpy as np
import shapefile
from matplotlib import pyplot as plt
from numpy import pi, sin, cos, tan, arcsin, arccos, degrees, ceil, radians, arctan2, hypot, cumsum
from shapefile import NULL, Shape, ShapefileException


def load_shaperecords(filename) -> Iterable[ShapeRecord]:
	""" load a shapefile, with a good error message, and transmuting from the annoying pyshp api to
	    this dictier one. also, optionally, import some small shapes from a similar higher-
	    fidelity shapefile, in case you want to have all the shaperecords present even if they're
	    too small to show.
	    https://www.naturalearthdata.com/
	"""
	try:
		sf = shapefile.Reader(f"shapefiles/{filename}")
	except ShapefileException:
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
	for shape, record in record_regions:
		shape = None
		for other_shape, other_record in shape_regions:
			if record[identifier] == other_record[identifier]:
				shape = other_shape
				break
		if shape is None:
			shape = Shape(NULL)
		new_regions.append((shape, record))
	return new_regions


def obliquify(lat1, lon1, lat0, lon0):
	""" go from relative to absolute coordinates """
	latf = arcsin(sin(lat0)*sin(lat1) - cos(lat0)*cos(lon1)*cos(lat1))
	innerFunc = sin(lat1)/cos(lat0)/cos(latf) - tan(lat0)*tan(latf)
	if lat0 == pi/2: # accounts for special case when lat0 = pi/2
		lonf = lon1+lon0
	elif lat0 == -pi/2: # accounts for special case when lat0 = -pi/2
		lonf = -lon1+lon0 + pi
	elif abs(innerFunc) > 1: # accounts for special case when cos(lat1) -> 0
		if (lon1 == 0 and lat1 < -lat0) or (lon1 != 0 and lat1 < lat0):
			lonf = lon0 + pi
		else:
			lonf = lon0
	elif sin(lon1) > 0:
		lonf = lon0 + arccos(innerFunc)
	else:
		lonf = lon0 - arccos(innerFunc)

	while lonf > pi:
		lonf -= 2*pi
	while lonf < -pi:
		lonf += 2*pi
		
	return latf, lonf


def plot(coords: list[tuple[float, float]], midx: Optional[list[int]] = None, close=False, fourmat='pr', clazz=None, ident=None, tabs=4) -> str:
	"""
	express a list of 2D points as an SVG <path> tag
	:param coords: the coordinate pairs of the vertices
	:param midx: the indices at which each part of the path starts (that is, the indices of the movetos in the path)
	:param close: whether to toss a Z on the end of the path string
	:param fourmat: the order and units of the coordinates in each pair.  probably "pr" if it's (ф,λ) coordinates in
	                radians, or "xd" if it's (λ,ф) in degrees.  the output will always be expressed as (λ,ф) in degrees,
	                so that's why I ask.
	:param clazz: the class attribute to give the <path>, if any
	:param ident: the id attribute to give the <path>, if any
	:param tabs: the number of tab characters to prepend to the <path> element
	:return: a string describing the <path> in SVG syntax, including some tab characters on the front and a newline on the back
	"""
	if midx is None:
		midx = [0]
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
			coord = (degrees(c) for c in coord)
		if 'p' in fourmat:
			y, x = coord
		elif 'x' in fourmat:
			x, y = coord
		else:
			raise ValueError(f"unrecognized format string: '{fourmat}'")

		tag += '{}{:.3f},{:.3f} '.format(letter, x, y)
	if close:
		tag += 'Z'
	tag += '" />\n'
	return tag.replace('.000', '')


def trim_edges(coast, coast_parts):
	"""remove the extra points placed along the edges of the Plate Carree map"""
	for i in range(len(coast)-1, -1, -1):
		if i in coast_parts:
			continue
		x0, y0 = coast[i-1] if i > 0 else coast[i]
		x1, y1 = coast[i]
		x2, y2 = coast[i+1] if i < len(coast)-1 else coast[i]
		if (abs(x0) > 179.99 or abs(y0) > 89.99) and \
			(abs(x1) > 179.99 or abs(y1) > 89.99) and \
			(abs(x2) > 179.99 or abs(y2) > 89.99):
			coast.pop(i)
			for j in range(len(coast_parts)):
				if coast_parts[j] > i:
					coast_parts[j] -= 1
	return coast, coast_parts


def fuse_edges(points: list[tuple[float, float]], part_indices: list[int]) -> tuple[list[tuple[float, float]], list[int]]:
	"""look for parts with matching vertices on either side of the antimeridian and combine them"""
	# first, we must apply the trimming algorithm.
	points, parts = trim_edges(points, part_indices)
	# break the data up into the individual parts
	parts = []
	for i in range(len(part_indices)):
		start = part_indices[i]
		end = part_indices[i + 1] if i + 1 < len(part_indices) else len(points)
		parts.append(points[start:end])
		if points[start] == points[end - 1]:
			parts[-1] = parts[-1][:-1]  # also remove the duplicate endpoint

	# look for pairs of parts with matching antimeridianal points, and stitch them together
	parts = fuse_edges_of_parts(parts)

	# put it all back into a single list of points and a list of moveto indices
	for part in parts:
		part.append(part[0])  # don't forget to undo the part where we removed the duplicate endpoint
	points = sum(parts, start=[])
	part_indices = cumsum([len(part) for part in parts])
	part_indices = [0] + list(part_indices[:-1])
	return points, part_indices


def fuse_edges_of_parts(parts: list[list[tuple[float, float]]]) -> list[list[tuple[float, float]]]:
	"""look for parts with matching vertices on either side of the antimeridian and combine them"""
	for i in range(len(parts)):
		# look at every point in part i
		for j in range(0, len(parts[i])):
			# search for any that are on the antimeridian, and after another on the antimeridian
			if abs(parts[i][j][0]) > 179.99 and abs(parts[i][j - 1][0]) > 179.99:
				for k in range(i):
					# look at every point in part k
					for l in range(len(parts[k])):
						# search for any that are on the antimeridian, and after another on the antimeridian
						if abs(parts[k][l][0]) > 179.99 and abs(parts[k][l - 1][0]) > 179.99:
							# if the two points' latitudes match
							if abs(parts[i][j - 1][1] - parts[k][l][1]) < 0.01:
								# do the stitching
								parts[i] = parts[i][:j] + parts[k][l:] + parts[k][:l] + parts[i][j:]
								parts.pop(k)
								# recurse in case we need to stitch together a few in a row
								return fuse_edges_of_parts(parts)
	# return the unmodified list if we found absolutely noting to fuse
	return parts


def lengthen_edges(coast):
	"""add more points to any long straight bits"""
	for i in range(len(coast)-2, -1, -1):
		x0, y0 = coast[i]
		x1, y1 = coast[i+1]
		if x0 == x1 and abs(y1-y0) > 1:
			step = 1 if y0 > y1 else -1
			for j in range(int(ceil(abs(y1-y0))-1)):
				coast.insert(i+1, (x1, y1+step*(j+1)))
		elif y0 == y1 and abs(x1-x0) > 1:
			step = 1 if x0 > x1 else -1
			for j in range(int(ceil(abs(x1-x0))-1)):
				coast.insert(i+1, (x1+step*(j+1), y1))
	return coast


def get_centroid(points, parts=None):
	"""Compute the centroid of the spherical shape"""
	minL = min([p[0] for p in points])
	minP = min([p[1] for p in points])
	maxL = max([p[0] for p in points])
	maxP = max([p[1] for p in points])

	if maxL-minL < 15 and maxP-minP < 60:
		return (maxL + minL)/2, (maxP + minP)/2
	elif parts: # if there are multiple parts, try guessing the centroid of just one; the biggest one by bounding box
		parts.append(len(points))
		max_area, best_part = 0, None
		for i in range(1, len(parts)):
			part = points[parts[i-1]:parts[i]]
			minL = min([p[0] for p in part])
			minP = min([p[1] for p in part])
			maxL = max([p[0] for p in part])
			maxP = max([p[1] for p in part])
			if (maxL-minL)*(maxP-minP) > max_area:
				max_area = (maxL-minL)*(maxP-minP)
				best_part = (parts[i-1], parts[i])
		if best_part is not None:
			return get_centroid(points[best_part[0]:best_part[1]])
		else:
			raise ValueError("no parts from which to take the centroid were found")

	lines = []
	for i in range(1, len(points)):
		lines += [(*points[i-1], *points[i])]
	xc, yc, zc = 0, 0, 0
	j = 0
	num_in = 0
	min_latnum = sin(radians(minP))
	max_latnum = sin(radians(maxP))
	min_lonnum = radians(minL)
	max_lonnum = radians(maxL)
	rng.seed(0)
	while j < 4000 or num_in < 10: # monte carlo
		j += 1
		latr = arcsin(rng.random()*(max_latnum-min_latnum)+min_latnum)
		lonr = rng.random()*(max_lonnum-min_lonnum) + min_lonnum
		latd, lond = degrees(latr), degrees(lonr)
		num_crosses = 0
		for l1, p1, l2, p2 in lines: # count the lines a northward ray crosses
			if ((l1 <= lond and l2 > lond) or (l2 <= lond and l1 > lond)) and (p1*(l2-lond)/(l2-l1) + p2*(l1-lond)/(l1-l2) > latd):
				num_crosses += 1

		if num_crosses%2 == 1: # if odd,
			num_in += 1
			xc += cos(latr)*cos(lonr) # this point is in.
			yc += cos(latr)*sin(lonr) # update the centroid
			zc += sin(latr)

	loncr = arctan2(yc, xc)
	latcr = arctan2(zc, hypot(xc, yc))
	return degrees(loncr), degrees(latcr)


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
	better_shaperecord = (shaperecord.shape, record_fields)
	return better_shaperecord


ShapeRecord = tuple[Shape, dict[str, Any]]
