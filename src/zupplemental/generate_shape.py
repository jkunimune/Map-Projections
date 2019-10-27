#read data from a shapefile into SVG format

import shapefile

from helpers import plot, trim_edges, lengthen_edges


def plot_shapes(filename, max_rank=float('inf'), clazz=None, trim_antarctica=False, flesh_out_antarctica=False, filter_field=None, filter_vals=['']):
	"""data from http://www.naturalearthdata.com/"""
	sf = shapefile.Reader("shapefiles/{}".format(filename))
	rank_idx, filter_idx = None, None
	for i, field in enumerate(sf.fields):
		if 'rank' in field[0]:
			rank_idx = i-1
		if field[0] == filter_field:
			filter_idx = i-1
	for record, shape in zip(sf.records(), sf.shapes()):
		if filter_idx is not None and record[filter_idx] not in filter_vals:
			continue # skip it if it is filtered out
		if rank_idx is None or record[rank_idx] <= max_rank:
			if trim_antarctica:
				if shape.points[0][1] < -60: #if it is Antarctica (this will have a few false positives, but that's fine)
					shape.points = trim_edges(shape.points, shape.parts)
			elif flesh_out_antarctica:
				if shape.points[0][1] < -60:
					shape.points = lengthen_edges(shape.points)

			plot(shape.points, midx=shape.parts, close=False, clazz=clazz, fourmat='xd')
