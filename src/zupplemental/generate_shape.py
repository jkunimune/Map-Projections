#read data from a shapefile into SVG format

import shapefile

from helpers import plot, trim_edges, lengthen_edges


def plot_shapes(filename, max_rank=float('inf'), trim_antarctica=False, flesh_out_antarctica=False):
	"""data from http://www.naturalearthdata.com/"""
	sf = shapefile.Reader("../data/{}".format(filename))
	for i, field in enumerate(sf.fields):
		if 'rank' in field[0]:
			rank_idx = i-1
	for record, shape in zip(sf.records(), sf.shapes()):
		if record[rank_idx] <= max_rank:
			if trim_antarctica:
				if shape.points[0][1] < -60: #if it is Antarctica (this will have a few false positives, but that's fine)
					shape.points = trim_edges(shape.points, shape.parts)
			elif flesh_out_antarctica:
				if shape.points[0][1] < -60:
					shape.points = lengthen_edges(shape.points)

			plot(shape.points, midx=shape.parts, close=False, fourmat='xd')
