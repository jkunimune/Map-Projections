#read data from a shapefile into SVG format

import shapefile

from helpers import plot, trim_edges, lengthen_edges


def plot_shape(data, source, max_rank, trim_antarctica=False, flesh_out_antarctica=False):
	"""data from http://www.naturalearthdata.com/"""
	sf = shapefile.Reader("../data/{}_{}".format(source, data))
	for i, field in enumerate(sf.fields):
		if 'rank' in field[0]:
			rank_idx = i-1
	for record, shape in zip(sf.records(), sf.shapes()):
		if record[rank_idx] <= max_rank:
			if trim_antarctica:
				if shape.points[0][1] < -60: #if it is Antarctica (this will have a few false positives, but that's fine)
					shape.points = trim_edges(shape.points)
			elif flesh_out_antarctica:
				if shape.points[0][1] < -60:
					shape.points = lengthen_edges(shape.points)

			plot(shape.points, midx=shape.parts, close=False, fourmat='xd')

def generate_fine_borders(source, max_rank=float('inf')):
	plot_shape('admin_0_map_units', source, max_rank)

def generate_disputes(source, max_rank=float('inf')):
	plot_shape('admin_0_boundary_lines_disputed_areas', source, max_rank)

def generate_administrivia(source, max_rank=float('inf')):
	plot_shape('admin_1_states_provinces_lines', source, max_rank)

def generate_land(source, max_rank=float('inf'), trim_antarctica=False, flesh_out_antarctica=False):
	plot_shape('land', source, max_rank, trim_antarctica=trim_antarctica, flesh_out_antarctica=flesh_out_antarctica)

def generate_coastlines(source, max_rank=float('inf'), trim_antarctica=False, flesh_out_antarctica=False):
	plot_shape('coastline', source, max_rank, trim_antarctica=trim_antarctica, flesh_out_antarctica=flesh_out_antarctica)

def generate_rivers(source, max_rank=float('inf')):
	plot_shape('rivers_lake_centerlines', source, max_rank)

def generate_lakes(source, max_rank=float('inf')):
	plot_shape('lakes', source, max_rank)
