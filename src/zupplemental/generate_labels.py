#read data from a shapefile into SVG format

import shapefile
import math

from helpers import line_break, get_centroid


def plot_texts(data, label_class, source, max_rank, regulate_case=False, secondary_attr=None, force_points=False, text_size=None):
	"""data from http://www.naturalearthdata.com/"""
	sf = shapefile.Reader("data/{}_{}".format(source, data))
	lat_idx = None
	for i, field in enumerate(sf.fields):
		if field[0] in ['scalerank', 'LABELRANK']:
			rank_idx = i-1
		if field[0] in ['name', 'NAME']:
			name_idx = i-1
		if field[0] == 'featurecla':
			type_idx = i-1
		if field[0] == secondary_attr:
			secd_idx = i-1
		if field[0] == 'lat_y':
			lat_idx = i-1
		if field[0] == 'long_x':
			lon_idx = i-1
	for record, shape in zip(sf.records(), sf.shapes()):
		rank = record[rank_idx]
		if rank <= max_rank:
			if lat_idx is not None:
				x, y = record[lon_idx], record[lat_idx]
			else:
				x, y = get_centroid(shape.points, shape.parts)
			label = record[name_idx]
			if not label:
				continue #we can't waste our time worrying about features with no names
			if label == 'EUROPE':
				rank += 1 #You're not a continent! Get over it!
			if regulate_case:
				label = label.upper()
			if secondary_attr is not None and record[7] == 'Dependency':
				label = "{} ({})".format(label, record[secd_idx]) #indicate dependencies
			label_size = ['sm', 'md', 'lg', 'xl'][min(3, int(max_rank-rank))] if text_size is None else text_size #lower ranks get bigger text
			print('\t<text class="label-{} label-{}" x="{:.03f}" y="{:.03f}">{}</text>'.format(
					label_class, label_size, 180+x, 90-y, label))
			if force_points or record[type_idx] in ['mountain', 'depression', 'pole', 'waterfall']: #add circles to the ones that need markers
				print('\t<circle cx="{:.03f}" cy="{:.03f}" r="{:.03f}" />'.format(180+x, 90-y, math.sqrt(max_rank-rank+1)*0.15))

def generate_political_labels(source, max_rank=4):
	plot_texts('admin_0_countries', 'pol', source, max_rank, secondary_attr='NOTE_ADM0', text_size='md')

def generate_topographical_labels(source, max_rank=2):
	plot_texts('geography_regions_points', 'geo', source, max_rank)
	plot_texts('geography_regions_polys', 'geo', source, max_rank, regulate_case=True)
	plot_texts('geography_regions_elevation_points', 'geo', source, max_rank)
	plot_texts('geography_marine_polys', 'sea', source, max_rank)

def generate_urban_labels(source, max_rank=1):
	plot_texts('populated_places', 'pol', source, max_rank, force_points=True)
