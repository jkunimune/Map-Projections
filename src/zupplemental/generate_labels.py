#read data from a shapefile into SVG format

import shapefile
import math

from helpers import line_break, get_centroid


def plot_texts(filename, label_class, max_rank, text_size, regulate_case=False, secondary_attr=None, force_points=False):
	"""data from http://www.naturalearthdata.com/"""
	sf = shapefile.Reader("../data/{}".format(filename))
	rank_idx, type_idx, secd_idx, lat_idx, lon_idx = None, None, None, None, None
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
		rank = record[rank_idx] if rank_idx is not None else 0
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
				label = "{} ({})".format(label, record[secd_idx]) if secd_idx is not None else label #indicate dependencies
			label_size = ['xs', 'sm', 'md', 'lg', 'xl'][max(0, int(text_size-rank))] #lower ranks get bigger text
			print('\t<text class="label-{} label-{}" x="{:.03f}" y="{:.03f}">{}</text>'.format(
					label_class, label_size, 180+x, 90-y, label.replace('&', '&amp;')))
			if force_points or (type_idx is not None and record[type_idx] in ['mountain', 'depression', 'pole', 'waterfall']): #add circles to the ones that need markers
				print('\t<circle cx="{:.03f}" cy="{:.03f}" r="{:.03f}" />'.format(180+x, 90-y, math.sqrt(text_size-rank+1)*0.15))

def label_shapes(filename, label_class, max_rank=float('inf'), text_size=3):
	plot_texts(filename, label_class, max_rank, text_size, secondary_attr='NOTE_ADM0')

def label_points(filename, label_class, max_rank=float('inf'), text_size=3):
	plot_texts(filename, label_class, max_rank, text_size, force_points=True)

def generate_topographical_labels(source, max_rank=float('inf'), text_size=3):
	plot_texts(source+'_geography_regions_points', 'geo', max_rank, text_size)
	plot_texts(source+'_geography_regions_polys', 'geo', max_rank, text_size, regulate_case=True)
	plot_texts(source+'_geography_regions_elevation_points', 'geo', max_rank, text_size)
	plot_texts(source+'_geography_marine_polys', 'sea', max_rank, text_size)
