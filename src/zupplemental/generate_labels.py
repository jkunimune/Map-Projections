# read data from a shapefile into SVG format
import numpy as np
import shapefile
from matplotlib import pyplot as plt
from numpy import sqrt

from helpers import get_centroid, normalize_shaperecord


def plot_texts(
		filename, label_class, max_rank, text_size,
		regulate_case=False, secondary_attr=None, force_points=False) -> str:
	"""data from https://www.naturalearthdata.com/"""
	try:
		sf = shapefile.Reader("shapefiles/{}".format(filename))
	except shapefile.ShapefileException:
		raise FileNotFoundError(f"The shapefile {filename} is missing; please download it and put it in "
		                        f"src/zupplemental/shapefiles/")

	# look for the secondary attribute
	secd_idx = None
	if secondary_attr is not None:
		for i, field in enumerate(sf.fields):
			if field[0] == secondary_attr:
				secd_idx = i-1  # find the secondary attribute index
		if secd_idx is None:
			raise ValueError(f"I didn't find a {secondary_attr} field in this shapefile")

	result = ""
	for region in sf.shapeRecords():
		if secondary_attr is not None:
			secondary_attr_value = region.record[secd_idx]
		else:
			secondary_attr_value = None
		region = normalize_shaperecord(region)
		try:
			rank = region.record.labelrank
		except AttributeError:
			rank = region.record.scalerank
		if "opulated" in filename:
			rank = region.record.scalerank  # labelrank is broken for populated_places
		if region.record.name == 'EUROPE':
			rank += 1  # You're not a continent! Get over it!
		if rank is None or rank > max_rank:
			continue  # skip anything too obscure to note
		if not region.record.name:
			continue  # we can't waste our time worrying about features with no names
		try:
			x, y = region.record.long_x, region.record.lat_y
		except AttributeError:
			try:
				x, y = region.record.longitude, region.record.latitude
			except AttributeError:
				x, y = get_centroid(region.shape.points, region.shape.parts)
		try:
			tipe = region.record.type.lower()
		except AttributeError:
			tipe = region.record.featurecla.lower()

		if secondary_attr is not None and tipe == 'dependency':
			label = f"{region.record.name} ({secondary_attr_value})"  # indicate dependencies' sovereigns in parentheses
		else:
			label = region.record.name
		if regulate_case:
			label = region.record.name.upper()
		label = label.replace("&", "&amp;")
		label_size = ['xs', 'sm', 'md', 'lg', 'xl'][max(0, int(text_size-rank))]  # lower ranks get bigger text

		# create the label as a <text> tag
		result += f'\t<text class="label-{label_class} label-{label_size}" x="{180+x:.03f}" y="{90-y:.03f}">{label}</text>\n'
		# add circles to the ones that need markers
		if force_points or (tipe in ['mountain', 'depression', 'pole', 'waterfall']):
			circle_size = sqrt(text_size - rank + 1) if text_size != 0 else 1
			result += f'\t<circle cx="{180+x:.03f}" cy="{90-y:.03f}" r="{circle_size*0.15:.03f}" />\n'
	return result

def label_shapes(filename, label_class, max_rank=float('inf'), text_size=3, secondary_attr=None) -> str:
	return plot_texts(filename, label_class, max_rank, text_size, secondary_attr=secondary_attr)

def label_points(filename, label_class, max_rank=float('inf'), text_size=3) -> str:
	return plot_texts(filename, label_class, max_rank, text_size, force_points=True)

def generate_topographical_labels(source, max_rank=float('inf'), text_size=3) -> str:
	result = ""
	result += plot_texts(source+'_geography_regions_points', 'geo', max_rank, text_size)
	result += plot_texts(source+'_geography_regions_polys', 'geo', max_rank, text_size, regulate_case=True)
	result += plot_texts(source+'_geography_regions_elevation_points', 'geo', max_rank, text_size)
	result += plot_texts(source+'_geography_marine_polys', 'sea', max_rank, text_size)
	return result
