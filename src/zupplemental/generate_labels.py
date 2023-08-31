# read data from a shapefile into SVG format

import shapefile
from numpy import sqrt

from helpers import get_centroid


def plot_texts(
		filename, label_class, max_rank, text_size,
		regulate_case=False, secondary_attr=None, force_points=False) -> str:
	"""data from https://www.naturalearthdata.com/"""
	try:
		sf = shapefile.Reader("shapefiles/{}".format(filename))
	except shapefile.ShapefileException:
		raise FileNotFoundError(f"The shapefile {filename} is missing; please download it and put it in "
		                        f"src/zupplemental/shapefiles/")
	secd_idx = None
	for i, field in enumerate(sf.fields):
		if field[0] == secondary_attr:
			secd_idx = i-1  # find the secondary attribute index

	result = ""
	for region in sf.shapeRecords():
		try:
			rank = region.record.LABELRANK
		except AttributeError:
			rank = 0
		if rank > max_rank:
			continue  # skip anything too obscure to note
		label = region.record.NAME
		if not label:
			continue  # we can't waste our time worrying about features with no names
		try:
			x, y = region.record.long_x, region.record.lat_y
		except AttributeError:
			x, y = get_centroid(region.shape.points, region.shape.parts)
		try:
			tipe = region.record.TYPE
		except AttributeError:
			tipe = "thing"

		if region.record.NAME == 'EUROPE':
			rank += 1  # You're not a continent! Get over it!
		if regulate_case:
			label = region.record.NAME.upper()
		if secondary_attr is not None and tipe == 'Dependency':
			label = f"{region.record.NAME} ({region.record[secd_idx]})"  # indicate dependencies' sovereigns in parentheses
		label = label.replace("&", "&amp;")
		label_size = ['xs', 'sm', 'md', 'lg', 'xl'][max(0, int(text_size-rank))]  # lower ranks get bigger text

		# create the label as a <text> tag
		result += f'\t<text class="label-{label_class} label-{label_size}" x="{180+x:.03f}" y="{90-y:.03f}">{label}</text>\n'
		# add circles to the ones that need markers
		if force_points or (region.record.featurecla in ['mountain', 'depression', 'pole', 'waterfall']):
			result += f'\t<circle cx="{180+x:.03f}" cy="{90-y:.03f}" r="{sqrt(text_size-rank+1)*0.15:.03f}" />\n'
	return result

def label_shapes(filename, label_class, max_rank=float('inf'), text_size=3) -> str:
	return plot_texts(filename, label_class, max_rank, text_size, secondary_attr='NOTE_ADM0')

def label_points(filename, label_class, max_rank=float('inf'), text_size=3) -> str:
	return plot_texts(filename, label_class, max_rank, text_size, force_points=True)

def generate_topographical_labels(source, max_rank=float('inf'), text_size=3) -> str:
	result = ""
	result += plot_texts(source+'_geography_regions_points', 'geo', max_rank, text_size)
	result += plot_texts(source+'_geography_regions_polys', 'geo', max_rank, text_size, regulate_case=True)
	result += plot_texts(source+'_geography_regions_elevation_points', 'geo', max_rank, text_size)
	result += plot_texts(source+'_geography_marine_polys', 'sea', max_rank, text_size)
	return result
