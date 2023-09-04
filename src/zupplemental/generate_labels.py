# read data from a shapefile into SVG format
from numpy import sqrt

from helpers import get_centroid, load_shaperecords


def plot_texts(
		filename, label_class, max_rank, text_size,
		regulate_case=False, secondary_attr=None, force_points=False) -> str:
	result = ""
	for region in load_shaperecords(filename):
		if secondary_attr is not None:
			secondary_attr_value = region.record[secondary_attr]
		else:
			secondary_attr_value = None
		if "labelrank" in region.record and "opulated" not in filename:
			rank = region.record["labelrank"]
		else:
			rank = region.record["scalerank"]  # (labelrank is broken for populated_places)
		if region.record["name"] == 'EUROPE':
			rank += 1  # You're not a continent! Get over it!
		if rank is None or rank > max_rank:
			continue  # skip anything too obscure to note
		if not region.record["name"]:
			continue  # we can't waste our time worrying about features with no names
		if "long_x" in region.record:
			x, y = region.record["long_x"], region.record["lat_y"]
		elif "longitude" in region.record:
			x, y = region.record["longitude"], region.record["latitude"]
		else:
			x, y = get_centroid(region.shape.points, region.shape.parts)
		if "type" in region.record:
			tipe = region.record["type"].lower()
		else:
			tipe = region.record["featurecla"].lower()
		if tipe == "depression":
			continue  # skip depressions because I don't like them

		if secondary_attr is not None and tipe == 'dependency':
			label = f"{region.record['name']} ({secondary_attr_value})"  # indicate dependencies' sovereigns in parentheses
		else:
			label = region.record["name"]
		if regulate_case:
			label = region.record["name"].upper()
		label = label.replace("&", "&amp;")
		label_size = ['xs', 'sm', 'md', 'lg', 'xl'][max(0, int(text_size-rank))]  # lower ranks get bigger text

		# create the label as a <text> tag
		result += f'\t<text class="label-{label_class} label-{label_size}" x="{180+x:.03f}" y="{90-y:.03f}">{label}</text>\n'
		# add circles to the ones that need markers
		if force_points or (tipe in ['mountain', 'pole', 'waterfall']):
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
