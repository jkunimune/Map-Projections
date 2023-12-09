# read data from a shapefile into SVG format

from helpers import get_centroid, load_shaperecords


def plot_texts(
		filename, max_rank, label_type=None, circle_size=None,
		regulate_case=False, secondary_attr=None, force_points=False) -> str:
	result = ""
	for shape, record in load_shaperecords(filename):
		if secondary_attr is not None:
			secondary_attr_value = record[secondary_attr]
		else:
			secondary_attr_value = None
		if "labelrank" in record and "opulated" not in filename:
			rank = record["labelrank"]
		else:
			rank = record["scalerank"]  # (labelrank is broken for populated_places)
		if record["name"] == 'EUROPE':
			rank += 1  # You're not a continent! Get over it!
		if rank is None or rank > max_rank:
			continue  # skip anything too obscure to note
		if not record["name"]:
			continue  # we can't waste our time worrying about features with no names
		if "long_x" in record:
			x, y = record["long_x"], record["lat_y"]
		elif "longitude" in record:
			x, y = record["longitude"], record["latitude"]
		else:
			x, y = get_centroid(shape.points, shape.parts)
		if "type" in record:
			tipe = record["type"].lower()
		else:
			tipe = record["featurecla"].lower()
		if tipe == "depression":
			continue  # skip depressions because I don't like them

		if secondary_attr is not None and tipe == 'dependency':
			label = f"{record['name']} ({secondary_attr_value})"  # indicate dependencies' sovereigns in parentheses
		else:
			label = record["name"]
		if regulate_case:
			label = record["name"].upper()
		label = label.replace("&", "&amp;")

		# set an appropriate label class
		if label_type is not None:
			label_class = f"label-{label_type}"
		elif tipe in ["range/mtn", "foothills"]:
			label_class = "label-mountain"
		elif tipe in ["peninsula", "isthmus", "pen/cape"]:
			label_class = "label-peninsula"
		elif tipe in ["sea", "bay", "gulf", "sound"]:
			label_class = "label-sea"
		elif tipe in ["strait", "channel"]:
			label_class = "label-strait"
		elif tipe in ["tundra", "desert", "plain"]:
			label_class = "label-biome"
		else:
			label_class = f"label-{tipe.replace(' ', '-')}"

		# create the label as a <text> tag
		result += f'\t<text class="{label_class}" x="{180 + x:.03f}" y="{90 - y:.03f}">{label}</text>\n'
		# add circles to the ones that need markers
		if force_points or (tipe in ['mountain', 'pole', 'waterfall']):
			result += f'\t<circle cx="{180+x:.03f}" cy="{90-y:.03f}" r="{circle_size:.03f}" />\n'
	return result

def label_shapes(filename, max_rank=float('inf'), label_type=None, circle_size=3, secondary_attr=None) -> str:
	return plot_texts(filename, max_rank, label_type, circle_size, secondary_attr=secondary_attr)

def label_points(filename, max_rank=float('inf'), label_type=None, circle_size=3) -> str:
	return plot_texts(filename, max_rank, label_type, circle_size, force_points=True)

def generate_topographical_labels(source, circle_size, max_rank=float('inf')) -> str:
	result = ""
	result += plot_texts(source+'_geography_regions_points', max_rank, circle_size=circle_size)
	result += plot_texts(source+'_geography_regions_polys', max_rank, circle_size=circle_size, regulate_case=True)
	result += plot_texts(source+'_geography_regions_elevation_points', max_rank, circle_size=circle_size)
	result += plot_texts(source+'_geography_marine_polys', max_rank, circle_size=circle_size)
	return result
