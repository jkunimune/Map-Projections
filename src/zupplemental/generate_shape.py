# read data from a shapefile into SVG format

import shapefile

from helpers import plot, trim_edges, lengthen_edges, normalize_shaperecord


def plot_shapes(filename, max_rank=float('inf'), clazz=None, trim_antarctica=False, flesh_out_antarctica=False, mark_antarctica=False, filter_field=None, filter_values=['']) -> str:
	"""data from https://www.naturalearthdata.com/"""
	try:
		sf = shapefile.Reader("shapefiles/{}".format(filename))
	except shapefile.ShapefileException:
		raise FileNotFoundError(f"The shapefile {filename} is missing; please download it and put it in "
		                        f"src/zupplemental/shapefiles/")
	filter_idx = None
	if filter_field is not None:
		for i, field in enumerate(sf.fields):
			if field[0] == filter_field:
				filter_idx = i-1
		if filter_idx is None:
			raise ValueError(f"I didn't find a {filter_field} field in this shapefile")

	result = ""

	for region in sf.shapeRecords():
		if filter_field is not None:
			filter_value = region.record[filter_idx]
			if filter_value not in filter_values:
				continue  # skip it if it is filtered out
		region = normalize_shaperecord(region)
		try:
			rank = region.record.labelrank
		except AttributeError:
			rank = region.record.scalerank
		if len(region.shape.points) == 0:
			continue  # skip it if it is empty
		if rank is None or rank > max_rank:
			continue  # skip it if it is not notable
		clazz_for_this_section = clazz
		# check if it's Antarctica (this will have a few false positives, but that's fine)
		is_antarctica = region.shape.points[0][1] < -60
		if is_antarctica:
			if trim_antarctica:
				region.shape.points = trim_edges(region.shape.points, region.shape.parts)
			elif flesh_out_antarctica:
				region.shape.points = lengthen_edges(region.shape.points)
			if mark_antarctica:
				if clazz_for_this_section is None:
					clazz_for_this_section = "antarctic"
				elif "antarctic" not in clazz_for_this_section:
					clazz_for_this_section += " antarctic"

		result += plot(region.shape.points, midx=region.shape.parts, close=False,
		               clazz=clazz_for_this_section, fourmat='xd')

	return result
