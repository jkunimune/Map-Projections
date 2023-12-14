# read data from a shapefile into SVG format

from helpers import plot, trim_edges, lengthen_edges, load_shaperecords


def plot_shapes(
		filename, max_rank=float('inf'), clazz=None,
		trim_antarctica=False, flesh_out_antarctica=False, mark_antarctica=False,
		filter_field=None, filter_values=None) -> str:
	result = ""
	for shape, record in load_shaperecords(filename):
		if filter_field is not None:
			assert filter_values is not None
			filter_value = record[filter_field]
			if filter_value not in filter_values:
				continue  # skip it if it is filtered out
		if "labelrank" in record:
			rank = record["labelrank"]
		else:
			rank = record["scalerank"]
		if len(shape.points) == 0:
			continue  # skip it if it is empty
		if rank is None or rank > max_rank:
			continue  # skip it if it is not notable
		clazz_for_this_section = clazz
		# check if it's Antarctica (this will have a few false positives, but that's fine)
		is_antarctica = shape.points[0][1] < -60
		if is_antarctica:
			if trim_antarctica:
				shape.points = trim_edges(shape.points, shape.parts)
			elif flesh_out_antarctica:
				shape.points = lengthen_edges(shape.points)
			if mark_antarctica:
				if clazz_for_this_section is None:
					clazz_for_this_section = "antarctic"
				elif "antarctic" not in clazz_for_this_section:
					clazz_for_this_section += " antarctic"

		result += plot(shape.points, midx=shape.parts, close=False,
		               clazz=clazz_for_this_section, fourmat='xd')

	return result
