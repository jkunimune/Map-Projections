from math import inf

from numpy import pi, sqrt
from shapely import Polygon

from helpers import plot, trim_edges, load_shaperecords, ShapeRecord

SIZE_CLASSES = [
	'lg', 'md', 'sm', None, None, None]
CIRCLE_RADIUS = .7


def plot_political_shapes(filename, which="all", only_border=False, add_circles=False,
                          trim_antarctica=False, add_title=False) -> str:
	""" it's like plot_shapes but it also can make circles and handles, like, dependencies and stuff
	    :param filename: the name of the natural earth dataset to use (minus the .shp)
	    :param which: either "all" to do all countries, "big" to exclude small countries,
	                  or "small" to exclude big countries
	    :param only_border: redraw the border by copying an existing element of the same ID, and
	                        clip to that existing shape
	    :param add_circles: add circles to each small country
	    :param trim_antarctica: whether to adjust antarctica's shape
	    :param add_title: add mouseover text
	"""
	# first sort records into a dictionary by admin0 A3 code
	sovereignties: dict[str, list[ShapeRecord]] = {}
	for region in load_shaperecords(filename):
		if trim_antarctica:
			if region.record["sov_a3"] == 'ATA':  # if it is Antarctica, trim it
				region.shape.points = trim_edges(region.shape.points, region.shape.parts)
		sovereignties[region.record["sov_a3"]] = sovereignties.get(region.record["sov_a3"], []) + [region]

	# next, go thru and plot the borders
	result = ""
	for sovereignty_code, regions in sorted(sovereignties.items()):
		sovereign_code = complete_sovereign_code_if_necessary(sovereignty_code, regions)
		sovereign_content = ''
		for region in regions:
			region_code = region.record["adm0_a3"]
			is_sovereign = region_code == sovereign_code
			title = region.record["name"]

			small = True
			max_size = -inf
			for i in range(len(region.shape.parts)):
				if i + 1 < len(region.shape.parts):
					part = region.shape.points[region.shape.parts[i]:region.shape.parts[i + 1]]
				else:
					part = region.shape.points[region.shape.parts[i]:]
				# if Polygon(part).buffer(-CIRCLE_RADIUS).area == 0:
				if Polygon(part).area > pi*CIRCLE_RADIUS**2:
					small = False
				max_size = max(Polygon(part).area, max_size)

			if which == "small" and not small:
				continue
			elif which == "big" and small:
				continue

			if not only_border:
				clazz = region_code if not is_sovereign else None
				sovereign_content += plot(region.shape.points, midx=region.shape.parts, close=False,
				                          fourmat='xd', tabs=4, clazz=clazz, ident=region_code+"-shape",
				                          title=title if add_title else None)

			else:
				sovereign_content += (
					f'\t\t\t\t<clipPath id="{region_code}-clipPath">\n'
					f'\t\t\t\t\t<use href="#{region_code}-shape" />\n'
					f'\t\t\t\t</clipPath>\n'
					f'\t\t\t\t<use href="#{region_code}-shape" style="clip-path:url(#{region_code}-clipPath);" />\n'
				)

			if add_circles and small:
				capital_λ, capital_ф = float(region.record["label_x"]), float(region.record["label_y"])
				if is_sovereign:
					radius = CIRCLE_RADIUS
				else:
					radius = CIRCLE_RADIUS/sqrt(2)
				sovereign_content += f'\t\t\t\t<circle class="{region_code}" ' \
				                     f'cx="{capital_λ}" cy="{capital_ф}" r="{radius}" />\n'
				if add_title:
					sovereign_content = sovereign_content[:-4] + f'><title>{title}</title></circle>\n'

		if len(sovereign_content) > 0:
			result += f'\t\t\t<g class="{sovereign_code}">\n' + \
			          sovereign_content + \
			          f'\t\t\t</g>\n'

	return result


def complete_sovereign_code_if_necessary(sovereignty_code: str, regions: list[ShapeRecord]) -> str:
	""" so, under the NaturalEarth dataset system, countries are grouped together under
	    sovereignties.  and each sovereignty has one country within it which is in charge.  let's call it
	    the sovereign.  to encode this, the dataset gives every region a sovereign code and an admin0 code,
	    where the admin0 code is from some ISO standard, and the sovereign code = if it's the only region
	    under this sovereign { the admin0 code } else { the admin0 code of the sovereign with its last
	    letter replaced with a 1 };  I find this kind of weird; I can't find any basis for it in ISO
	    standards.  I get it, but I would rather the sovereign code just be the same as the admin0 code of
	    the sovereign.  so that's what this function does; it figures out what letter should go where that
	    1 is.
	"""
	sovereign_code = None
	for region in regions:
		if region.record["adm0_a3"][:2] == sovereignty_code[:2]:
			if sovereign_code is not None:
				if sovereignty_code == "KA1":
					return "KAZ"
				else:
					raise ValueError(f"there are two possible sovereign codes for {sovereignty_code} and I "
					                 f"don't know which to use: {sovereign_code} and {region.record['adm0_a3']}")
			else:
				sovereign_code = region.record['adm0_a3']
	if sovereign_code is None:
		raise ValueError(f"there are no possible sovereign codes for {sovereignty_code}")
	return sovereign_code
