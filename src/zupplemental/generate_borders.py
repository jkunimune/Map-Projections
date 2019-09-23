import shapefile

from helpers import plot, trim_edges


SIZE_CLASSES = ['lg', 'md', 'sm', None, None, None]


def generate_borders(source, borders_only=False, trim_antarctica=False):
	"""data from http://www.naturalearthdata.com/"""
	sf = shapefile.Reader("../data/{}_admin_0_countries".format(source))
	sovereigns = {} #key is 3-char code, value is list of (record,shape)
	for record, shape in zip(sf.records(), sf.shapes()):
		if trim_antarctica:
			if record[4] == 'ATA': #if it is Antarctica
				shape.points = trim_edges(shape.points, shape.parts)

		sovereigns[record[4]] = sovereigns.get(record[4], []) + [(record, shape)]

	for code in sorted(sovereigns.keys()):
		for record, shape in sovereigns[code]:
			if record[3] == record[8]:
				country_code = record[12]
		print('\t\t\t<g id="{}">'.format(country_code))
		for record, shape in sovereigns[code]:
			x1, y1, x2, y2 = shape.bbox
			region_code = record[12]
			if borders_only:
				print('\t\t\t\t<clipPath id="{0}-clipPath">'.format(region_code))
				print('\t\t\t\t\t<use href="#{0}-shape" />'.format(region_code))
				print('\t\t\t\t</clipPath>')
				print('\t\t\t\t<use href="#{0}-shape" style="fill:none; clip-path:url(#{0}-clipPath);" />'.format(region_code))
			else:
				plot(shape.points, midx=shape.parts, close=False, fourmat='xd', tabs=4, ident="{0}-shape".format(region_code))
		print('\t\t\t</g>')
