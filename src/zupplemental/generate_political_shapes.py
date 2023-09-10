from math import inf

import shapefile
from numpy import pi, sqrt
from shapely import Polygon

from helpers import plot, trim_edges, ShapeRecord, \
	load_shapes_from_one_place_and_records_from_another, load_shaperecords

SIZE_CLASSES = [
	'lg', 'md', 'sm', None, None, None]
CIRCLE_RADIUS = .7

SOVEREIGN_CODES = {
	"AU1": "AUS", "CH1": "CHN", "CU1": "CUB", "DN1": "DNK", "FI1": "FIN", "FR1": "FRA",
	"IS1": "ISR", "KA1": "KAZ", "GB1": "GBR", "NL1": "NLD", "NZ1": "NZL", "US1": "USA",
}
ISO_A3_TO_A2 = {
	"AFG": "AF", "AGO": "AO", "ALB": "AL", "AND": "AD", "ARE": "AE", "ARG": "AR", "ARM": "AM",
	"ATA": "AQ", "ATG": "AG", "AUS": "AU", "AUT": "AT", "AZE": "AZ", "BDI": "BI", "BEL": "BE",
	"BEN": "BJ", "BFA": "BF", "BGD": "BD", "BGR": "BG", "BHR": "BH", "BHS": "BS", "BIH": "BA",
	"BLR": "BY", "BLZ": "BZ", "BOL": "BO", "BRA": "BR", "BRB": "BB", "BRN": "BN", "BTN": "BT",
	"BWA": "BW", "CAF": "CF", "CAN": "CA", "CHN": "CN", "CHE": "CH", "CHL": "CL", "CIV": "CI",
	"CMR": "CM", "COD": "CD", "COG": "CG", "COL": "CO", "COM": "KM", "CPV": "CV", "CRI": "CR",
	"CUB": "CU", "CYN": "CY", "CYP": "CY", "CZE": "CZ", "DEU": "DE", "DJI": "DJ", "DMA": "DM",
	"DNK": "DK", "DOM": "DO", "DZA": "DZ", "ECU": "EC", "EGY": "EG", "ERI": "ER", "ESP": "ES",
	"EST": "EE", "ETH": "ET", "FIN": "FI", "FJI": "FJ", "FRA": "FR", "FSM": "FM", "GAB": "GA",
	"GBR": "GB", "GEO": "GE", "GHA": "GH", "GIN": "GN", "GMB": "GM", "GNB": "GW", "GNQ": "GQ",
	"GRC": "GR", "GRD": "GD", "GTM": "GT", "GUY": "GY", "HND": "HN", "HRV": "HR", "HTI": "HT",
	"HUN": "HU", "IDN": "ID", "IND": "IN", "IRL": "IE", "IRN": "IR", "IRQ": "IQ", "ISL": "IS",
	"ISR": "IL", "ITA": "IT", "JAM": "JM", "JOR": "JO", "JPN": "JP", "KAZ": "KZ", "KAS": "__",
	"KEN": "KE", "KGZ": "KG", "KHM": "KH", "KIR": "KI", "KNA": "KN", "KOR": "KR", "KOS": "XK",
	"KWT": "KW", "LAO": "LA", "LBN": "LB", "LBR": "LR", "LBY": "LY", "LCA": "LC", "LIE": "LI",
	"LKA": "LK", "LSO": "LS", "LTU": "LT", "LUX": "LU", "LVA": "LV", "MAR": "MA", "MCO": "MC",
	"MDA": "MD", "MDG": "MG", "MDV": "MV", "MEX": "MX", "MHL": "MH", "MKD": "MK", "MLI": "ML",
	"MLT": "MT", "MMR": "MM", "MNE": "ME", "MNG": "MN", "MOZ": "MZ", "MRT": "MR", "MUS": "MU",
	"MWI": "MW", "MYS": "MY", "NAM": "NA", "NER": "NE", "NGA": "NG", "NIC": "NI", "NLD": "NL",
	"NOR": "NO", "NPL": "NP", "NRU": "NR", "NZL": "NZ", "OMN": "OM", "PAK": "PK", "PAN": "PA",
	"PER": "PE", "PGA": "PG", "PHL": "PH", "PLW": "PW", "PNG": "PG", "POL": "PL", "PRK": "KP",
	"PRT": "PT", "PRY": "PY", "QAT": "QA", "ROU": "RO", "RUS": "RU", "RWA": "RW", "SAH": "EH",
	"SAU": "SA", "SDN": "SD", "SDS": "SS", "SEN": "SN", "SGP": "SG", "SLB": "SB", "SLE": "SL",
	"SLV": "SV", "SMR": "SM", "SOL": "__", "SOM": "SO", "SRB": "RS", "STP": "ST", "SUR": "SR",
	"SVK": "SK", "SVN": "SI", "SWE": "SE", "SWZ": "SZ", "SYC": "SC", "SYR": "SY", "TCD": "TD",
	"TGO": "TG", "THA": "TH", "TJK": "TJ", "TKM": "TM", "TLS": "TL", "TON": "TO", "TTO": "TT",
	"TUN": "TN", "TUR": "TR", "TUV": "TV", "TWN": "TW", "TZA": "TZ", "UGA": "UG", "UKR": "UA",
	"URY": "UY", "USA": "US", "UZB": "UZ", "VAT": "VA", "VCT": "VC", "VEN": "VE", "VNM": "VN",
	"VUT": "VU", "WSM": "WS", "YEM": "YE", "ZAF": "ZA", "ZMB": "ZM", "ZWE": "ZW", }


def plot_political_shapes(filename, which="all", only_border=False, add_circles=False,
                          trim_antarctica=False, add_title=False, include_circles_from=None) -> str:
	""" it's like plot_shapes but it also can make circles and handles, like, dependencies and stuff
	    :param filename: the name of the natural earth dataset to use (minus the .shp)
	    :param which: either "all" to do all countries, "big" to exclude small countries,
	                  or "small" to exclude big countries
	    :param only_border: redraw the border by copying an existing element of the same ID, and
	                        clip to that existing shape
	    :param add_circles: add circles to each small country
	    :param trim_antarctica: whether to adjust antarctica's shape
	    :param add_title: add mouseover text
	    :param include_circles_from: an additional filename from which to borrow small objects. anything
	                                 present in this shapefile but not in the main one will be included
	                                 in the output as only a circle (assuming add_circles is True)
	"""
	if include_circles_from is None:
		regions = load_shaperecords(filename)
	else:
		regions = load_shapes_from_one_place_and_records_from_another(filename, include_circles_from)
	# go thru the records and sort them into a dictionary by ISO 3166 codes
	hierarchially_arranged_regions: dict[tuple[str, ...], ShapeRecord] = {}
	for region in regions:
		# if it Antarctica, trim it
		if trim_antarctica:
			if region.record["sov_a3"] == 'ATA':
				region.shape.points = trim_edges(region.shape.points, region.shape.parts)
		sovereign_code = region.record["sov_a3"]
		if sovereign_code in SOVEREIGN_CODES:
			sovereign_code = SOVEREIGN_CODES[sovereign_code]  # convert US1 to USA
		if "iso_3166_2" in region.record:  # either key them by sovereign, admin0, admin1
			sovereign_code = ISO_A3_TO_A2[sovereign_code]
			province_code = region.record["iso_3166_2"]
			if province_code.startswith("-99"):
				province_code = "__" + province_code[3:]
			if province_code.endswith("~"):
				province_code = province_code[:-1]
			country_code = province_code[:province_code.index("-")]
			key = (sovereign_code, country_code, province_code)
		else:  # or by sovereign, admin0
			country_code = region.record["adm0_a3"]
			key = (sovereign_code, country_code)
		if key[0] == key[1]:  # remove duplicate layers
			key = key[1:]
		hierarchially_arranged_regions[key] = region

	# next, go thru and plot the borders
	current_state = []
	result = ""
	# for each item
	for key in sorted(hierarchially_arranged_regions.keys()):
		region = hierarchially_arranged_regions[key]

		# decide whether it's "small"
		is_small = True
		if region.shape.shapeType != shapefile.NULL:
			max_size = -inf
			for i in range(len(region.shape.parts)):
				if i + 1 < len(region.shape.parts):
					part = region.shape.points[region.shape.parts[i]:region.shape.parts[i + 1]]
				else:
					part = region.shape.points[region.shape.parts[i]:]
				# if Polygon(part).buffer(-CIRCLE_RADIUS).area == 0:
				if Polygon(part).area > pi*CIRCLE_RADIUS**2:
					is_small = False
				max_size = max(Polygon(part).area, max_size)

		if which == "small" and not is_small:
			continue
		elif which == "big" and is_small:
			continue

		# make some other decisions
		is_sovereign = len(key) == 1  # this won't work for admin-1-states-provinces but that's fine
		title = region.record["name"]
		is_inhabited = (region.record.get("pop_est", inf) > 500 and  # Vatican is inhabited but US Minor Outlying I. are not
		                region.record["type"] != "Lease" and  # don't circle Baykonur or Guantanamo
		                "Base" not in region.record["admin"])  # don't circle military bases

		# exit any <g>s we're no longer in
		while current_state and (len(current_state) > len(key) or current_state[-1] != key[len(current_state) - 1]):
			result += '\t'*(2 + len(current_state)) + f'</g>\n'
			current_state.pop()
		# enter any new <g>s
		while len(current_state) < len(key):
			current_state.append(key[len(current_state)])
			result += '\t'*(2 + len(current_state)) + f'<g class="{current_state[-1]}">\n'

		# then put in whatever type of content is appropriate:
		# the normal polygon
		if not only_border:
			result += plot(region.shape.points, midx=region.shape.parts, close=False,
			               fourmat='xd', tabs=3 + len(current_state), ident=key[-1],
			               title=title if add_title else None)
		# or the clipped and copied thick border
		else:
			indentation = '\t'*(3 + len(current_state))
			result += (
				f'{indentation}<clipPath id="{key[-1]}-clipPath">\n'
				f'{indentation}<use href="#{key[-1]}" />\n'
				f'{indentation}</clipPath>\n'
				f'{indentation}<use href="#{key[-1]}" style="clip-path:url(#{key[-1]}-clipPath);" />\n'
			)
		# and potentially also a circle
		if add_circles and is_small and is_inhabited:
			x_center, y_center = float(region.record["label_x"]), float(region.record["label_y"])
			if is_sovereign:
				radius = CIRCLE_RADIUS
			else:
				radius = CIRCLE_RADIUS/sqrt(2)
			indentation = '\t'*(3 + len(current_state))
			result += f'{indentation}<circle id="{key[-1]}-circle" cx="{x_center}" cy="{y_center}" r="{radius}" />\n'
			if add_title:  # TODO: move circle just into group
				result = result[:-4] + f'><title>{title}</title></circle>\n'

	# exit all <g>s before returning
	while len(current_state) > 0:
		result += '\t'*(2 + len(current_state)) + f'</g>\n'
		current_state.pop()
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
