import re
from typing import Optional, Any

import shapefile
from numpy import pi, sqrt, cos, radians, inf
from shapely import Polygon

from helpers import plot, trim_edges, ShapeRecord, \
	load_shapes_from_one_place_and_records_from_another, load_shaperecords, fuse_edges

SIZE_CLASSES = [
	'lg', 'md', 'sm', None, None, None]
CIRCLE_RADIUS = .7

ISO_A3_CODES_THAT_DONT_NEED_BE_UNIQUE = {
	"AU1": "AUS", "CH1": "CHN", "CU1": "CUB", "DN1": "DNK", "FI1": "FIN", "FR1": "FRA",
	"IS1": "ISR", "KA1": "KAZ", "GB1": "GBR", "NL1": "NLD", "NZ1": "NZL", "PN1": "PNG",
	"PR1": "PRT", "US1": "USA",
	"FXX": "FRA", "NOW": "NOR", "PNX": "PNG", "SRS": "SRB", "SYX": "SYR",
}
ISO_A3_TO_A2 = {  # why did I think hard-coding this table was a good idea...
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
	"ISR": "IL", "ITA": "IT", "JAM": "JM", "JOR": "JO", "JPN": "JP", "KAZ": "KZ", "KEN": "KE",
	"KGZ": "KG", "KHM": "KH", "KIR": "KI", "KNA": "KN", "KOR": "KR", "KOS": "XK", "KWT": "KW",
	"LAO": "LA", "LBN": "LB", "LBR": "LR", "LBY": "LY", "LCA": "LC", "LIE": "LI", "LKA": "LK",
	"LSO": "LS", "LTU": "LT", "LUX": "LU", "LVA": "LV", "MAR": "MA", "MCO": "MC", "MDA": "MD",
	"MDG": "MG", "MDV": "MV", "MEX": "MX", "MHL": "MH", "MKD": "MK", "MLI": "ML", "MLT": "MT",
	"MMR": "MM", "MNE": "ME", "MNG": "MN", "MOZ": "MZ", "MRT": "MR", "MUS": "MU", "MWI": "MW",
	"MYS": "MY", "NAM": "NA", "NER": "NE", "NGA": "NG", "NIC": "NI", "NLD": "NL", "NOR": "NO",
	"NPL": "NP", "NRU": "NR", "NZL": "NZ", "OMN": "OM", "PAK": "PK", "PAN": "PA", "PER": "PE",
	"PGA": "SP", "PHL": "PH", "PLW": "PW", "PNG": "PG", "POL": "PL", "PRK": "KP", "PRT": "PT",
	"PRY": "PY", "QAT": "QA", "ROU": "RO", "RUS": "RU", "RWA": "RW", "SAH": "EH", "SAU": "SA",
	"SDN": "SD", "SDS": "SS", "SEN": "SN", "SGP": "SG", "SLB": "SB", "SLE": "SL", "SLV": "SV",
	"SMR": "SM", "SOM": "SO", "SRB": "RS", "STP": "ST", "SUR": "SR", "SVK": "SK", "SVN": "SI",
	"SWE": "SE", "SWZ": "SZ", "SYC": "SC", "SYR": "SY", "TCD": "TD", "TGO": "TG", "THA": "TH",
	"TJK": "TJ", "TKM": "TM", "TLS": "TL", "TON": "TO", "TTO": "TT", "TUN": "TN", "TUR": "TR",
	"TUV": "TV", "TWN": "TW", "TZA": "TZ", "UGA": "UG", "UKR": "UA", "URY": "UY", "USA": "US",
	"UZB": "UZ", "VAT": "VA", "VCT": "VC", "VEN": "VE", "VNM": "VN", "VUT": "VU", "WSM": "WS",
	"YEM": "YE", "ZAF": "ZA", "ZMB": "ZM", "ZWE": "ZW",
	"NFK": "NF", "HKG": "HK", "MAC": "MO", "GRL": "GL", "FRO": "FO", "ALD": "AX",
	"ATF": "TF", "BLM": "BL", "MAF": "MF", "NCL": "NC", "PYF": "PF", "SPM": "PM", "WLF": "WF",
	"AIA": "AI", "BMU": "BM", "CYM": "KY", "FLK": "FK", "GGY": "GG", "GIB": "GI", "IMN": "IM",
	"IOT": "IO", "JEY": "JE", "MSR": "MS", "SHN": "SH", "TCA": "TC", "VGB": "VG", "PSX": "PS",
	"ABW": "AW", "CUW": "CW", "SXM": "SX", "COK": "CK", "NIU": "NU", "ASM": "AS", "GUM": "GU",
	"MNP": "MP", "PRI": "PR", "VIR": "VI", "HMD": "HM", "PCN": "PN", "SGS": "GS",
	"GUF": "GF", "BVT": "BV", "CCK": "CC", "CXR": "CX", "GLP": "GP", "MTQ": "MQ", "MYT": "YT",
	"REU": "RE", "NLY": "BQ", "NSV": "SJ", "TKL": "TK", "SYU": "XD",
	"CLP": "CP", "UMI": "UM",
	"EUE": "EU",
}
SUPRANATIONS = {
	"EUE": {
		"AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GRC", "HUN",
		"IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "ESP",
		"SWE",
	}
}


def plot_political_shapes(filename, mode="normal",
                          trim_antarctica=False, fuse_russia=False,
                          add_title=False, include_circles_from=None, filter_field: Optional[str] = None,
                          filter_values: Optional[list[Any]] = None, filter_mode="in") -> str:
	""" it's like plot_shapes but it also can make circles and handles, like, dependencies and stuff
	    :param filename: the name of the natural earth dataset to use (minus the .shp)
	    :param mode: one of "normal" to draw countries' shapes as one would expect,
	                        "trace" to redraw each border by copying an existing element of the same
	                        ID and clip to that existing shape, or
	                        "circle" to do circles at the center of mass specificly for small countries
	                        "bubble" to also do circles but their area is proportional to population
	    :param trim_antarctica: whether to adjust antarctica's shape
	    :param fuse_russia: whether to reattach the bit of Russia that's in the western hemisphere
	    :param add_title: add mouseover text
	    :param include_circles_from: an additional filename from which to borrow small objects. anything
	                                 present in this shapefile but not in the main one will be included
	                                 in the output as only a circle (assuming add_circles is True)
	    :param filter_field: a record key that will be used to select a subset of records to be used
	    :param filter_values: a list of values that will be compared against each record's
	                          filter_field to determine whether to use it
	    :param filter_mode: if "in", use only records whose filter_field value is in filter_values.
	                        if "out", use only records whose filter_field value is *not* in filter_values.
	"""
	if include_circles_from is None:
		regions = load_shaperecords(filename)
	else:
		regions = load_shapes_from_one_place_and_records_from_another(filename, include_circles_from)
	# go thru the records and sort them into a dictionary by ISO 3166 codes
	hierarchically_arranged_regions: dict[tuple[tuple[str, ...], str], ShapeRecord] = {}
	for shape, record in regions:
		# if it is filtered out, skip it
		if filter_field is not None:
			assert filter_values is not None
			filter_value = record[filter_field]
			if filter_mode == "in":
				if filter_value not in filter_values:
					continue
			elif filter_mode == "out":
				if filter_value in filter_values:
					continue
			else:
				raise ValueError(f"unrecognized filter_mode: '{filter_mode}' (must be 'in' or 'out')")

		# if it Antarctica, trim it
		if trim_antarctica:
			if record["sov_a3"] == "ATA":
				shape.points, shape.parts = trim_edges(shape.points, shape.parts)
		# fuse edges across the antimeridian if desired
		if fuse_russia:
			if record["sov_a3"] == "RUS":
				shape.points, shape.parts = fuse_edges(shape.points, shape.parts)

		# normalize some alpha3 codes
		for key in ["sov_a3", "adm0_a3", "gu_a3"]:
			if record[key] in ISO_A3_CODES_THAT_DONT_NEED_BE_UNIQUE:
				record[key] = ISO_A3_CODES_THAT_DONT_NEED_BE_UNIQUE[record[key]]

		# then figure out the keying
		sovereign_code = record["sov_a3"]
		country_code = record["adm0_a3"]
		geounit_code = record["gu_a3"]
		unique_identifier = geounit_code
		# key them by sovereign, admin0, geounit
		hierarchical_identifier = (sovereign_code, country_code, geounit_code)
		# and then by admin1 if applicable
		if "iso_3166_2" in record:
			province_code = record["iso_3166_2"]
			if province_code.startswith("-99"):
				province_code = "__" + province_code[3:]
			if province_code.endswith("~"):
				province_code = province_code[:-1]
			unique_identifier = record["adm1_code"]
			hierarchical_identifier = hierarchical_identifier + (province_code,)
		# put nations under supranations when applicable
		for code, member_codes in SUPRANATIONS.items():
			if sovereign_code in member_codes:
				hierarchical_identifier = (code,) + hierarchical_identifier
		# remove redundant layers
		for i in range(len(hierarchical_identifier) - 1, 0, -1):
			if hierarchical_identifier[i] == hierarchical_identifier[i - 1]:
				hierarchical_identifier = hierarchical_identifier[:i] + hierarchical_identifier[i + 1:]
			elif hierarchical_identifier[i - 1] in ISO_A3_TO_A2:
				if hierarchical_identifier[i] == ISO_A3_TO_A2[hierarchical_identifier[i - 1]]:
					hierarchical_identifier = hierarchical_identifier[:i] + hierarchical_identifier[i + 1:]
		# key it in this dictionary by both its classes and its id
		key = (hierarchical_identifier, unique_identifier)
		hierarchically_arranged_regions[key] = (shape, record)

	# next, go thru and plot the borders
	current_state = []
	result = ""
	already_titled = set()
	# for each item
	for key in sorted(hierarchically_arranged_regions.keys()):
		hierarchy, identifier = key
		shape, record = hierarchically_arranged_regions[key]

		# decide whether it's "small"
		is_small = True
		for i in range(len(shape.parts)):
			if i + 1 < len(shape.parts):
				part = shape.points[shape.parts[i]:shape.parts[i + 1]]
			else:
				part = shape.points[shape.parts[i]:]
			area = Polygon(part).area*cos(radians(part[0][1]))
			if area > pi*CIRCLE_RADIUS**2:
				is_small = False

		# make some other decisions
		has_geometry = shape.shapeType != shapefile.NULL
		is_sovereign = len(hierarchy) == 1 or (len(hierarchy) == 2 and hierarchy[0] in SUPRANATIONS)
		is_inhabited = (record.get("pop_est", inf) > 500 and  # Vatican is inhabited but US Minor Outlying I. are not
		                record["type"] != "Lease" and  # don't circle Baykonur or Guantanamo
		                "Base" not in record["admin"] and  # don't circle military bases
		                (record["type"] != "Indeterminate" or record["pop_est"] > 100_000))  # don't circle Cyprus No Mans Land or Siachen Glacier
		if not is_sovereign and "note_adm0" in record and record["note_adm0"] != "":
			label = f"{record['name_long']} ({record['note_adm0']})"  # indicate dependencies' sovereigns in parentheses
		elif "name_long" in record:
			label = record["name_long"]
		else:
			label = record["name"]

		# exit any <g>s we're no longer in
		while current_state and (len(current_state) > len(hierarchy) or current_state[-1] != hierarchy[len(current_state) - 1]):
			result += '\t'*(3 + len(current_state)) + f'</g>\n'
			current_state.pop()
		# enter any new <g>s
		while len(current_state) < len(hierarchy):
			current_state.append(hierarchy[len(current_state)])
			clazz = current_state[-1]
			if current_state[-1] in ISO_A3_TO_A2.keys():
				clazz += " " + ISO_A3_TO_A2[clazz]
			elif len(clazz) < 4 and record["iso_a2"] != "-99":
				print(f"warning: no ISO 3166 alpha2 code found for {clazz}: {record['name']}, {record['iso_a2']}")
			result += '\t'*(3 + len(current_state)) + f'<g class="{clazz}">\n'
		indentation = '\t'*(4 + len(current_state))

		# then put in whatever type of content is appropriate:
		any_content = False
		# the normal polygon
		if mode == "normal":
			if has_geometry:
				result += plot(shape.points, midx=shape.parts, close=False,
				               fourmat='xd', tabs=4 + len(current_state), ident=identifier)
				any_content = True
		# or the clipped and copied thick border
		elif mode == "trace":
			if has_geometry:
				result += (
					f'{indentation}<clipPath id="{identifier}-clipPath">\n'
					f'{indentation}<use href="#{identifier}" />\n'
					f'{indentation}</clipPath>\n'
					f'{indentation}<use href="#{identifier}" style="clip-path:url(#{identifier}-clipPath);" />\n'
				)
				any_content = True
		# or just a circle
		elif mode == "circle":
			if is_small and is_inhabited:
				x_center, y_center = float(record["label_x"]), float(record["label_y"])
				if is_sovereign:
					radius = CIRCLE_RADIUS
				else:
					radius = round(CIRCLE_RADIUS/sqrt(2), 2)
				result += f'{indentation}<circle id="{identifier}-circle" cx="{x_center:.3f}" cy="{y_center:.3f}" r="{radius}" />\n'
				any_content = True
		# or a circle with its area set to the population
		elif mode == "bubble":
			x_center = record["label_x"] if "label_x" in record else record["longitude"]
			y_center = record["label_y"] if "label_y" in record else record["latitude"]
			area = record["pop_est"] if "pop_est" in record else 1e8
			result += f'{indentation}<circle id="{identifier}-bubble" cx="{x_center:.3f}" cy="{y_center:.3f}" r="{2e-4*sqrt(area):.3f}" />\n'
			any_content = True

		# also a title if that's desired
		if add_title and any_content and tuple(hierarchy) not in already_titled:
			if has_geometry or is_inhabited:
				result += f'{indentation}<title>{label}</title>\n'
				already_titled.add(tuple(hierarchy))  # occasionally a thing can get two titles if the hierarchy isn't unique; only label the first one

	# exit all <g>s before returning
	while len(current_state) > 0:
		result += '\t'*(3 + len(current_state)) + f'</g>\n'
		current_state.pop()

	# remove any groups with no elements
	for _ in range(3):
		result = re.sub(r'(\t)*<g class="([A-Za-z _-]+)">(\s*)</g>\n',
		                '', result)
	# simplify any groups with only a single element
	result = re.sub(r'<g class="([A-Za-z _-]+)">(\s*)<([a-z]+) ([^\n]*)>(\s*)</g>',
	                '<\\3 class="\\1" \\4>', result)

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
	for shape, record in regions:
		if record["adm0_a3"][:2] == sovereignty_code[:2]:
			if sovereign_code is not None:
				if sovereignty_code == "KA1":
					return "KAZ"
				else:
					raise ValueError(f"there are two possible sovereign codes for {sovereignty_code} and I "
					                 f"don't know which to use: {sovereign_code} and {record['adm0_a3']}")
			else:
				sovereign_code = record['adm0_a3']
	if sovereign_code is None:
		raise ValueError(f"there are no possible sovereign codes for {sovereignty_code}")
	return sovereign_code
