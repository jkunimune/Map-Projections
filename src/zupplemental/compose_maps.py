#compose_maps.py
#make ALL the maps

import math

from generate_borders import generate_borders
from generate_graticule import generate_graticule, generate_backdrop
from generate_indicatrices import generate_indicatrices
from generate_orthodromes import generate_orthodromes
from generate_shape import plot_shapes
from generate_labels import generate_topographical_labels, label_shapes, label_points


def compose_landmasses():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="land">')
	plot_shapes('ne_50m_land', trim_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="water">')
	plot_shapes('ne_50m_lakes', max_rank=4)
	print('\t\t</g>')
	print('\t</g>')

def compose_graticule():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="graticule">')
	generate_graticule(5, 1, include_tropics=True, adjust_poles=True)
	print('\t\t</g>')
	print('\t</g>')

def compose_graticule2():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="graticule">')
	generate_graticule(15, .25, include_tropics=True, adjust_poles=True, double_dateline=True)
	print('\t\t</g>')
	print('\t</g>')

def compose_compound():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="land">')
	plot_shapes('ne_50m_land', trim_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="river">')
	plot_shapes('ne_50m_rivers_lake_centerlines', max_rank=4)
	print('\t\t</g>')
	print('\t\t<g class="lakes">')
	plot_shapes('ne_50m_lakes', max_rank=4)
	print('\t\t</g>')
	print('\t\t<g class="graticule">')
	generate_graticule(15, 1, include_tropics=True, adjust_poles=True)
	print('\t\t</g>')
	print('\t</g>')

def compose_indicatrices():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="land">')
	plot_shapes('ne_50m_land', trim_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="lakes">')
	plot_shapes('ne_50m_lakes', max_rank=4)
	print('\t\t</g>')
	print('\t\t<g class="tissot">')
	generate_indicatrices(15, math.radians(3.75), resolution=180, adjust_poles=True)
	print('\t\t</g>')
	print('\t</g>')

def compose_indicatrices2(ctr_meridian):
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="water">')
	generate_backdrop(.5, ctr_meridian=ctr_meridian)
	print('\t\t</g>')
	print('\t\t<g class="land">')
	plot_shapes('ne_110m_land', flesh_out_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="lakes">')
	plot_shapes('ne_110m_lakes')
	print('\t\t</g>')
	print('\t\t<g class="graticule">')
	generate_graticule(10, .5, double_dateline=(ctr_meridian==0))
	print('\t\t</g>')
	print('\t\t<g class="tissot">')
	generate_indicatrices(30, 500/6371, ctr_meridian=ctr_meridian, adjust_poles=True, resolution=120, side_res=5, pole_res=120)
	print('\t\t</g>')
	print('\t</g>')

def compose_political():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="country">')
	generate_borders('ne_50m', trim_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="lakes">')
	plot_shapes('ne_50m_lakes', max_rank=4)
	print('\t\t</g>')
	print('\t</g>')
	label_shapes('ne_50m_admin_0_countries', "pol")

def compose_orthodromes():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="lines">')
	generate_orthodromes()
	print('\t\t</g>')
	print('\t</g>')

def compose_everything():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="country">')
	generate_borders('ne_10m', trim_antarctica=True, borders_only=False)
	print('\t\t<g class="border">')
	generate_borders('ne_10m', trim_antarctica=True, borders_only=True)
	print('\t\t</g>')
	print('\t\t</g>')
	print('\t\t<g class="sovereign">')
	plot_shapes('ne_10m_admin_0_map_units')
	print('\t\t</g>')
	print('\t\t<g class="admin">')
	plot_shapes('ne_10m_admin_1_states_provinces_lines', filter_field='adm0_a3',
		filter_vals=['RUS','CAN','CHN','USA','BRA','AUS','IND','ARG','KAZ'])
	print('\t\t</g>')
	print('\t\t<g class="dispute">')
	plot_shapes('ne_10m_admin_0_boundary_lines_disputed_areas')
	print('\t\t</g>')
	print('\t\t<g class="coastline">')
	plot_shapes('ne_10m_coastline', trim_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="river">')
	plot_shapes('ne_10m_rivers_lake_centerlines', max_rank=5)
	print('\t\t</g>')
	print('\t\t<g class="lake">')
	plot_shapes('ne_10m_lakes', max_rank=4)
	print('\t\t</g>')
	print('\t\t<g class="graticule">')
	generate_graticule(5, 1, include_tropics=True, adjust_poles=True)
	plot_shapes('ne_10m_geographic_lines', clazz="dateline", filter_field='name', filter_vals=["International Date Line"])
	print('\t\t</g>')
	print('\t</g>')
	generate_topographical_labels('ne_50m', max_rank=2, text_size=4)
	label_shapes('ne_10m_lakes', "sea", max_rank=1, text_size=1)
	label_shapes('ne_10m_admin_0_countries', "pol", text_size=4)
	label_points('cities_capital', "cap", text_size=1)
	label_points('cities_other', "cit", text_size=0)


if __name__ == '__main__':
	# compose_landmasses()
	# compose_graticule()
	# compose_compound()
	# compose_indicatrices()
	# compose_indicatrices2(-0)
	# compose_political()
	# compose_orthodromes()
	compose_everything()
