#compose_maps.py
#make ALL the maps

import math

from generate_borders import generate_borders
from generate_graticule import generate_graticule, generate_backdrop
from generate_indicatrices import generate_indicatrices
from generate_orthodromes import generate_orthodromes
from generate_shape import generate_administrivia, generate_disputes, generate_fine_borders, generate_land, generate_coastlines, generate_rivers, generate_lakes
from generate_labels import generate_topographical_labels, generate_political_labels, generate_urban_labels


def compose_landmasses():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="land">')
	generate_land('ne_50m', trim_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="water">')
	generate_lakes('ne_50m', max_rank=4)
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
	generate_land('ne_50m', trim_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="river">')
	generate_rivers('ne_50m', max_rank=4)
	print('\t\t</g>')
	print('\t\t<g class="lakes">')
	generate_lakes('ne_50m', max_rank=4)
	print('\t\t</g>')
	print('\t\t<g class="graticule">')
	generate_graticule(15, 1, include_tropics=True, adjust_poles=True)
	print('\t\t</g>')
	print('\t</g>')

def compose_indicatrices():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="land">')
	generate_land('ne_50m', trim_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="lakes">')
	generate_lakes('ne_50m', max_rank=4)
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
	generate_land('ne_110m', flesh_out_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="lakes">')
	generate_lakes('ne_110m')
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
	generate_lakes('ne_50m', max_rank=4)
	print('\t\t</g>')
	print('\t</g>')
	generate_political_labels('ne_50m')

def compose_orthodromes():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="lines">')
	generate_orthodromes()
	print('\t\t</g>')
	print('\t</g>')

def compose_everything():
	print('\t<g transform="matrix(1,0,0,-1,180,90)">')
	print('\t\t<g class="country">')
	generate_borders('ne_50m', trim_antarctica=True, borders_only=False)
	print('\t\t</g>')
	print('\t\t<g class="sovereign">')
	generate_fine_borders('ne_50m')
	print('\t\t</g>')
	print('\t\t<g class="adminis">')
	generate_administrivia('ne_50m')
	print('\t\t</g>')
	print('\t\t<g class="border">')
	generate_borders('ne_50m', trim_antarctica=True, borders_only=True)
	print('\t\t</g>')
	print('\t\t<g class="dispute">')
	generate_disputes('ne_50m')
	print('\t\t</g>')
	print('\t\t<g class="coastline">')
	generate_coastlines('ne_50m', trim_antarctica=True)
	print('\t\t</g>')
	print('\t\t<g class="river">')
	generate_rivers('ne_50m')
	print('\t\t</g>')
	print('\t\t<g class="lakes">')
	generate_lakes('ne_50m')
	print('\t\t</g>')
	print('\t\t<g class="graticule">')
	generate_graticule(5, 1, include_tropics=True, adjust_poles=True)
	print('\t\t</g>')
	print('\t</g>')
	generate_topographical_labels('ne_50m')
	generate_political_labels('ne_50m')
	generate_urban_labels('ne_50m')


if __name__ == '__main__':
	# compose_landmasses()
	# compose_graticule()
	# compose_compound()
	# compose_indicatrices()
	# compose_indicatrices2(-0)
	# compose_political()
	# compose_orthodromes()
	compose_everything()
