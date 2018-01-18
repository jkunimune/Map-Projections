import math


GRATICULE = 10
INCLUDE_TROPICS = False
ADJUST_POLES = False #set to True to reduce the number of meridians as you approach the pole

AXIAL_TILT = 23.43694
ANTI_TILT = 66.56306


def plot_meridian(lamd, cut=0, clazz=None):
	class_attr = 'class="{}" '.format(clazz) if clazz is not None else ''
	tag = '\t\t\t<path {}d="M{},{}'.format(class_attr, lamd, cut-90) + 'v1'*(180-2*cut) + '" />'
	print(tag)


def plot_parallel(phid, cut=0, clazz=None):
	class_attr = 'class="{}" '.format(clazz) if clazz is not None else ''
	tag = '\t\t\t<path {}d="M{},{}'.format(class_attr, cut-180, phid) + 'h1'*(360-2*cut) + '" />'
	print(tag)


if __name__ == '__main__':
	NUM_BASE = 90//GRATICULE
	cuts = [0]*(4*NUM_BASE)
	if ADJUST_POLES:
		old_num = 1
		for p in range(0, 90, GRATICULE):
			new_num = old_num*int(1/math.cos(math.radians(p+GRATICULE/2))/old_num)
			while NUM_BASE%new_num != 0: 	new_num -= 1
			if new_num >= 2*old_num:
				for i in range(len(cuts)):
					if i%old_num == 0 and i%new_num != 0:
						cuts[i] = 90-p
				old_num = new_num

	for x in range(GRATICULE, 180, GRATICULE):
		plot_meridian(-x, cut=cuts[x//GRATICULE%len(cuts)])
		plot_meridian(x, cut=cuts[x//GRATICULE%len(cuts)])
	for y in range(GRATICULE, 90, GRATICULE):
		plot_parallel(-y)
		plot_parallel(y)
	for x in [-180, 0, 180]:
		plot_meridian(x, clazz="prime-m")
	plot_parallel(0, clazz="equator")
	if INCLUDE_TROPICS:
		plot_parallel(-AXIAL_TILT, clazz="tropics")
		plot_parallel(AXIAL_TILT, clazz="tropics")
		plot_parallel(-ANTI_TILT, clazz="circles")
		plot_parallel(ANTI_TILT, clazz="circles")
