import math


AXIAL_TILT = 23.43694
ANTI_TILT = 66.56306 #I have both of these because it's the easiest way to deal with roundoff


def plot_meridian(lamd, granularity, cut=0, clazz=None) -> str:
	class_attr = 'class="{}" '.format(clazz) if clazz is not None else ''
	tag = '\t\t\t<path {}d="M{},{}'.format(class_attr, lamd, cut-90) + 'v{}'.format(granularity)*round((180-2*cut)/granularity) + '" />\n'
	return tag


def plot_parallel(phid, granularity, cut=0, clazz=None) -> str:
	class_attr = 'class="{}" '.format(clazz) if clazz is not None else ''
	tag = '\t\t\t<path {}d="M{},{}'.format(class_attr, cut-180, phid) + 'h{}'.format(granularity)*round((360-2*cut)/granularity) + '" />\n'
	return tag


def generate_graticule(spacing, granularity, include_tropics=False, adjust_poles=False, double_dateline=False) -> str:
	"""Generate a mesh of latitude and longitude lines"""
	result = ""
	NUM_BASE = 90//spacing
	cuts = [0]*(4*NUM_BASE)
	if adjust_poles: #if this is True, reduce the number of meridians as you approach the pole
		old_num = 1
		for p in range(0, 90, spacing):
			new_num = old_num*int(1/math.cos(math.radians(p+spacing/2))/old_num)
			while NUM_BASE%new_num != 0:
				new_num -= 1
			if new_num >= 2*old_num:
				for i in range(len(cuts)):
					if i%old_num == 0 and i%new_num != 0:
						cuts[i] = 90-p
				old_num = new_num

	for x in range(spacing, 180, spacing):
		result += plot_meridian(-x, granularity, cut=cuts[x//spacing%len(cuts)])
		result += plot_meridian(x, granularity, cut=cuts[x//spacing%len(cuts)])
	for y in range(spacing, 90, spacing):
		result += plot_parallel(-y, granularity)
		result += plot_parallel(y, granularity)
	for x in [-180, 0, 180] if double_dateline else [-180, 0]:
		result += plot_meridian(x, granularity, clazz="prime-m")
	result += plot_parallel(0, granularity, clazz="equator")
	if include_tropics:
		result += plot_parallel(-AXIAL_TILT, granularity, clazz="tropics")
		result += plot_parallel(AXIAL_TILT, granularity, clazz="tropics")
		result += plot_parallel(-ANTI_TILT, granularity, clazz="circles")
		result += plot_parallel(ANTI_TILT, granularity, clazz="circles")

	return result


def generate_backdrop(resolution, ctr_meridian=0):
	"""generate a backdrop type thing that will only work in standard aspect"""
	left_side = (ctr_meridian+.001)%360 - 180
	middle = (ctr_meridian+180)%360 - 180
	right_side = (ctr_meridian-.001)%360 - 180
	tag = '\t\t\t<path d="M{:.3f},-90'.format(left_side) + 'v{}'.format(resolution)*int(180//resolution) +\
			'L{:.3f},90L{:.3f},90'.format(middle, right_side) + 'v-{}'.format(resolution)*int(180//resolution) +\
			'L{:.3f},-90 Z" />\n'.format(middle)
	return tag
