import math


AXIAL_TILT = 23.43694
ANTI_TILT = 66.56306 # I have both of these because it's the easiest way to deal with roundoff


def plot_meridian(lamd, granularity, cut=0, clazz=None) -> str:
	class_attr = 'class="{}" '.format(clazz) if clazz is not None else ''
	tag = '\t\t\t<path {}d="M{},{}'.format(class_attr, lamd, cut-90) + 'v{}'.format(granularity)*round((180-2*cut)/granularity) + '" />\n'
	return tag


def plot_parallel(phid, granularity, cut=0, clazz=None) -> str:
	class_attr = 'class="{}" '.format(clazz) if clazz is not None else ''
	tag = '\t\t\t<path {}d="M{},{}'.format(class_attr, cut-180, phid) + 'h{}'.format(granularity)*round((360-2*cut)/granularity) + '" />\n'
	return tag


def generate_graticule(spacing, granularity, include_tropics=False, adjust_poles=False) -> str:
	""" Generate a mesh of latitude and longitude lines
	    :param spacing: the number of degrees between adjacent lines
	    :param granularity: the number of degrees between vertices on a line
	    :param include_tropics: whether to also include the tropic lines and polar circles
	    :param adjust_poles: whether to make most meridians stop short of the poles to reduce clutter
	"""
	result = ""
	NUM_BASE = 90//spacing
	cuts = [0]*(4*NUM_BASE)
	if adjust_poles: # if this is True, reduce the number of meridians as you approach the pole
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
	for x in [-180, 0]:
		result += plot_meridian(x, granularity, clazz="prime-m")
	result += plot_parallel(0, granularity, clazz="equator")
	if include_tropics:
		result += plot_parallel(-AXIAL_TILT, granularity, clazz="tropics")
		result += plot_parallel(AXIAL_TILT, granularity, clazz="tropics")
		result += plot_parallel(-ANTI_TILT, granularity, clazz="circles")
		result += plot_parallel(ANTI_TILT, granularity, clazz="circles")

	return result
