import math

from helpers import obliquify, plot


# IND_RAD = math.radians(3.75)
# IND_RAD = 750/6371


def plot_indicatrix(phi0, lam0, r, res):
	points = []
	for l in range(0, 360, 360//res):
		points.append(obliquify(math.pi/2-r, math.radians(l), phi0, lam0))
	plot(points)


def plot_side_indicatrix(phi0, r, res):
	points = []
	midx = [0]
	for l in range(0, 181, 360//res):
		pr, lr = obliquify(math.pi/2-r, math.radians(l), phi0, 0)
		points.append((pr, lr - math.pi))
	midx.append(len(points))
	for l in range(180, 361, 360//res):
		pr, lr = obliquify(math.pi/2-r, math.radians(l), phi0, 0)
		points.append((pr, lr + math.pi))
	plot(points, midx=midx, close=False)


def plot_pole_indicatrix(north, r, res):
	if north:
		pp, pr = math.pi/2, math.pi/2-r
	else:
		pp, pr = -math.pi/2, -math.pi/2+r

	points = [(pp, -math.pi)]
	for x in range(-180, 181, 360//res):
		points.append((pr, math.radians(x)))
	points.append((pp, math.pi))
	plot(points)


def generate_indicatrices(spacing, radius, adjust_poles=False, resolution=60, ctr_meridian=0):
	plot_pole_indicatrix(True, radius, resolution)
	for y in range(-90+spacing, 90, spacing):
		if adjust_poles:
			step = spacing * round(1/math.cos(math.radians(y)))
		else:
			step = spacing
		for x in range(-180, 180, step):
			if abs(x+ctr_meridian) == 180:
				plot_side_indicatrix(math.radians(y), radius, resolution)
			else:
				plot_indicatrix(math.radians(y), math.radians(x+ctr_meridian), radius, resolution)
	plot_pole_indicatrix(False, radius, resolution)
