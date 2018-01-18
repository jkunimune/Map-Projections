import math

from helpers import obliquify, plot


IND_NUM = 360
# IND_RAD = math.radians(3.75)
IND_RAD = 500/6371
SPACING = 30
ADJUST_SPACING = False


def plot_indicatrix(phi0, lam0, r):
	points = []
	for l in range(0, 360, 360//IND_NUM):
		points.append(obliquify(math.pi/2-r, math.radians(l), phi0, lam0))
	plot(points)


def plot_side_indicatrix(phi0, r):
	points = []
	midx = [0]
	for l in range(0, 181, 360//IND_NUM):
		points.append(obliquify(math.pi/2-r, math.radians(l), phi0, -math.pi))
	midx.append(len(points))
	for l in range(180, 361, 360//IND_NUM):
		pr, lr = obliquify(math.pi/2-r, math.radians(l), phi0, -math.pi)
		points.append((pr, lr+2*math.pi))
	plot(points, midx=midx, close=False)


def plot_pole_indicatrix(north, r):
	if north:
		pp, pr = math.pi/2, math.pi/2-r
	else:
		pp, pr = -math.pi/2, -math.pi/2+r

	points = [(pp, -math.pi)]
	for x in range(-180, 181, 360//IND_NUM):
		points.append((pr, math.radians(x)))
	points.append((pp, math.pi))
	plot(points)


if __name__ == '__main__':
	plot_pole_indicatrix(True, IND_RAD)
	for y in range(-90+SPACING, 90, SPACING):
		if ADJUST_SPACING:
			step = SPACING * round(1/math.cos(math.radians(y)))
		else:
			step = SPACING
		for x in range(-180, 180, step):
			if abs(x) == 180:
				plot_side_indicatrix(math.radians(y), math.pi/36)
			else:
				plot_indicatrix(math.radians(y), math.radians(x), IND_RAD)
	plot_pole_indicatrix(False, IND_RAD)
