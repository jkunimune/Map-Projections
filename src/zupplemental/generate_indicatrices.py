import math

from helpers import obliquify, plot


IND_NUM = 120
# IND_RAD = math.radians(3.75)
IND_RAD = 500/6371


def plot_indicatrix(phi0, lam0, r):
	points = []
	for l in range(0, 360, 360//IND_NUM):
		points.append(obliquify(math.pi/2-r, math.radians(l), phi0, lam0))
	plot(points)


def plot_side_indicatrix(phi0, r):
	points = []
	midx = [0]
	for l in range(0, 181, 360//IND_NUM):
		points.append(obliquify(math.pi/2-r, math.radians(l), phi0, 0))
	midx.append(len(points))
	for l in range(180, 361, 360//IND_NUM):
		points.append(obliquify(math.pi/2-r, math.radians(l), phi0, 0))
	plot(points)


def plot_pole_indicatrix(north, r):
	if north:
		pp, pr = 2*math.pi, 2*math.pi-r
	else:
		pp, pr = 0, r

	points = [(pp, 0)]
	for x in range(0, 361, 360//IND_NUM):
		points.append((pr, math.radians(x)))
	points.append((pp, 2*math.pi))
	plot(points)


if __name__ == '__main__':
	plot_pole_indicatrix(True, IND_RAD)
	for y in range(-75, 90, 15):
		step = 15 * round(1/math.cos(math.radians(y)))
		if x == 0:
	 			plot_side_indicatrix(math.radians(y), math.pi/36)
	 		else:
	 			plot_indicatrix(math.radians(y), math.radians(x), IND_RAD)
	plot_pole_indicatrix(False, IND_RAD)
