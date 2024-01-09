import math

from helpers import obliquify, plot


# IND_RAD = math.radians(3.75)
# IND_RAD = 750/6371


def plot_indicatrix(phi0, lam0, r, res) -> str:
	points = []
	for l in range(0, 360, 360//res):
		points.append(obliquify(math.pi/2-r, math.radians(l), phi0, lam0))
	return plot(points, close=True)


def plot_side_indicatrix(phi0, r, res, side_res=1) -> str:
	points = []
	midx = [0]
	for l in range(0, 181, 360//res):
		pr, lr = obliquify(math.pi/2-r, math.radians(l), phi0, 0)
		points.append((pr, lr - math.pi))
	for i in range(side_res):
		points.append((phi0+r-(2*r)/side_res*(i+1), -math.pi))
	midx.append(len(points))
	for l in range(180, 361, 360//res):
		pr, lr = obliquify(math.pi/2-r, math.radians(l), phi0, 0)
		points.append((pr, lr + math.pi))
	for i in range(side_res):
		points.append((phi0-r+(2*r)/side_res*(i+1), math.pi))
	return plot(points, midx=midx)


def plot_pole_indicatrix(north, r, res, ctr_meridian=0., side_res=1, pole_res=1) -> str:
	if north:
		pp, pr = math.pi/2, math.pi/2-r
	else:
		pp, pr = -math.pi/2, -math.pi/2+r

	points = []
	for i in range(side_res):
		points.append((pp + (pr-pp)/side_res*i, ctr_meridian-math.pi))
	for x in range(-180, 181, 360//res):
		points.append((pr, ctr_meridian+math.radians(x)))
	for i in range(side_res):
		points.append((pr + (pp-pr)/side_res*(i+1), ctr_meridian+math.pi))
	for x in range(180, -181, -360//pole_res):
		points.append((pp, ctr_meridian+math.radians(x)))
	return plot(points, close=True)


def generate_indicatrices(spacing, radius, adjust_poles=False, resolution=60, ctr_meridian=0, side_res=1, pole_res=1) -> str:
	result = ""
	result += plot_pole_indicatrix(
		True, radius, resolution, ctr_meridian=math.radians(ctr_meridian), side_res=side_res, pole_res=pole_res)
	for y in range(-90+spacing, 90, spacing):
		if adjust_poles:
			step = spacing * round(1/math.cos(math.radians(y)))
		else:
			step = spacing
		for x in range(-180, 180, step):
			if abs(x+ctr_meridian) == 180:
				result += plot_side_indicatrix(
					math.radians(y), radius, resolution, side_res=side_res)
			else:
				result += plot_indicatrix(
					math.radians(y), math.radians(x+ctr_meridian), radius, resolution)
	result += plot_pole_indicatrix(
		False, radius, resolution, ctr_meridian=math.radians(ctr_meridian), side_res=side_res, pole_res=pole_res)
	return result
