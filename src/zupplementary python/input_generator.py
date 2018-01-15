import math
import numpy as np

PHI = (math.sqrt(5)+1)/2
ATH = math.atan(1/2)

IND_NUM = 120
IND_RAD = 500/6371


def obliquify(lat1, lon1, lat0, lon0):
	""" go from relative to absolute coordinates """
	latf = math.asin(math.sin(lat0)*math.sin(lat1) - math.cos(lat0)*math.cos(lon1)*math.cos(lat1))
	innerFunc = math.sin(lat1)/math.cos(lat0)/math.cos(latf) - math.tan(lat0)*math.tan(latf)
	if lat0 == math.pi/2: # accounts for special case when lat0 = pi/2
		lonf = lon1+lon0
	elif lat0 == -math.pi/2: # accounts for special case when lat0 = -pi/2
		lonf = -lon1+lon0 + math.pi
	elif abs(innerFunc) > 1: # accounts for special case when cos(lat1) -> 0
		if (lon1 == 0 and lat1 < -lat0) or (lon1 != 0 and lat1 < lat0):
			lonf = lon0 + math.pi
		else:
			lonf = lon0
	elif math.sin(lon1) > 0:
		lonf = lon0 + math.acos(innerFunc)
	else:
		lonf = lon0 - math.acos(innerFunc)
		
	return latf, lonf


def plot(coords):
	tag = '      <path d="'
	for phi, lam in coords:
		tag += 'L{:.3f},{:.3f} '.format(2*math.degrees(lam), 2*math.degrees(phi)+180)
	tag += 'Z"/>'
	print(tag.replace('L', 'M', 1).replace('.000',''))


def plot_orthodrome(phi0, lam0, tht0):
	tag = '    <path d="'
	for p in range(90, -90, -1):
		phin, lamn = obliquify(math.radians(p), tht0, phi0, lam0)
		tag += 'L{:.3f},{:.3f} '.format(math.degrees(lamn), math.degrees(phin))
	for p in range(-90, 90, 1):
		phin, lamn = obliquify(math.radians(p), tht0+math.pi, phi0, lam0)
		tag += 'L{:.3f},{:.3f} '.format(math.degrees(lamn), math.degrees(phin))
	tag += 'Z"/>'
	print(tag.replace('L', 'M', 1))


def plot_indicatrix(phi0, lam0, r):
	tag = '      <path d="'
	for l in range(0, 360, 360//IND_NUM):
		phin, lamn = obliquify(math.pi/2-r, math.radians(l), phi0, lam0)
		tag += 'L{:.3f},{:.3f} '.format(2*math.degrees(lamn), 2*math.degrees(phin)+180)
	tag += 'Z"/>'
	print(tag.replace('L', 'M', 1).replace('.000',''))


def plot_side_indicatrix(phi0, r):
	tag0 = '      <path d="'
	for l in range(0, 181, 360//IND_NUM):
		phin, lamn = obliquify(math.pi/2-r, math.radians(l), phi0, 0)
		tag0 += 'L{:.3f},{:.3f} '.format(2*math.degrees(lamn), 2*math.degrees(phin)+180)
	tag1 = ''
	for l in range(180, 361, 360//IND_NUM):
		phin, lamn = obliquify(math.pi/2-r, math.radians(l), phi0, 0)
		tag1 += 'L{:.3f},{:.3f} '.format(2*math.degrees(lamn)+720, 2*math.degrees(phin)+180)
	tag1 += 'Z"/>'
	print(tag0.replace('L', 'M', 1).replace('.000','') + tag1.replace('L', 'M', 1).replace('.000',''))


def plot_pole_indicatrix(north, r):
	if north:
		yp, yr = 360, 360-2*math.degrees(r)
	else:
		yp, yr = 0, 2*math.degrees(r)

	tag = '      <path d="'
	tag += 'M{},{} '.format(0, yp)
	for x in range(0, 721, 360//IND_NUM):
		tag += 'L{},{:.3f} '.format(x, yr)
	tag += 'L{},{} Z"/>'.format(720, yp)
	print(tag)


def plot_meridian(lamd):
	tag = '      <path d="'
	for phi in range(-90, 91):
		tag += 'L{},{} '.format(2*lamd, 2*phi + 180)
	tag += '"/>'
	print(tag.replace('L', 'M', 1))


def plot_parallel(phid):
	tag = '      <path d="'
	for lam in range(0, 361):
		tag += 'L{},{} '.format(2*lam, 2*phid + 180)
	tag += '"/>'
	print(tag.replace('L', 'M', 1))


if __name__ == '__main__':
	# for l in range(0, 180, 36):
	# 	plot_orthodrome(math.pi/2, 0, math.radians(l))
	# for l in range(0, 360, 72):
	# 	plot_orthodrome(ATH, math.radians(l), math.pi*.4)
	# 	plot_orthodrome(ATH, math.radians(l), math.pi*.8)

	plot_pole_indicatrix(False, IND_RAD)
	for y in [60, 30, 0, -30, -60]:
		for x in range(0, 360, 30):
			# if x == 0:
			# 	plot_side_indicatrix(math.radians(y), math.pi/36)
			# else:
				plot_indicatrix(math.radians(y), math.radians(x), IND_RAD)
	plot_pole_indicatrix(True, IND_RAD)

	# for x in range(0, 361, 10):
	# 	plot_meridian(x)
	# for y in range(-80, 90, 10):
	# 	plot_parallel(y)
