from math import pi, sqrt, radians, atan

import numpy as np

from helpers import obliquify, plot

PHI = (sqrt(5)+1)/2


def plot_orthodrome(phi0, lam0, tht0):
	points = []
	for p in np.linspace(-90, 90, 1001):
		points.append(obliquify(radians(p), tht0, phi0, lam0))
	plot(points, close=False)
	

def generate_orthodromes():
	"""generate an icosohedral orthodromic mesh, like the Brilliant logo (#notsponsored)"""
	for l in range(-180, 180, 36):
		plot_orthodrome(pi/2, 0, radians(l))
	for l in range(0, 360, 72):
		plot_orthodrome(atan(1/2), radians(l), pi*0.2)
		plot_orthodrome(atan(1/2), radians(l), pi*0.4)
		plot_orthodrome(atan(1/2), radians(l), pi*1.2)
		plot_orthodrome(atan(1/2), radians(l), pi*1.4)
