from math import pi, sqrt, radians, atan

import numpy as np

from helpers import obliquify, plot

PHI = (sqrt(5)+1)/2


def plot_orthodrome(phi0, lam0, tht0) -> str:
	points = []
	for p in np.linspace(-90, 90, 1001):
		points.append(obliquify(radians(p), tht0, phi0, lam0))
	return plot(points)
	

def generate_orthodromes() -> str:
	"""generate an icosahedral orthodromic mesh, like the Brilliant logo (#notsponsored)"""
	result = ""
	for l in range(-180, 180, 36):
		result += plot_orthodrome(pi/2, 0, radians(l))
	for l in range(0, 360, 72):
		result += plot_orthodrome(atan(1/2), radians(l), pi*0.2)
		result += plot_orthodrome(atan(1/2), radians(l), pi*0.4)
		result += plot_orthodrome(atan(1/2), radians(l), pi*1.2)
		result += plot_orthodrome(atan(1/2), radians(l), pi*1.4)
	return result
