import math
import numpy as np

from helpers import obliquify, plot

PHI = (math.sqrt(5)+1)/2
ATH = math.atan(1/2)


def plot_orthodrome(phi0, lam0, tht0):
	points = []
	for p in range(90, -90, -1):
		points.append(obliquify(math.radians(p), tht0, phi0, lam0))
	for p in range(-90, 90, 1):
		points.append(obliquify(math.radians(p), tht0+math.pi, phi0, lam0))
	plot(points)


if __name__ == '__main__':
	for l in range(0, 180, 36):
		plot_orthodrome(math.pi/2, 0, math.radians(l))
	for l in range(0, 360, 72):
		plot_orthodrome(ATH, math.radians(l), math.pi*.4)
		plot_orthodrome(ATH, math.radians(l), math.pi*.8)
