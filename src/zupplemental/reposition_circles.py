import os
import re

import numpy as np
from matplotlib import pyplot as plt
from numpy import inf, hypot, log, pi, cos, sin, newaxis
from numpy.typing import NDArray
from scipy import optimize

DIRECTORY = "../../output"
TERMINATION_THRESHOLD = 1e-9


def main():
	for filename in os.listdir(DIRECTORY):
		if filename.endswith(".svg"):
			svg, centers_original, radii = read_svg_circles(os.path.join(DIRECTORY, filename))
			if radii.size > 0:
				print(filename)
				fig = plt.figure(facecolor="none", figsize=(8, 4))
				plt.gca().invert_yaxis()
				centers_final = reposition_circles(centers_original, radii)
				os.makedirs(os.path.join(DIRECTORY, "spaced"), exist_ok=True)
				write_svg_circles(os.path.join(DIRECTORY, "spaced", filename), svg, centers_final, radii)
				print("done")
				plt.close(fig)


def reposition_circles(centers_original: NDArray[float], radii: NDArray[float]) -> NDArray[float]:
	if radii.size == 0:
		raise ValueError("the given file did not contain any circles")

	centers = centers_original

	length_scale = np.min(radii)

	plot_circles(centers, centers_original, radii)

	# first ramp up a center-to-center repulsive force until no circles overlap
	repulsion_strength = length_scale**2
	while any_overlap(centers, radii):
		repulsion_strength *= 1.5

		def total_energy(x) -> float:
			current_centers = x.reshape((-1, 2))
			current_anchor_energy = calculate_anchor_energy(current_centers, centers_original)
			current_repulsive_energy = repulsion_strength*calculate_repulsive_energy(
				current_centers, radii, "center-to-center")
			return current_anchor_energy + current_repulsive_energy

		def total_negative_force(x) -> NDArray[float]:
			current_centers = x.reshape((-1, 2))
			current_anchor_force = calculate_anchor_force(current_centers, centers_original)
			current_repulsive_force = repulsion_strength*calculate_repulsive_force(
				current_centers, radii, "center-to-center")
			return -np.ravel(current_anchor_force + current_repulsive_force)

		solution = optimize.minimize(fun=total_energy, jac=total_negative_force, x0=centers.ravel())
		centers = solution.x.reshape((-1, 2))

		plot_circles(centers, centers_original, radii)

	# then apply a gradually weakening edge-to-edge force until we reach the optimum
	change = inf
	while change > length_scale*TERMINATION_THRESHOLD:
		repulsion_strength /= 2

		def total_energy(x):
			current_centers = x.reshape((-1, 2))
			current_anchor_energy = calculate_anchor_energy(current_centers, centers_original)
			current_repulsive_energy = repulsion_strength*calculate_repulsive_energy(
				current_centers, radii, "inner")
			return current_anchor_energy + current_repulsive_energy

		def total_negative_force(x):
			current_centers = x.reshape((-1, 2))
			current_anchor_force = calculate_anchor_force(current_centers, centers_original)
			current_repulsive_force = repulsion_strength*calculate_repulsive_force(
				current_centers, radii, "inner")
			return -np.ravel(current_anchor_force + current_repulsive_force)

		solution = optimize.minimize(fun=total_energy, jac=total_negative_force, x0=centers.ravel())
		new_centers = solution.x.reshape((-1, 2))
		change = np.sum((new_centers - centers)**2)
		centers = new_centers

		plot_circles(centers, centers_original, radii)

	return centers


def plot_circles(centers: NDArray[float], centers_original: NDArray[float], radii: NDArray[float]) -> None:
	plt.clf()
	plt.axis("equal")
	θ = np.linspace(0, 2*pi)
	for i in range(len(centers)):
		plt.plot([centers[i, 0], centers_original[i, 0]],
		         [centers[i, 1], centers_original[i, 1]], "C1-")
		plt.plot(centers[i, 0] + radii[i]*cos(θ), centers[i, 1] + radii[i]*sin(θ), "k-")
	plt.pause(.50)


def any_overlap(center, radius):
	distance, _ = calculate_distance_matrix(center, radius, mode="inner")
	return np.any(distance < 0)


def calculate_anchor_force(centers: NDArray[float], centers_original: NDArray[float]):
	return centers_original - centers


def calculate_anchor_energy(centers: NDArray[float], centers_original: NDArray[float]):
	return 1/4*np.sum((centers_original - centers)**2)


def calculate_repulsive_force(center: NDArray[float], radii: NDArray[float], mode: str):
	length_scale = np.min(radii)
	distance_matrix, direction_matrix = calculate_distance_matrix(center, radii, mode)
	distance_matrix = np.where(distance_matrix != 0, distance_matrix, length_scale)
	force_matrix = direction_matrix/distance_matrix[:, :, newaxis]
	return np.sum(force_matrix, axis=0)


def calculate_repulsive_energy(centers: NDArray[float], radii: NDArray[float], mode: str) -> float:
	if mode == "inner" and any_overlap(centers, radii):
		return inf
	length_scale = np.min(radii)
	distance_matrix, _ = calculate_distance_matrix(centers, radii, mode)
	distance_matrix = np.where(distance_matrix != 0, distance_matrix, length_scale)
	energy_matrix = -log(distance_matrix/length_scale)
	return np.sum(energy_matrix)/2


def calculate_distance_matrix(centers: NDArray[float], radii: NDArray[float], mode: str
                              ) -> tuple[NDArray[float], NDArray[float]]:
	center_to_center_vector = centers[newaxis, :, :] - centers[:, newaxis, :]
	center_to_center_distance = hypot(center_to_center_vector[:, :, 0],
	                                  center_to_center_vector[:, :, 1])
	center_to_center_direction = center_to_center_vector/np.where(
		center_to_center_distance > 0, center_to_center_distance, 1)[:, :, newaxis]
	if mode == "center-to-center":
		return center_to_center_distance, center_to_center_direction
	elif mode == "inner":
		n = center_to_center_distance.shape[0]
		inner_distance = center_to_center_distance - radii[newaxis, :] - radii[:, newaxis]
		inner_distance[np.arange(n), np.arange(n)] = 0
		return inner_distance, center_to_center_direction
	else:
		raise ValueError(f"unrecognized mode: '{mode}'")


def read_svg_circles(filename: str) -> tuple[str, NDArray[float], NDArray[float]]:
	circle_pattern = r'cx="([-0-9.]+)" cy="([-0-9.]+)" r="([-0-9.]+)"'
	centers = []
	radii = []
	with open(filename) as svg_file:
		contents = svg_file.read()
	while re.search(circle_pattern, contents):
		i = len(radii)
		match = re.search(circle_pattern, contents)
		x, y, r = [float(value) for value in match.groups()]
		centers.append((x, y))
		radii.append(r)
		contents = re.sub(circle_pattern, f"<<{i}>>", contents, count=1)
	return contents, np.array(centers), np.array(radii)


def write_svg_circles(filename: str, contents: str, centers: NDArray[float], radii: NDArray[float]):
	for i in range(len(radii)):
		contents = contents.replace(
			f"<<{i}>>", f'cx="{centers[i, 0]}" cy="{centers[i, 1]}" r="{radii[i]}"')
	with open(f"{filename[:-4]}.svg", "w") as svg_file:
		svg_file.write(contents)


if __name__ == "__main__":
	main()
