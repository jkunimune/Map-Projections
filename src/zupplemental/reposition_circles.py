import numpy as np
from matplotlib import pyplot as plt
from numpy import inf, hypot, log, pi, cos, sin, newaxis
from numpy.typing import NDArray
from scipy import optimize


TERMINATION_THRESHOLD = 1e-9


def main():
	centers_original = np.random.random((12, 2))
	radii = np.random.uniform(.05, .20, 12)

	centers = centers_original

	length_scale = np.min(radii)

	plot_circles(centers, centers_original, radii)

	# first ramp up a center-to-center repulsive force until no circles overlap
	repulsion_strength = length_scale**2
	while any_overlap(centers, radii):
		repulsion_strength *= 2

		def total_energy(x) -> float:
			current_centers = x.reshape((-1, 2))
			current_anchor_energy = calculate_anchor_energy(current_centers, centers_original)
			current_repulsive_energy = repulsion_strength*calculate_repulsive_energy(
				current_centers, radii, "center-to-center")
			return current_anchor_energy + current_repulsive_energy

		# def total_force(x) -> NDArray[float]:
		# 	current_centers = x.reshape((-1, 2))
		# 	current_anchor_force = calculate_anchor_force(current_centers, centers_original)
		# 	current_repulsive_force = repulsion_strength*calculate_repulsive_force(
		# 		current_centers, radii, "center-to-center")
		# 	return current_anchor_force + current_repulsive_force

		solution = optimize.minimize(total_energy, x0=centers.ravel())
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

		# def total_force(x):
		# 	current_centers = x.reshape((-1, 2))
		# 	current_anchor_force = calculate_anchor_force(current_centers, centers_original)
		# 	current_repulsive_force = repulsion_strength*calculate_repulsive_force(
		# 		current_centers, radii, "inner")
		# 	return current_anchor_force + current_repulsive_force

		solution = optimize.minimize(total_energy, x0=centers.ravel())
		new_centers = solution.x.reshape((-1, 2))
		change = np.sum((new_centers - centers)**2)
		centers = new_centers

		plot_circles(centers, centers_original, radii)

	print("done")
	plt.show()


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


# def calculate_anchor_force(centers: NDArray[float], centers_original: NDArray[float]):
# 	return centers_original - centers


def calculate_anchor_energy(centers: NDArray[float], centers_original: NDArray[float]):
	return 1/2*np.sum((centers_original - centers)**2)


# def calculate_repulsive_force(center: NDArray[float], radii: NDArray[float], mode: str):
# 	distance_matrix, direction_matrix = calculate_distance_matrix(center, radii, mode)
# 	distance_matrix = np.where(distance_matrix != 0, distance_matrix, REFERENCE_LENGTH_SCALE)
# 	force_matrix = direction_matrix/distance_matrix[:, :, newaxis]
# 	return np.sum(force_matrix, axis=0)


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


if __name__ == "__main__":
	main()
