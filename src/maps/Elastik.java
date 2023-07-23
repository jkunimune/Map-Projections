/**
 * MIT License
 * <p>
 * Copyright (c) 2018 Justin Kunimune
 * <p>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * <p>
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package maps;

import maps.Projection.Property;
import maps.Projection.Type;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import static java.lang.Double.NaN;
import static java.lang.Double.isFinite;
import static java.lang.Double.isNaN;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.hypot;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.toRadians;

/**
 * A class for completely arbitrary projections, where every square degree can be specified anywhere on the plane.
 *
 * @author Justin Kunimune
 */
public class Elastik {

	public static final ElastikProjection ELASTIK_I = new ElastikProjection(
			"Elastic Earth I", "A map optimised to display landmasses accurately and without interruption.",
			true, Type.OTHER, Property.COMPROMISE, true,
			"elastic-earth-I.csv");


	public static final ElastikProjection ELASTIK_II = new ElastikProjection(
			"Elastic Earth II", "A map optimised to display oceans and their drainage basins accurately and without interruption.",
			true, Type.OTHER, Property.COMPROMISE, true,
			"elastic-earth-II.csv");


	public static final ElastikProjection ELASTIK_III = new ElastikProjection(
			"Elastic Earth III", "A map optimised to show off the continents by compressing the oceans.",
			false, Type.OTHER, Property.COMPROMISE, true,
			"elastic-earth-III.csv");



	private static class ElastikProjection extends Projection {

		private final String filename; // the data filename
		private Polygon[] section_borders; // the unprojected bounds of each section
		private SplineSurface[][] sections; // the x and y projection information for each section
		private double[][][] inverse_raster; // the pixel values, for inverse mapping
		private double raster_left, raster_lower; // the extreme coordinates of the inverse raster
		private double raster_width, raster_height; // the spacial extent of the inverse raster


		public ElastikProjection(
				String title, String description, boolean interrupted, Type type, Property property,
				boolean based_on_land, String filename) {
			super(title, description, 0, 0,
			      true, false, true, !interrupted,
			      type, property, 4,
			      new String[0], new double[0][], !based_on_land);
			this.filename = filename;
		}


		/**
		 * convert a location on the globe to a location on the map plane
		 * @param ф the latitude, in radians
		 * @param λ the longitude, in radians
		 * @return the x value and y value, in the same units as this.width and this.height
		 */
		public double[] project(double ф, double λ) {
			for (int i = 0; i < sections.length; i ++) {
				if (section_borders[i].contains(ф, λ)) {
					double[] result = new double[2];
					for (int l = 0; l < 2; l ++)
						result[l] = sections[i][l].evaluate(ф, λ);
					return result;
				}
			}
			return new double[] {NaN, NaN};
		}


		/**
		 * convert a location on the map plane to a location on the globe
		 * @param x the x value, in the same units as this.width
		 * @param y the y value, in the same units as this.height
		 * @return the latitude and longitude, in radians, or null if the point is not on the mesh
		 */
		public double[] inverse(double x, double y) {
			// start by interpolating on the raster
			double[] guess = inverse_by_interpolation(x, y);
			// then use Newton-Raphson iteration to arrive at an exact solution for each section
			double[][] exact_solutions = new double[sections.length][];
			for (int i = 0; i < sections.length; i ++) {
				exact_solutions[i] = inverse_by_iteration(x, y, i, guess);
				// if any solution is non-null and also in-bounds, use it immediately
				if (exact_solutions[i] != null &&
				    section_borders[i].contains(exact_solutions[i][0], exact_solutions[i][1]))
					return exact_solutions[i];
			}
			// otherwise, arbitrarily choose one of the non-null solutions
			for (int i = 0; i < sections.length; i ++)
				if (exact_solutions[i] != null)  // add 2π to the longitude to mark it as out-of-bounds
					return new double[] {exact_solutions[i][0], exact_solutions[i][1] + 2*PI};
			// if no solutions are non-null, shikatanai.
			return null;
		}

		/**
		 * convert a location on the map plane to an approximate location on the globe by
		 * bilinearly interpolating on the inverse raster
		 * @param x the x value, in the same units as this.width
		 * @param y the y value, in the same units as this.height
		 * @return the latitude and longitude, in radians
		 */
		private double[] inverse_by_interpolation(double x, double y) {
			// find the correct bin
			double i_partial = (y - raster_lower)/raster_height*(inverse_raster.length - 1);
			int i = (int) floor(min(i_partial, inverse_raster.length - 2));
			double j_partial = (x - raster_left)/raster_width*(inverse_raster[i].length - 1);
			int j = (int) floor(min(j_partial, inverse_raster[i].length - 2));

			// calculate the linear weights
			double right_weight = j_partial - j;
			double left_weight = 1 - right_weight;
			double upper_weight = i_partial - i;
			double lower_weight = 1 - upper_weight;

			// perform the interpolation on our 3D cartesian inverse raster
			double[] result = new double[3];
			for (int l = 0; l < 3; l++) {
				result[l] = left_weight*(lower_weight*inverse_raster[i][j][l] +
				                         upper_weight*inverse_raster[i + 1][j][l]) +
				            right_weight*(lower_weight*inverse_raster[i][j + 1][l] +
				                          upper_weight*inverse_raster[i + 1][j + 1][l]);
			}

			// convert to spherical coordinates
			return new double[] {
					atan(result[2]/hypot(result[0], result[1])),
					atan2(result[1], result[0])};
		}

		/**
		 * convert a location on the map plane to a location on the globe by iteratively searching
		 * for the coordinates that project to x and y on a given section.  the algorithm is
		 * Levenberg-Marquardt; see <i>J. Res. NIST</i> 103, p. 633 (1998).
		 * @param x the x value, in the same units as this.width
		 * @param y the y value, in the same units as this.height
		 * @param i the index of the section to search
		 * @param initial_guess an initial input that projects to the correct vicinity
		 * @return the latitude and longitude, in radians, or null if no solution can be found
		 */
		private double[] inverse_by_iteration(double x, double y, int i, double[] initial_guess) {
			final double finite_difference = 1e-5; // radians
			final double second_finite_difference = 0.1; // dimensionless
			final double distance_tolerance = pow(1e-2, 2); // km^2
			final double cosine_tolerance = 1e-3; // dimensionless
			final double backstep_factor = 2.0;
			final double backstep_relaxation_factor = 12.0;
			final int max_num_steps = 40;
			final int max_num_step_sizes = 40;
			double[] target = {x, y};

			// instantiate the state variables
			double[] guess = initial_guess.clone();
			double[] residual = new double[2];
			for (int l = 0; l < 2; l ++)
				residual[l] = sections[i][l].evaluate(guess[0], guess[1]) - target[l];
			double distance = 1/2.*LinAlg.square(residual);
			if (isNaN(distance))  // if this section doesn't cover the initial guess, it's probably pointless
				return null;

			// until we find the solution...
			int num_steps = 0;
			double step_limiter = 0;
			while (true) {

				// check the stopping conditions
				if (distance < distance_tolerance)
					return guess;  // solution is found
				if (num_steps > max_num_steps)
					return null;  // too many iterations have elapsed

				// evaluate the jacobian using finite differences
				double[][] jacobian = new double[2][2];
				for (int l = 0; l < 2; l ++) {
					double north_residual = sections[i][l].evaluate(
							guess[0] + finite_difference, guess[1]) - target[l];
					jacobian[l][0] = (north_residual - residual[l])/finite_difference;
					double east_residual = sections[i][l].evaluate(
							guess[0], guess[1] + finite_difference) - target[l];
					jacobian[l][1] = (east_residual - residual[l])/finite_difference;
					if (isNaN(jacobian[l][0]) || isNaN(jacobian[l][1]))
						return null;  // if we're butting up against the edge, that means there's probably no solution
				}

				double[][] jacobian_transpose = LinAlg.transpose(jacobian);
				double[][] hessian = LinAlg.dot(jacobian_transpose, jacobian);
				double[] gradient = LinAlg.dot(jacobian_transpose, residual);

				// check the other stopping condition
				double sum_cosines_squared = 0;
				for (int l = 0; l < 2; l ++)
					sum_cosines_squared += pow(gradient[l], 2)/LinAlg.square(jacobian_transpose[l])*distance;
				if (sum_cosines_squared < cosine_tolerance)
					return guess;  // local minimum is found

				// perform a backtracking line-search
				double λ_factor = max(.01, pow(cos(guess[0]), 2));
				int num_step_sizes = 0;
				while (true) {

					// calculate the step
					double[][] amplified_hessian = new double[][] {
							{ hessian[0][0] + step_limiter, hessian[0][1] },
							{ hessian[1][0], hessian[1][1] + step_limiter*λ_factor } };
					double[] step = LinAlg.solve_linear_equation(amplified_hessian, gradient);

					// calculate the second-order correction to the step
					double[] expected_change = LinAlg.dot(jacobian, step);
					double[] curvature = new double[2];
					for (int l = 0; l < 2; l ++) {
						double fore_residual = sections[i][l].evaluate(
								guess[0] + step[0]*second_finite_difference,
								guess[1] + step[1]*second_finite_difference) - target[l];
						curvature[l] = 2./second_finite_difference*(
								(fore_residual - residual[l])/second_finite_difference
								- expected_change[l]);
					}
					double[] step_correction = LinAlg.solve_linear_equation(
							amplified_hessian,
							LinAlg.dot(jacobian_transpose, curvature));

					// check the correction magnitude condition before applying the correction
					if (LinAlg.square(step_correction)/LinAlg.square(step) < .25) {
						// take the step
						double[] new_guess = new double[2];
						for (int l = 0; l < 2; l++)
							new_guess[l] = guess[l] - step[l] - 1/2.*step_correction[l];

						// put the new location back in bounds
						new_guess[0] = SplineSurface.fix_latitude(new_guess[0]);
						new_guess[1] = SplineSurface.fix_longitude(new_guess[1]);

						// re-evaluate the function
						double[] new_residual = new double[2];
						for (int l = 0; l < 2; l++)
							new_residual[l] = sections[i][l].evaluate(new_guess[0], new_guess[1]) - target[l];
						double new_distance = 1/2.*LinAlg.square(new_residual);

						// check the line-search stopping conditions
						if (!isNaN(new_distance) && new_distance < distance) {
							// valid step is found
							guess = new_guess;
							residual = new_residual;
							distance = new_distance;
							step_limiter /= backstep_relaxation_factor;
							break;
						}
						if (num_step_sizes > max_num_step_sizes) {
							// too many iterations have elapsed
							return null;
						}
					}

					if (step_limiter == 0)
						step_limiter += (sqrt(1.4*(1 - pow(hessian[0][1], 2)/(hessian[0][0]*hessian[1][1])) + 1) - 1)
						                *min(hessian[0][0], hessian[1][1]/λ_factor);
					else
						step_limiter *= backstep_factor;
					num_step_sizes ++;
				}

				num_steps ++;
			}
		}


		/**
		 * load the data file. this is only called if and when the user tries to use this Projection, so the data file
		 * won’t be loaded until they’re needed.
		 * @param params an ignored parameter
		 * @throws IllegalArgumentException if the data files can’t be found or are in the wrong format
		 */
		@Override
		public void initialize(double... params) throws IllegalArgumentException {
			// start by checking if we’ve already done all this (it only needs to be done once)
			if (sections != null)
				return;

			BufferedReader in = null;
			try {
				in = new BufferedReader(new FileReader(String.format("data/%s", filename))); // parsing the input mesh is pretty simple

				// load the basic projection information
				String line = in.readLine();  // read the header
				int num_sections = parseInt(line.substring(line.length() - 12, line.length() - 11));  // get the number of sections
				section_borders = new Polygon[num_sections];
				sections = new SplineSurface[num_sections][2];

				// load each section
				for (int i = 0; i < num_sections; i ++) {
					in.readLine();  // read the section header (we don’t need any information from it)
					line = in.readLine();  // read the section border header
					int num_vertices = parseInt(line.substring(8, line.length() - 11));  // get the length of the section border
					double[][] border = new double[2][num_vertices];
					for (int j = 0; j < num_vertices; j ++) {
						String[] row = in.readLine().split(",\\s*");
						border[0][j] = toRadians(parseDouble(row[0]));
						border[1][j] = toRadians(parseDouble(row[1]));
					}
					section_borders[i] = new Polygon(border[0], border[1]);
					line = in.readLine();  // read this section points header
					String[] row = line.substring(18, line.length() - 9).split("x");  // get the size of the point grid
					int num_фs = parseInt(row[0]);
					int num_λs = parseInt(row[1]);
					double[][][] points = new double[2][num_фs][num_λs];
					for (int j = 0; j < num_фs; j ++) {
						line = in.readLine();  // read this row of coordinates
						row = line.split(",\\s*");
						for (int k = 0; k < num_λs; k ++) {
							points[0][j][k] = parseDouble(row[2*k]);
							points[1][j][k] = parseDouble(row[2*k + 1]);
						}
					}
					for (int l = 0; l < 2; l ++)
						sections[i][l] = new SplineSurface(points[l]);
				}

				// load the projected border
				line = in.readLine();  // read the projected border header
				int num_vertices = parseInt(line.substring(18, line.length() - 11));  // get the length of the border
				double left = 0, right = 0;
				double lower = 0, upper = 0;
				for (int j = 0; j < num_vertices; j ++) {
					line = in.readLine();  // read the border vertex coordinates
					String[] row = line.split(",\\s*");
					double x = parseDouble(row[0]);
					double y = parseDouble(row[1]);
					left = min(left, x); // use this to update the bounding box
					right = max(right, x);
					lower = min(lower, y);
					upper = max(upper, y);
				}

				// load the inverse raster
				line = in.readLine();  // read the inverse raster header
				String[] row = line.substring(26, line.length() - 9).split("x");  // get the size of the raster
				int num_xs = parseInt(row[0]);
				int num_ys = parseInt(row[1]);
				row = in.readLine().split(",\\s*");  // read the bounding box
				raster_left = parseDouble(row[0]);
				raster_width = parseDouble(row[1]) - raster_left;
				raster_lower = parseDouble(row[2]);
				raster_height = parseDouble(row[3]) - raster_lower;
				inverse_raster = new double[num_ys][num_xs][3];
				for (int j = 0; j < num_ys; j ++) {
					line = in.readLine();  // read each row of coordinates
					row = line.split(",\\s*");
					for (int k = 0; k < num_xs; k ++) {
						double ф = toRadians(parseDouble(row[2*k]));
						double λ = toRadians(parseDouble(row[2*k + 1]));
						inverse_raster[j][k][0] = cos(ф)*cos(λ);  // save the inverse raster in cartesian
						inverse_raster[j][k][1] = cos(ф)*sin(λ);
						inverse_raster[j][k][2] = sin(ф);
					}
				}

				width = right - left;
				height = upper - lower;
			}
			catch (IOException | NullPointerException | ArrayIndexOutOfBoundsException | StringIndexOutOfBoundsException | NumberFormatException e) {
				sections = null;
				inverse_raster = null;
				raster_left = 0.;
				raster_lower = 0.;
				raster_width = 0.;
				raster_height = 0.;
				e.printStackTrace();
				throw new IllegalArgumentException("Missing or corrupt data file for " + this.getName());
			}
			finally {
				try {
					if (in != null) in.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}

			// adjust the splines so that they're fully continuus over the shared nodes
			for (int i_A = 0; i_A < sections.length; i_A ++) {
				// look at each adjacent pair of sections
				int i_B = (i_A + 1)%sections.length;
				// look at each node
				for (int j = 0; j < sections[i_A][0].values.length; j ++) {
					for (int k = 0; k < sections[i_A][0].values[j].length; k ++) {
						// if the nodes are real and in the same location
						if (!isNaN(sections[i_A][0].values[j][k]) &&
						    sections[i_A][0].values[j][k] == sections[i_B][0].values[j][k] &&
						    sections[i_A][1].values[j][k] == sections[i_B][1].values[j][k]) {
							// set them to have the same gradients as well
							SplineSurface.set_gradients_equal(sections, i_A, j, k, i_B, j, k);
						}
					}
				}
			}

			// adjust the splines so they're fully continuus over the antimeridian
			for (int i = 0; i < sections.length; i ++) {
				// look at each latitude
				for (int j = 0; j < sections[i][0].values.length; j ++) {
					// if the nodes are real, they must be in the same location
					int n = sections[i][0].values[j].length;
					if (!isNaN(sections[i][0].values[j][0]) && !isNaN(sections[i][0].values[j][n - 1])) {
						// so set the gradients on the far left and far right to be equal
						SplineSurface.set_gradients_equal(sections, i, j, 0, i, j, n - 1);
					}
				}
			}
		}
	}


	/**
	 * a polygon on a globe with a function for determining if a point is contained by it or not.
	 * the points are assumed to go counterclockwise, so if a point sees the vertices going left,
	 * that point is contained.
	 */
	private static class Polygon {
		private final double[] ф_vertices;
		private final double[] λ_vertices;

		/**
		 * @param ф_vertices the list of vertex latitudes. the last value should be the same as the first one.
		 * @param λ_vertices the list of vertex longitudes. it should be the same length as ф_vertices, and the
		 *                   last value should be the same as the first one.
		 */
		public Polygon(double[] ф_vertices, double[] λ_vertices) {
			this.ф_vertices = ф_vertices;
			this.λ_vertices = λ_vertices;
		}

		/**
		 * determine whether the given point falls within this polygon or not. this is accomplished by drawing a
		 * a north-south line out from the point, finding where it first crosses a polygon edge, and checking if the
		 * edge is going clockwise or widdershins. it is contained iff the nearest edge is going widdershins.
		 * @param ф the latitude of the point
		 * @param λ the longitude of the point
		 * @return true if the point is contained by the polygon and false otherwise
		 */
		public boolean contains(double ф, double λ) {
			int crossings = 0;
			// for each edge of the polygon
			for (int i = 1; i < ф_vertices.length; i ++) {
				if (abs(λ_vertices[i] - λ_vertices[i - 1]) > PI)
					continue;  // skip edges that are wrapping around the backside

				// see if our north-south line crosses it
				boolean crosses = (λ_vertices[i - 1] < λ) != (λ_vertices[i] < λ) ||
				                  (λ == -PI && min(λ_vertices[i - 1], λ_vertices[i]) == -PI);
				if (crosses) {
					// calculate *where* it crosses
					double ф_intersect;
					if (ф_vertices[i - 1] != ф_vertices[i]) {
						double Δλ0 = λ - λ_vertices[i - 1];
						double Δλ1 = λ_vertices[i] - λ;
						ф_intersect = (ф_vertices[i - 1]*Δλ1 + ф_vertices[i]*Δλ0)/(Δλ0 + Δλ1);
					}
					else {
						ф_intersect = ф_vertices[i];  // be mindful of efficiency and avoiding roundoff
					}
					double Δф = abs(ф_intersect - ф);
					if (Δф == 0)
						return true;  // if it's on the line, count it as in
					// otherwise, count this crossing
					if ((ф_intersect < ф) == (λ_vertices[i - 1] < λ_vertices[i]))
						crossings += 1;
					else
						crossings -= 1;
				}
			}
			return crossings > 0;
		}
	}


	/**
	 * the coefficients and functions needed for two-dimensional nan-aware spline interpolation
	 */
	private static class SplineSurface {
		private final double[][] values; // the value of the spline at each node
		private final double[][] gradients_dф; // the derivative of the spline along index 0 at each node
		private final double[][] gradients_dλ; // the derivative of the spline along index 1 at each node

		/**
		 * @param values the value of the function at each node, where the nodes are assumed to span -π/2 <= ф <= π/2
		 *               along the zeroth axis and -π <= λ <= π on the first axis
		 */
		public SplineSurface(double[][] values) {
			this.values = values;

			// set the gradients at the nodes to finite-difference estimates
			final int m = values.length;
			final int n = values[0].length;
			this.gradients_dф = new double[m][n];
			this.gradients_dλ = new double[m][n];

			for (int i = 0; i < m; i ++) {
				for (int j = 0; j < n; j ++) {
					if (isFinite(values[i][j])) {

						if (i - 1 >= 0 && isFinite(values[i - 1][j])) {
							if (i + 1 < m && isFinite(values[i + 1][j]))
								gradients_dф[i][j] = (values[i + 1][j] - values[i - 1][j])/2.;
							else {
								if (i - 2 >= 0 && isFinite(values[i - 2][j]))
									gradients_dф[i][j] = (3*values[i][j] - 4*values[i - 1][j] + values[i - 2][j])/2.;
								else
									gradients_dф[i][j] = values[i][j] - values[i - 1][j];
							}
						}
						else {
							if (i + 1 < m && isFinite(values[i + 1][j])) {
								if (i + 2 < m && isFinite(values[i + 2][j]))
									gradients_dф[i][j] = (-3*values[i][j] + 4*values[i + 1][j] - values[i + 2][j])/2.;
								else
									gradients_dф[i][j] = values[i + 1][j] - values[i][j];
							}
							else
								gradients_dф[i][j] = NaN;
						}

						if (j - 1 >= 0 && isFinite(values[i][j - 1])) {
							if (j + 1 < n && isFinite(values[i][j + 1]))
								gradients_dλ[i][j] = (values[i][j + 1] - values[i][j - 1])/2.;
							else {
								if (j - 2 >= 0 && isFinite(values[i][j - 2]))
									gradients_dλ[i][j] = (3*values[i][j] - 4*values[i][j - 1] + values[i][j - 2])/2.;
								else
									gradients_dλ[i][j] = values[i][j] - values[i][j - 1];
							}
						}
						else {
							if (j + 1 < n && isFinite(values[i][j + 1])) {
								if (j + 2 < n && isFinite(values[i][j + 2]))
									gradients_dλ[i][j] = (-3*values[i][j] + 4*values[i][j + 1] - values[i][j + 2])/2.;
								else
									gradients_dλ[i][j] = values[i][j + 1] - values[i][j];
							}
							else
								gradients_dλ[i][j] = NaN;
						}

					}
					else {
						gradients_dф[i][j] = NaN;
						gradients_dλ[i][j] = NaN;
					}
				}
			}
		}

		/**
		 * interpolate this function to the given coordinates
		 * @param ф the input that moves along the zeroth axis of values
		 * @param λ the input that moves along the first axis of values
		 * @return the value of the spline
		 */
		public double evaluate(double ф, double λ) {
			ф = fix_latitude(ф);
			λ = fix_longitude(λ);
			double i_partial = (ф + PI/2)/PI*(values.length - 1);
			int i = (int)floor(min(i_partial, values.length - 2));
			double j_partial = (λ + PI)/(2*PI)*(values[i].length - 1);
			int j = (int)floor(min(j_partial, values[i].length - 2));
			if (isNaN(values[i][j]) || isNaN(values[i][j + 1]) ||
			    isNaN(values[i + 1][j]) || isNaN(values[i + 1][j + 1]))
				return NaN;

			double z_west = hermite_spline_1D(values[i][j], gradients_dф[i][j],
			                                  values[i + 1][j], gradients_dф[i + 1][j],
			                                  i_partial - i);
			double dzdλ_west = hermite_spline_1D(gradients_dλ[i][j], 0.,
			                                     gradients_dλ[i + 1][j], 0.,
			                                     i_partial - i);
			double z_east = hermite_spline_1D(values[i][j + 1], gradients_dф[i][j + 1],
			                                  values[i + 1][j + 1], gradients_dф[i + 1][j + 1],
			                                  i_partial - i);
			double dzdλ_east = hermite_spline_1D(gradients_dλ[i][j + 1], 0.,
			                                     gradients_dλ[i + 1][j + 1], 0.,
			                                     i_partial - i);
			return hermite_spline_1D(z_west, dzdλ_west, z_east, dzdλ_east,
			                         j_partial - j);
		}

		/**
		 * a polynomial with specified values and derivatives at two endpoints
		 * @param y_0 the value of the function at x=0
		 * @param dydx_0 the slope of the function at x=0
		 * @param y_1 the value of the function at x=1
		 * @param dydx_1 the slope of the function at x=1
		 * @param x the dimensionless function input in the range [0, 1]
		 * @return a smoothly varying value somewhere in the vicinity of y_0 and y_1
		 */
		private static double hermite_spline_1D(double y_0, double dydx_0, double y_1, double dydx_1, double x) {
			return y_0 +
			       x*(dydx_0 +
			          x*((3*(y_1 - y_0) - 2*dydx_0 - dydx_1) +
			             x*(dydx_0 + dydx_1 - 2*(y_1 - y_0))));
		}

		/**
		 * average two gradients at two locations in a stack of spline surfaces so that the splines are
		 * C^1 continuous at that point.
		 */
		private static void set_gradients_equal(
				SplineSurface[][] surfaces, int i_A, int j_A, int k_A, int i_B, int j_B, int k_B) {
			for (int l = 0; l < 2; l ++) {
				double mean_gradient_dф = (
						                          surfaces[i_A][l].gradients_dф[j_A][k_A] +
						                          surfaces[i_B][l].gradients_dф[j_B][k_B])/2;
				surfaces[i_A][l].gradients_dф[j_A][k_A] = mean_gradient_dф;
				surfaces[i_B][l].gradients_dф[j_B][k_B] = mean_gradient_dф;
				double mean_gradient_dλ = (
						                          surfaces[i_A][l].gradients_dλ[j_A][k_A] +
						                          surfaces[i_B][l].gradients_dλ[j_B][k_B])/2;
				surfaces[i_A][l].gradients_dλ[j_A][k_A] = mean_gradient_dλ;
				surfaces[i_B][l].gradients_dλ[j_B][k_B] = mean_gradient_dλ;
			}
		}

		/**
		 * recenter a latitude so that it is in the range [-PI/2, PI/2].  if it's outside of that,
		 * it will get reflected about the poles until it isn't.
		 * @param ф the latitude (radians)
		 * @return the equivalent in-bounds latitude
		 */
		private static double fix_latitude(double ф) {
			if (ф < -PI/2 || ф > PI/2) {
				ф = ф - floor((ф + PI/2)/(2*PI))*(2*PI);
				if (ф > PI/2)
					ф = PI - ф;
			}
			return ф;
		}

		/**
		 * recenter a longitude so that it is in the range [-PI, PI].  if it's outside of that, it
		 * will get modulused such that it isn't.
		 * @param λ the longitude (radians)
		 * @return the equivalent in-bounds longitude
		 */
		private static double fix_longitude(double λ) {
			if (λ < -PI || λ > PI)
				λ = λ - floor((λ + PI)/(2*PI))*(2*PI);
			return λ;
		}
	}


	/**
	 * some static matrix operations for the Levenberg-Marquardt implementation
	 */
	private static class LinAlg {
		/**
		 * evaluate the inner product of two vectors
		 * @return the 2-vector resulting from the dot-product of A and b
		 */
		private static double dot(double[] a, double[] b) {
			if (a.length != b.length)
				throw new IllegalArgumentException("the vectors must have the same length.");
			double product = 0;
			for (int i = 0; i < a.length; i++)
				product += a[i]*b[i];
			return product;
		}

		/**
		 * evaluate the linear equation Ab where A is a matrix and b is a column vector
		 * @return the 2-vector resulting from the dot-product of A and b
		 */
		private static double[] dot(double[][] A, double[] b) {
			if (A.length > 0 && A[0].length != b.length)
				throw new IllegalArgumentException("the number of columns in A must be the same as the length of b.");
			double[] product = new double[A.length];
			for (int i = 0; i < A.length; i++)
				for (int j = 0; j < A[i].length; j++)
					product[i] += A[i][j]*b[j];
			return product;
		}

		/**
		 * evaluate the matrix product of two matrices
		 * @return the 2-vector resulting from the dot-product of A and b
		 */
		private static double[][] dot(double[][] A, double[][] B) {
			if (A.length > 0 && A[0].length != B.length)
				throw new IllegalArgumentException("this function has only been implemented for 2D equations.");
			double[][] product = new double[A.length][B[0].length];
			for (int i = 0; i < A.length; i++)
				for (int j = 0; j < A[i].length; j++)
					for (int k = 0; k < B[j].length; k++)
						product[i][k] += A[i][j]*B[j][k];
			return product;
		}

		/**
		 * evaluate the square of the L-2 norm of a vector
		 * @return the dot-product of a with itself
		 */
		private static double square(double[] a) {
			return dot(a, a);
		}

		/**
		 * flip the rows of A with its columns
		 * @return the transpose of A
		 */
		private static double[][] transpose(double[][] A) {
			double[][] AT = new double[A[0].length][A.length];
			for (int i = 0; i < A.length; i ++)
				for (int j = 0; j < A[i].length; j ++)
					AT[j][i] = A[i][j];
			return AT;
		}

		/**
		 * find the vector x such that Ax = b.  this function only works on 2d systems.
		 * @param A a 2×2 matrix to be inverted and multiplied
		 * @param b a 2-vector to be multiplied by the inverse of A
		 * @return the 2-vector x that solves the linear equation Ax = b
		 */
		private static double[] solve_linear_equation(double[][] A, double[] b) {
			if (A.length != 2 || b.length != 2)
				throw new IllegalArgumentException("this function has only been implemented for 2D equations.");
			double det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
			return new double[]{
					(A[1][1]*b[0] - A[0][1]*b[1])/det,
					(A[0][0]*b[1] - A[1][0]*b[0])/det};
		}
	}
}
