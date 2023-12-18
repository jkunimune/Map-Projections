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
import utils.Shape;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import static java.lang.Double.NaN;
import static java.lang.Double.isFinite;
import static java.lang.Double.isNaN;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.Math.PI;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.hypot;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.toRadians;
import static java.lang.String.format;

/**
 * A class for completely arbitrary projections, where every square degree can be specified anywhere on the plane.
 *
 * @author Justin Kunimune
 */
public class Elastic {

	public static final ElasticProjection ELASTIC_I = new ElasticProjection(
			"Elastic I", "A map optimised to display landmasses accurately and without interruption",
			true, Type.OTHER, Property.COMPROMISE, true,
			"elastic-I.txt");


	public static final ElasticProjection ELASTIC_II = new ElasticProjection(
			"Elastic II", "A map optimised to display oceans and their drainage basins accurately and without interruption",
			true, Type.OTHER, Property.COMPROMISE, true,
			"elastic-II.txt");


	public static final ElasticProjection ELASTIC_III = new ElasticProjection(
			"Elastic III", "A map optimised to show off the continents by compressing the oceans",
			false, Type.OTHER, Property.COMPROMISE, true,
			"elastic-III.txt");



	private static class ElasticProjection extends Projection {

		private final String filename; // the data filename
		private Region[] section_borders; // the unprojected bounds of each section
		private SplineSurface[][] sections; // the x and y projection information for each section
		private double[][][] inverse_raster; // the pixel values, for inverse mapping
		private double raster_left, raster_lower; // the extreme coordinates of the inverse raster
		private double raster_width, raster_height; // the spacial extent of the inverse raster


		public ElasticProjection(
				String title, String description, boolean interrupted, Type type, Property property,
				boolean based_on_land, String filename) {
			super(title, description, null,
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
		@Override
		public double[] project(double ф, double λ) {
			for (int i = 0; i < sections.length; i ++) {
				// find the section that contains this point
				if (section_borders[i].contains(ф, λ)) {
					// and project using that section
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
		@Override
		public double[] inverse(double x, double y) {
			// start by interpolating on the raster
			double[] guess = inverse_by_interpolation(x, y);
			if (guess == null)
				return null;

			// then use Newton-Raphson iteration to arrive at an exact solution for each section
			double[][] exact_solutions = new double[sections.length][];
			for (int section_index = 0; section_index < sections.length; section_index ++) {
				exact_solutions[section_index] = inverse_by_iteration(x, y, section_index, guess);
				// if any solution is non-null and also in-bounds, use it immediately
				if (exact_solutions[section_index] != null &&
				    section_borders[section_index].contains(
							exact_solutions[section_index][0], exact_solutions[section_index][1]))
					return exact_solutions[section_index];
			}

			// otherwise, arbitrarily choose one of the non-null solutions
			for (int section_index = 0; section_index < sections.length; section_index ++)
				if (exact_solutions[section_index] != null)  // add 2π to the longitude to mark it as out-of-bounds
					return new double[] {exact_solutions[section_index][0],
					                     exact_solutions[section_index][1] + 2*PI};
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
			double i_partial = (x - raster_left)/raster_width*(inverse_raster.length - 1);
			int i = (int) floor(min(i_partial, inverse_raster.length - 2));
			double j_partial = (y - raster_lower)/raster_height*(inverse_raster[i].length - 1);
			int j = (int) floor(min(j_partial, inverse_raster[i].length - 2));
			
			// check that we have everything we need to interpolate
			if (isNaN(inverse_raster[i][j][0]) || isNaN(inverse_raster[i][j + 1][0]) ||
			    isNaN(inverse_raster[i + 1][j][0]) || isNaN(inverse_raster[i + 1][j + 1][0]))
				return null;

			// calculate the linear weights
			double right_weight = i_partial - i;
			double left_weight = 1 - right_weight;
			double upper_weight = j_partial - j;
			double lower_weight = 1 - upper_weight;

			// perform the interpolation on our 3D cartesian inverse raster
			double[] result = new double[3];
			for (int l = 0; l < 3; l ++) {
				result[l] = left_weight*(lower_weight*inverse_raster[i][j][l] +
				                         upper_weight*inverse_raster[i][j + 1][l]) +
				            right_weight*(lower_weight*inverse_raster[i + 1][j][l] +
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
		 * Levenberg-Marquardt; see <i>J. Res. NIST</i> 103, p. 633 (1998).  this procedure could
		 * probably be a lot simpler (or omitted entirely) depending on how high quality an inverse
		 * you want.
		 * @param x the x value, in the same units as this.width
		 * @param y the y value, in the same units as this.height
		 * @param section_index the index of the section to search
		 * @param initial_guess an initial input that projects to the correct vicinity
		 * @return the latitude and longitude, in radians, or null if no solution can be found
		 */
		private double[] inverse_by_iteration(double x, double y, int section_index, double[] initial_guess) {
			final double finite_difference = 1e-7; // radians
			final double second_finite_difference = 0.1; // dimensionless
			final double strict_cost_tolerance = pow(.1, 2); // km^2
			final double lenient_cost_tolerance = pow(1000, 2); // km^2
			final double cosine_tolerance = 1e-2; // dimensionless
			final double backstep_factor = 2.0;
			final double backstep_relaxation_factor = 12.0;
			final int max_num_steps = 40;
			final int max_num_step_sizes = 40;
			double[] target = {x, y};
			
			SplineSurface[] section = sections[section_index];

			// instantiate the state variables
			double[] guess = initial_guess.clone();
			double[] residual = new double[2];
			for (int l = 0; l < 2; l ++)
				residual[l] = section[l].evaluate(guess[0], guess[1]) - target[l];
			double cost = 1/2.*LinAlg.square(residual);
			if (isNaN(cost))  // if this section doesn't cover the initial guess, it's probably pointless
				return null;

			// until we find the solution...
			int num_steps = 0;
			double step_limiter = 0;
			while (true) {

				// check the stopping conditions
				if (cost < strict_cost_tolerance)
					return guess;  // solution is found
				if (num_steps > max_num_steps)
					return null;  // too many outer iterations have elapsed

				// evaluate the jacobian using finite differences
				double[][] jacobian = new double[2][2];
				for (int l = 0; l < 2; l ++) {
					double north_residual = section[l].evaluate(
							guess[0] + finite_difference, guess[1]) - target[l];
					jacobian[l][0] = (north_residual - residual[l])/finite_difference;
					double east_residual = section[l].evaluate(
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
					sum_cosines_squared += pow(gradient[l], 2)/(LinAlg.square(jacobian_transpose[l])*cost);
				if (sum_cosines_squared < cosine_tolerance && cost < lenient_cost_tolerance)
					return guess;  // local minimum is found

				// perform a backtracking line-search
				double λ_factor = max(.01, pow(cos(guess[0]), 2));
				int num_step_sizes = 0;
				while (true) {

					// calculate the step
					double[][] amplified_hessian = new double[][] {
							{ hessian[0][0] + step_limiter, hessian[0][1] },
							{ hessian[1][0], hessian[1][1] + step_limiter*λ_factor } };
					double[] step = LinAlg.dot(
							-1, LinAlg.solve_linear_equation(amplified_hessian, gradient));

					// calculate the second-order correction to the step
					double[] expected_slope = LinAlg.dot(jacobian, step);
					double[] curvature = new double[2];
					for (int l = 0; l < 2; l ++) {
						double fore_residual = section[l].evaluate(
								guess[0] + step[0]*second_finite_difference,
								guess[1] + step[1]*second_finite_difference) - target[l];
						curvature[l] = 2./second_finite_difference*(
								(fore_residual - residual[l])/second_finite_difference
								- expected_slope[l]);
					}
					double[] step_correction = LinAlg.dot(
							-1/2., LinAlg.solve_linear_equation(
									amplified_hessian,
									LinAlg.dot(jacobian_transpose, curvature)));

					// check the correction magnitude condition before applying the correction
					if (LinAlg.square(step_correction)/LinAlg.square(step) < .25) {
						// take the step
						double[] new_guess = new double[2];
						for (int l = 0; l < 2; l++)
							new_guess[l] = guess[l] + step[l] + step_correction[l];

						// put the new location back in bounds
						new_guess = SplineSurface.fix_latitude_and_longitude(new_guess[0], new_guess[1]);

						// re-evaluate the function
						double[] new_residual = new double[2];
						for (int l = 0; l < 2; l++)
							new_residual[l] = section[l].evaluate(new_guess[0], new_guess[1]) - target[l];
						double new_cost = 1/2.*LinAlg.square(new_residual);

						// check the line-search stopping conditions
						if (!isNaN(new_cost) && new_cost < cost) {
							// valid step is found
							guess = new_guess;
							residual = new_residual;
							cost = new_cost;
							step_limiter /= backstep_relaxation_factor;
							break;
						}
						if (num_step_sizes > max_num_step_sizes) {
							// too many inner iterations have elapsed
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
		 * won’t be loaded until it's needed.
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
				in = new BufferedReader(new FileReader(format("res/%s", filename))); // parsing the input mesh is pretty simple

				// load the basic projection information
				String line = in.readLine();  // read the header
				int num_sections = parseInt(line.substring(line.length() - 12, line.length() - 11));  // get the number of sections
				section_borders = new Region[num_sections];
				sections = new SplineSurface[num_sections][2];

				// load each section
				for (int i = 0; i < num_sections; i ++) {
					in.readLine();  // read the section header (we don’t need any information from it)
					line = in.readLine();  // read the section border header
					int num_vertices = parseInt(line.substring(10, line.length() - 11));  // get the length of the section border
					double[][] border = new double[2][num_vertices];
					for (int j = 0; j < num_vertices; j ++) {
						String[] row = in.readLine().split(",\\s*");
						border[0][j] = toRadians(parseDouble(row[0]));
						border[1][j] = toRadians(parseDouble(row[1]));
					}
					section_borders[i] = new Region(border[0], border[1]);

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
				int num_vertices = parseInt(line.substring(20, line.length() - 11));  // get the length of the border
				double[][] coordinates = new double[num_vertices][2];
				for (int i = 0; i < num_vertices; i ++) {
					line = in.readLine();  // read the border vertex coordinates
					String[] row = line.split(",\\s*");
					coordinates[i][0] = parseDouble(row[0]);
					coordinates[i][1] = parseDouble(row[1]);
				}
				shape = Shape.polygon(coordinates);  // wrap the vertices in a Shape object

				// load the inverse raster
				line = in.readLine();  // read this section inverse points header
				String[] row = line.substring(26, line.length() - 9).split("x");  // get the size of the point grid
				int num_xs = parseInt(row[0]);
				int num_ys = parseInt(row[1]);
				row = in.readLine().split(",\\s*");  // read the bounding box
				raster_left = parseDouble(row[0]);
				raster_lower = parseDouble(row[1]);
				raster_width = parseDouble(row[2]) - raster_left;
				raster_height = parseDouble(row[3]) - raster_lower;
				inverse_raster = new double[num_xs][num_ys][3];
				for (int i = 0; i < num_xs; i ++) {
					line = in.readLine();  // read each row of coordinates
					row = line.split(",\\s*");
					for (int j = 0; j < num_ys; j ++) {
						double ф = toRadians(parseDouble(row[2*j]));
						double λ = toRadians(parseDouble(row[2*j + 1]));
						inverse_raster[i][j][0] = cos(ф)*cos(λ);  // save the inverse raster in cartesian
						inverse_raster[i][j][1] = cos(ф)*sin(λ);
						inverse_raster[i][j][2] = sin(ф);
					}
				}

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
					int n = sections[i][0].values[j].length;
					// set the gradients on the far left and far right to be equal
					SplineSurface.set_gradients_equal(sections, i, j, 0, i, j, n - 1);
				}
			}
		}
	}


	/**
	 * a polygon on a globe with a function for determining if a point is contained by it or not.
	 * the points are assumed to go counterclockwise, so if a point sees the vertices going left,
	 * that point is contained.
	 */
	private static class Region {
		private final double[] ф_vertices;
		private final double[] λ_vertices;

		/**
		 * @param ф_vertices the list of vertex latitudes. the last value should be the same as the first one.
		 * @param λ_vertices the list of vertex longitudes. it should be the same length as ф_vertices, and the
		 *                   last value should be the same as the first one.
		 */
		public Region(double[] ф_vertices, double[] λ_vertices) {
			this.ф_vertices = ф_vertices;
			this.λ_vertices = λ_vertices;
		}

		/**
		 * determine whether the given point falls within this polygon or not. the method used is a typical even/odd
		 * algorithm, with the caveat that the directionality of the polygon is significant.  if the vertices go
		 * counterclockwise, the contained points are everything inside the polygon as one would expect, but if the
		 * vertices go clockwise, the contained points are actually everything <i>outside</i> the polygon.
		 * @param ф the latitude of the point
		 * @param λ the longitude of the point
		 * @return true if the point is contained by the polygon and false otherwise
		 */
		public boolean contains(double ф, double λ) {
			// this is a slightly modified version of your typical ray casting scheme, but instead of just counting
			// the crossings, we weigh them by their distance.  it's like the polygon is a current-carrying wire;
			// points inside have a positive "magnetic field" and points outside have a negative "magnetic field".
			double magnetic_field = 0;

			// for each edge of the polygon
			for (int i = 1; i < ф_vertices.length; i ++) {
				// that is not due north-south
				boolean vertical = λ_vertices[i - 1] == λ_vertices[i];
				if (!vertical) {
					// see if our north-south line crosses it
					boolean crosses = (λ_vertices[i - 1] < λ) != (λ_vertices[i] < λ) ||
					                  (λ == -PI && min(λ_vertices[i - 1], λ_vertices[i]) == -PI);
					if (crosses) {
						// calculate *where* it crosses
						double ф_intersect;
						if (ф_vertices[i - 1] == ф_vertices[i])
							ф_intersect = ф_vertices[i];  // be mindful of efficiency and avoiding roundoff errors
						else if (λ_vertices[i - 1] == λ)
							ф_intersect = ф_vertices[i - 1];
						else if (λ_vertices[i] == λ)
							ф_intersect = ф_vertices[i];
						else {
							double Δλ0 = λ - λ_vertices[i - 1];
							double Δλ1 = λ_vertices[i] - λ;
							ф_intersect = (ф_vertices[i - 1]*Δλ1 + ф_vertices[i]*Δλ0)/(Δλ0 + Δλ1);
						}
						// if it's on the line, count this point as in
						if (ф_intersect == ф)
							return true;
						// otherwise, check the direction
						else
							magnetic_field += signum(λ_vertices[i - 1] - λ_vertices[i])/
							                  (ф_intersect - ф);
					}
				}
			}
			return magnetic_field > 0;
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
		 * @param values the value of the function at each node, where the nodes are assumed to span
		 *               -π/2 <= ф <= π/2 along the zeroth axis and -π <= λ <= π on the first axis
		 */
		public SplineSurface(double[][] values) {
			// set the values at the nodes to the specification
			this.values = values;
			// set the gradients at the nodes to finite-difference estimates
			double[][][] gradients = estimate_gradients(values);
			this.gradients_dф = gradients[0];
			this.gradients_dλ = gradients[1];
		}

		/**
		 * use finite-differences to estimate what the gradient at each node should
		 * be to ensure the smoothest map projection
		 * @param values the value of the function at each node
		 * @return the approximate latitude-derivative of the function at each node, and
		 *         the approximate longitude-derivate of the function at each node
		 */
		private static double[][][] estimate_gradients(double[][] values) {
			// instantiate both gradient arrays
			final int m = values.length;
			final int n = values[0].length;
			double[][] gradients_dф = new double[m][n];
			double[][] gradients_dλ = new double[m][n];

			// use twoth- or first-order finite differences to estimate gradients
			for (int i = 0; i < m; i ++) {
				for (int j = 0; j < n; j ++) {
					if (isFinite(values[i][j])) {

						// derivative /dф
						if (i - 1 >= 0 && isFinite(values[i - 1][j])) {
							if (i + 1 < m && isFinite(values[i + 1][j]))
								// central second-order
								gradients_dф[i][j] = (values[i + 1][j] - values[i - 1][j])/2.;
							else
								// backward first-order
								gradients_dф[i][j] = values[i][j] - values[i - 1][j];
						}
						else if (i + 1 < m && isFinite(values[i + 1][j]))
							// forward first-order
							gradients_dф[i][j] = values[i + 1][j] - values[i][j];

						// derivative /dλ
						if (j - 1 >= 0 && isFinite(values[i][j - 1])) {
							if (j + 1 < n && isFinite(values[i][j + 1]))
								// central second-order
								gradients_dλ[i][j] = (values[i][j + 1] - values[i][j - 1])/2.;
							else
								// backward first-order
								gradients_dλ[i][j] = values[i][j] - values[i][j - 1];
						}
						else if (j + 1 < n && isFinite(values[i][j + 1]))
							// forward first-order
							gradients_dλ[i][j] = values[i][j + 1] - values[i][j];
					}
				}
			}

			// wrap the two arrays together as a 3D array so they can both be returnd
			return new double[][][] {gradients_dф, gradients_dλ};
		}

		/**
		 * interpolate this function to the given coordinates
		 * @param ф the input that moves along the zeroth axis of values
		 * @param λ the input that moves along the first axis of values
		 * @return the value of the spline
		 */
		public double evaluate(double ф, double λ) {
			double[] fixed_coordinates = fix_latitude_and_longitude(ф, λ);
			ф = fixed_coordinates[0];
			λ = fixed_coordinates[1];

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
		 * recenter a set of spherical coordinates so that the latitude is in the range [-PI/2, PI/2]
		 * and the longitude is in the range [-PI, PI].  if it's outside of the latitude bounds, it
		 * will get reflected about the poles until it isn't.  if it's out of the longitude bounds,
		 * it will get moduloed such that it isn't.
		 * @param ф the latitude (radians)
		 * @param λ the longitude (radians)
		 * @return the equivalent in-bounds latitude and longitude in an array
		 */
		private static double[] fix_latitude_and_longitude(double ф, double λ) {
			if (ф < -PI/2 || ф > PI/2) {
				ф = ф - floor((ф + PI/2)/(2*PI))*(2*PI);
				if (ф > PI/2) {
					ф = PI - ф;
					λ = λ + PI;
				}
			}
			if (λ < -PI || λ > PI)
				λ = λ - floor((λ + PI)/(2*PI))*(2*PI);
			return new double[] {ф, λ};
		}
	}


	/**
	 * some static matrix operations for the Levenberg-Marquardt implementation
	 */
	private static class LinAlg {
		/**
		 * evaluate the product of a scalar with a vector
		 * @return the vector scaled by a factor of a
		 */
		private static double[] dot(double a, double[] b) {
			double[] product = new double[b.length];
			for (int i = 0; i < b.length; i ++)
				product[i] = a*b[i];
			return product;
		}

		/**
		 * evaluate the inner product of two vectors
		 * @return the 2-vector resulting from the dot-product of a and b
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
