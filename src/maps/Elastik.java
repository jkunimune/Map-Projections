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
import static java.lang.Double.POSITIVE_INFINITY;
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
import static java.lang.Math.min;
import static java.lang.Math.sin;
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
		private double[][] projected_border; // the x and y values of the full projected map outline
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
		 * @param ф the latitude in radians
		 * @param λ the longitude in radians
		 * @return the x value and y value in the same units as this.width and this.height
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
		 * @param x the x value in the same units as this.width
		 * @param y the y value in the same units as this.height
		 * @return the latitude and longitude in radians
		 */
		public double[] inverse(double x, double y) {
			// find the correct bin
			double i_partial = (y - raster_lower)/raster_height*(inverse_raster.length - 1);
			int i = (int)floor(min(i_partial, inverse_raster.length - 2));
			double j_partial = (x - raster_left)/raster_width*(inverse_raster[i].length - 1);
			int j = (int)floor(min(j_partial, inverse_raster[i].length - 2));

			// calculate the linear weights
			double right_weight = j_partial - j;
			double left_weight = 1 - right_weight;
			double upper_weight = i_partial - i;
			double lower_weight = 1 - upper_weight;

			// perform the interpolation on our 3D cartesian inverse raster
			double[] result = new double[3];
			for (int l = 0; l < 3; l ++) {
				result[l] = left_weight*(lower_weight*inverse_raster[i][j][l] +
				                         upper_weight*inverse_raster[i + 1][j][l]) +
				            right_weight*(lower_weight*inverse_raster[i][j + 1][l] +
				                          upper_weight*inverse_raster[i + 1][j + 1][l]);
			}

			// convert to spherical coordinates and return
			return new double[] {
					atan(result[2]/hypot(result[0], result[1])),
			        atan2(result[1], result[0]) };
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
				projected_border = new double[num_vertices][2];
				for (int j = 0; j < projected_border.length; j ++) {
					line = in.readLine();  // read the border vertex coordinates
					String[] row = line.split(",\\s*");
					projected_border[j][0] = parseDouble(row[0]);
					projected_border[j][1] = parseDouble(row[1]);
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

				// calculate width and height
				double left = 0, right = 0;
				double lower = 0, upper = 0;
				for (double[] point : projected_border) {
					if (point[0] < left)
						left = point[0];
					if (point[0] > right)
						right = point[0];
					if (point[1] < lower)
						lower = point[1];
					if (point[1] > upper)
						upper = point[1];
				}
				width = right - left;
				height = upper - lower;
			}
			catch (IOException | NullPointerException | ArrayIndexOutOfBoundsException | StringIndexOutOfBoundsException | NumberFormatException e) {
				sections = null;
				projected_border = null;
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
							for (int l = 0; l < 2; l ++) {
								double mean_gradient_dф = (
										sections[i_A][l].gradients_dф[j][k] +
										sections[i_B][l].gradients_dф[j][k])/2;
								sections[i_A][l].gradients_dф[j][k] = mean_gradient_dф;
								sections[i_B][l].gradients_dф[j][k] = mean_gradient_dф;
								double mean_gradient_dλ = (
										sections[i_A][l].gradients_dλ[j][k] +
										sections[i_B][l].gradients_dλ[j][k])/2;
								sections[i_A][l].gradients_dλ[j][k] = mean_gradient_dλ;
								sections[i_B][l].gradients_dλ[j][k] = mean_gradient_dλ;
							}
						}
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
		 * @param λ_vertices the list of vertex longitudes. it should be the same length as xVertices, and the
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
			double nearestSegment = POSITIVE_INFINITY;
			boolean contained = false;
			// for each edge of the polygon
			for (int i = 1; i < ф_vertices.length; i ++) {
				if (abs(λ_vertices[i] - λ_vertices[i - 1]) > PI)
					continue;  // skip edges that are wrapping around the backside
				// if our north-south line crosses it
				if (λ_vertices[i - 1] != λ_vertices[i]) {
					if ((λ_vertices[i - 1] < λ) != (λ_vertices[i] < λ) ||
					    λ_vertices[i - 1] == λ || λ_vertices[i] == λ) {
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
						// if this is closer than any intersections we've found before
						double Δф = abs(ф_intersect - ф);
						if (Δф < nearestSegment) {
							// record it
							nearestSegment = Δф;
							// and save the polarity
							if (Δф != 0)
								contained = (ф_intersect < ф) == (λ_vertices[i - 1] < λ_vertices[i]);
							else
								contained = true;  // if it's on the line, count it as in
						}
					}
				}
			}
			return contained;
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

			// set the gradients at the nodes
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

		public double evaluate(double ф, double λ) {
			double i_partial = (ф + PI/2)/PI*(values.length - 1);
			int i = (int)floor(min(i_partial, values.length - 2));
			double j_partial = (λ + PI)/(2*PI)*(values[i].length - 1);
			int j = (int)floor(min(j_partial, values[i].length - 2));
			if (isNaN(values[i][j]) || isNaN(values[i][j + 1]) ||
			    isNaN(values[i + 1][j]) || isNaN(values[i + 1][j + 1]))
				return NaN;

			double z_west = hermiteSpline1D(values[i][j], gradients_dф[i][j],
			                                values[i + 1][j], gradients_dф[i + 1][j],
			                                i_partial - i);
			double dzdλ_west = hermiteSpline1D(gradients_dλ[i][j], 0.,
			                                   gradients_dλ[i + 1][j], 0.,
			                                   i_partial - i);
			double z_east = hermiteSpline1D(values[i][j + 1], gradients_dф[i][j + 1],
			                                values[i + 1][j + 1], gradients_dф[i + 1][j + 1],
			                                i_partial - i);
			double dzdλ_east = hermiteSpline1D(gradients_dλ[i][j + 1], 0.,
			                                   gradients_dλ[i + 1][j + 1], 0.,
			                                   i_partial - i);
			return hermiteSpline1D(z_west, dzdλ_west, z_east, dzdλ_east,
			                       j_partial - j);
		}

		private static double hermiteSpline1D(double y_0, double dydx_0, double y_1, double dydx_1, double x) {
			return y_0 +
			       x*(dydx_0 +
			          x*((3*(y_1 - y_0) - 2*dydx_0 - dydx_1) +
			             x*(dydx_0 + dydx_1 - 2*(y_1 - y_0))));
		}
	}
}
