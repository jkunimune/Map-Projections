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

import utils.linalg.Matrix;
import utils.linalg.Vector;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.ceil;
import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.floorMod;
import static java.lang.Math.hypot;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sin;
import static utils.Math2.coerceAngle;

/**
 * A class for completely arbitrary projections, specified by an unstructured triangular mesh.
 *
 * @author Justin Kunimune
 */
public class Liquid {

	public static final LiquidProjection LIQUID_EARTH = new LiquidProjection(
			"Liquid Earth", "An equal-area map in the shape of the Equal Earth projection, optimized around the continents",
			EqualEarth.EQUAL_EARTH, true, "liquid");


	private static class LiquidProjection extends Projection {
		
		private final String filename; // the data filename
		private final Projection baseProjection;
		private int[][] triangles;
		private Vector[] verticesInitial;
		private Vector[] verticesTransformed;
		private double[] фBins;
		private double[] λBins;
		private List<int[]>[][] lookupTableInitial;
		private List<int[]>[][] lookupTableTransformed;
		
		
		public LiquidProjection(
				String title, String description, Projection baseProjection,
				boolean based_on_land, String filename) {
			super(title, description, null,
			      true, false, true, baseProjection.isContinuous(),
			      Type.OTHER, baseProjection.getProperty(), 4,
			      new String[0], new double[0][], !based_on_land);
			this.filename = filename;
			this.baseProjection = baseProjection;
		}
		
		
		/**
		 * convert a location on the globe to a location on the map plane
		 *
		 * @param ф the latitude, in radians
		 * @param λ the longitude, in radians
		 * @return the x value and y value, in the same units as this.width and this.height
		 */
		@Override
		public double[] project(double ф, double λ) {
			int[] triangle = findContainingTriangle(
					this.фBins, this.λBins, this.lookupTableInitial, this.verticesInitial, ф, λ);
			int a = triangle[0];
			int b = triangle[1];
			int c = triangle[2];
			double[] transformedCoordinates = mapFromTriangleToTriangle(
					this.verticesInitial[a], this.verticesInitial[b], this.verticesInitial[c],
					this.verticesTransformed[a], this.verticesTransformed[b], this.verticesTransformed[c],
					ф, λ);
			return this.baseProjection.project(transformedCoordinates[0], transformedCoordinates[1]);
		}
		
		
		/**
		 * convert a location on the map plane to a location on the globe
		 *
		 * @param x the x value, in the same units as this.width
		 * @param y the y value, in the same units as this.height
		 * @return the latitude and longitude, in radians, or null if the point is not on the mesh
		 */
		@Override
		public double[] inverse(double x, double y) {
			double[] transformedCoordinates = this.baseProjection.inverse(x, y);
			if (transformedCoordinates == null)
				return null;
			int[] triangle = findContainingTriangle(
					this.фBins, this.λBins, this.lookupTableTransformed,
					this.verticesTransformed,
					transformedCoordinates[0], transformedCoordinates[1]);
			int a = triangle[0];
			int b = triangle[1];
			int c = triangle[2];
			return mapFromTriangleToTriangle(
					this.verticesTransformed[a], this.verticesTransformed[b], this.verticesTransformed[c],
					this.verticesInitial[a], this.verticesInitial[b], this.verticesInitial[c],
					transformedCoordinates[0], transformedCoordinates[1]);
		}
		
		
		private static int[] findContainingTriangle(
				double[] фBins, double[] λBins, List<int[]>[][] lookupTable,
				Vector[] vertices, double ф, double λ) {
			λ = coerceAngle(λ);
			Vector point = new Vector(cos(ф)*cos(λ), cos(ф)*sin(λ), sin(ф));
			int i = (int) floor((ф - фBins[0])/(фBins[1] - фBins[0]));
			int j = (int) floor((λ - λBins[0])/(λBins[1] - λBins[0]));
			triangleLoop:
			for (int[] triangle: lookupTable[i][j]) {
				for (int k = 0; k < 3; k ++) {
					assert vertices[triangle[k]].dot(vertices[triangle[(k + 1)%3]].cross(vertices[triangle[(k + 2)%3]])) > 0;
					double sign = vertices[triangle[k]].dot(vertices[triangle[(k + 1)%3]].cross(point));
					if (sign < 0)
						continue triangleLoop;
				}
				return triangle;
			}
			throw new RuntimeException(String.format(
					"the point <%.3fN, %.3fE> is not contained by any triangle.", ф, λ));
		}
		
		
		/**
		 * map a point from one spherical triangle to another, such that the edges map linearly between each other.
		 * this function assumes that all of the passed vectors have unit length
		 */
		private static double[] mapFromTriangleToTriangle(
				Vector a1, Vector b1, Vector c1,
				Vector a2, Vector b2, Vector c2,
				double ф, double λ) {
			// convert the point to cartesian coordinates
			Vector initialPoint = new Vector(cos(ф)*cos(λ), cos(ф)*sin(λ), sin(ф));
			
			// gnomonicly project the point to the input 3D triangle plane
			Vector normal1 = a1.minus(b1).cross(c1.minus(a1));
			Vector point1 = initialPoint.times(a1.dot(normal1)/initialPoint.dot(normal1)); // the gnomonic projection of the point onto the input triangle's plane
			
			// convert the projected 3D point to triangular coordinates
			double s = point1.dot(a1.cross(b1))/c1.dot(a1.cross(b1));
			double t = point1.dot(b1.cross(c1))/a1.dot(b1.cross(c1));
			double r = point1.dot(c1.cross(a1))/b1.dot(c1.cross(a1));
			assert abs(1 - s - t - r) < 1e-9: String.format("%.4f + %.4f + %.4f != 1", s, t, r);
			
			// use the triangular coordinates to place the point in the output triangle plane
			Vector normal2 = a2.minus(b2).cross(c2.minus(a2));
			Matrix basis = Matrix.fromColumns(
					a2.cross(b2).times(1/c2.dot(a2.cross(b2))),
					b2.cross(c2).times(1/a2.dot(b2.cross(c2))),
					normal2.times(1/a2.dot(normal2))
			).transpose();
			Vector point2 = Vector.fromMatrix(basis.inverse().times(new Vector(s, t, 1.)));
			
			// convert the point back to spherical coordinates (this implicitly applies an inverse gnomonic projection)
			return new double[] {
					atan(point2.getElement(2)/hypot(point2.getElement(0), point2.getElement(1))),
					atan2(point2.getElement(1), point2.getElement(0)),
			};
		}
		
		
		/**
		 * load the data file. this is only called if and when the user tries to use this Projection, so the data file
		 * won’t be loaded until it's needed.
		 *
		 * @param params an ignored parameter
		 * @throws IllegalArgumentException if the data files can’t be found or are in the wrong format
		 */
		@Override
		public void initialize(double... params) throws IllegalArgumentException {
			// start by checking if we’ve already done all this (it only needs to be done once)
			if (this.lookupTableTransformed != null)
				return;
			
			// load the files
			double[][] verticesInitial, verticesTransformed;
			int[][] triangles;
			try {
				verticesInitial = readNumpyFloatArray("res/" + this.filename + "_vertices_initial.npy");
				verticesTransformed = readNumpyFloatArray("res/" + this.filename + "_vertices_transformed.npy");
				triangles = readNumpyIntArray("res/" + this.filename + "_triangles.npy");
			} catch (IOException e) {
				throw new IllegalArgumentException("couldn't read the config file for this Liquid Projection: " + e);
			}
			if (verticesInitial.length != verticesTransformed.length)
				throw new IllegalArgumentException("the initial vertices file and the transformed vertices files don't match.");
			
			// convert the arrays to the proper formats
			this.verticesInitial = new Vector[verticesInitial.length];
			for (int i = 0; i < verticesInitial.length; i ++)
				this.verticesInitial[i] = new Vector(verticesInitial[i]);
			this.verticesTransformed = new Vector[verticesTransformed.length];
			for (int i = 0; i < verticesTransformed.length; i ++)
				this.verticesTransformed[i] = new Vector(verticesTransformed[i]);
			this.triangles = triangles;
			
			// bin the triangles
			this.фBins = linspace(-PI/2, PI/2, 50);
			this.λBins = linspace(-PI, PI, 100);
			this.lookupTableInitial = cacheTriangleLocations(
					this.triangles, this.verticesInitial, this.фBins, this.λBins);
			this.lookupTableTransformed = cacheTriangleLocations(
					this.triangles, this.verticesTransformed, this.фBins, this.λBins);
			
			// finally, initialize the base projection
			this.baseProjection.initialize();
			this.shape = baseProjection.getShape();
		}
		
		
		private static List<int[]>[][] cacheTriangleLocations(
				int[][] triangles, Vector[] vertices, double[] фBins, double[] λBins) {
			@SuppressWarnings("unchecked")
			List<int[]>[][] cache = (LinkedList<int[]>[][]) new LinkedList[фBins.length - 1][λBins.length - 1];
			// initialize the cache
			for (int i = 0; i < фBins.length - 1; i ++) {
				for (int j = 0; j < λBins.length - 1; j ++) {
					cache[i][j] = new LinkedList<>();
				}
			}
			// look thru the triangles
			for (int[] triangle: triangles) {
				double фMin = POSITIVE_INFINITY;
				double фMax = NEGATIVE_INFINITY;
				double λMin = POSITIVE_INFINITY;
				double λMax = NEGATIVE_INFINITY;
				double λPrevius = 0;
				for (int k = 0; k < 3; k ++) {
					double x = vertices[triangle[k]].getElement(0);
					double y = vertices[triangle[k]].getElement(1);
					double z = vertices[triangle[k]].getElement(2);
					double ф = atan(z/hypot(x, y));
					double λ = atan2(y, x);
					λ = coerceAngle(λ - λPrevius) + λPrevius; // put λ within ±π of the last λ if you can
					if (ф < фMin)
						фMin = ф;
					if (ф > фMax)
						фMax = ф;
					if (λ < λMin)
						λMin = λ;
					if (λ > λMax)
						λMax = λ;
					λPrevius = λ;
				}
				// really wide triangles might wrap around the antimeridian, so widen them manually
				if (λMax - λMin >= PI) {
					λMin = -PI;
					λMax = PI;
				}
				// wide triangles can also bulge poleward, so extend them toward the pole manually
				if (фMax > 0)
					фMax = min(фMax + (λMax - λMin)/2, PI/2); // this is a pretty conservative estimate of the amount of bulging
				if (фMin < 0)
					фMin = max(фMin - (λMax - λMin)/2, -PI/2);
				// now you can identify the cells that this triangle spans
				int iStart = (int) floor((фMin - фBins[0])/(фBins[1] - фBins[0])); // inclusive
				int iStop = (int) ceil((фMax - фBins[0])/(фBins[1] - фBins[0])); // exclusive
				int jStart = (int) floor((λMin - λBins[0])/(λBins[1] - λBins[0])); // inclusive
				int jStop = (int) ceil((λMax - λBins[0])/(λBins[1] - λBins[0])); // exclusive
				for (int i = iStart; i < iStop; i ++)
					for (int j = jStart; j < jStop; j ++)
						cache[i][floorMod(j, cache[i].length)].add(triangle);
			}
			return cache;
		}
		
		
		/**
		 * read the header of a NPY file and return the shape.
		 */
		private static double[][] readNumpyFloatArray(String filename) throws IOException {
			try (
					FileInputStream fileStream = new FileInputStream(filename);
					DataInputStream in = new DataInputStream(fileStream)
			) {
				int[] shape = readNumpyHeader(in);
				double[][] data = new double[shape[0]][shape[1]];
				for (int i = 0; i < data.length; i ++)
					for (int j = 0; j < data[i].length; j ++)
						data[i][j] = Double.longBitsToDouble(readLittleEndianInteger(in, 8));
				return data;
			}
		}
		
		
		/**
		 * read the header of a NPY file and return the shape.
		 */
		private static int[][] readNumpyIntArray(String filename) throws IOException {
			try (
					FileInputStream fileStream = new FileInputStream(filename);
					DataInputStream in = new DataInputStream(fileStream)
			) {
				int[] shape = readNumpyHeader(in);
				int[][] data = new int[shape[0]][shape[1]];
				for (int i = 0; i < data.length; i ++)
					for (int j = 0; j < data[i].length; j ++)
						data[i][j] = (int) readLittleEndianInteger(in, 4);
				return data;
			}
		}
		
		
		private static int[] readNumpyHeader(DataInputStream in) throws IOException {
			if (in.skip(6) != 6)
				throw new IllegalArgumentException("I can't find the NPY header.");
			byte[] version = new byte[2];
			if (in.read(version) != version.length)
				throw new IllegalArgumentException("I couldn't read the version number.");
			int header_len;
			if (version[0] == 1)
				header_len = (int) readLittleEndianInteger(in, 2);
			else
				header_len = (int) readLittleEndianInteger(in, 4);
			byte[] headerBytes = new byte[header_len];
			if (in.read(headerBytes) != headerBytes.length)
				throw new IllegalArgumentException("I couldn't read the whole header.");
			String header = new String(headerBytes);
			Pattern shapePattern = Pattern.compile("'shape': \\(([0-9]+), ([0-9]+)\\)");
			Matcher headerMatcher = shapePattern.matcher(header);
			if (!headerMatcher.find())
				throw new IllegalArgumentException("No match found");
			return new int[] {
				Integer.parseInt(headerMatcher.group(1)),
				Integer.parseInt(headerMatcher.group(2)),
			};
		}
		
		
		private static long readLittleEndianInteger(DataInputStream in, int numBytes) throws IOException {
			long result = 0;
			for (int i = 0; i < numBytes; i ++) {
				long nextByte = in.read();
				if (nextByte == -1)
					throw new IOException("end of file hit");
				else
					result += nextByte << 8*i;
			}
			return result;
		}
		
		
		private static double[] linspace(double start, double end, int numSteps) {
			double[] result = new double[numSteps + 1];
			for (int i = 0; i <= numSteps; i ++)
				result[i] = start + i*(end - start)/numSteps;
			return result;
		}
	}
}
