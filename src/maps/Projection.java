/**
 * MIT License
 * 
 * Copyright (c) 2017 Justin Kunimune
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package maps;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.function.BooleanSupplier;
import java.util.function.DoubleConsumer;

import image.SVGMap.Command;
import image.SVGMap.Path;
import utils.BoundingBox;

import static java.lang.Double.NaN;
import static java.lang.Double.isNaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.cos;
import static java.lang.Math.hypot;
import static java.lang.Math.log;
import static java.lang.Math.round;
import static java.lang.Math.sin;
import static java.lang.Math.tan;
import static utils.Math2.coerceAngle;
import static utils.Math2.hypot;
import static utils.Math2.linInterp;
import static utils.Math2.max;
import static utils.Math2.mean;
import static utils.Math2.outOfBoundsInSameDirection;
import static utils.Math2.rms;
import static utils.Math2.stdDev;

/**
 * An object that transforms coordinates between spheres and planes.
 * 
 * @author jkunimune
 */
public abstract class Projection {
	
	public static final double[] NORTH_POLE = {PI/2, 0, 0};
	
	
	private final String name; //typically the name of the dude credited for it
	private final String description; //a noun clause or sentence about it
	
	private final String[] paramNames; //the name of each parameter
	private final double[][] paramValues; //the bounds and default value of each parameter
	private final boolean hasAspect; //is it spherically symmetrical?
	
	private final boolean finite; //does it display the entire world?
	private final boolean invertable; //is the inverse solution closed-form?
	private final boolean solveable; //is the solution closed-form?
	private final boolean continuous; //does a random continuous path cross outside of the map?
	private final Type type; //the geometry of the projection
	private final Property property; //what it is good for
	private final int rating; //how good I think it is
	protected BoundingBox bounds; //bounding box of the full map
	
	
	
	protected Projection(
			String name, BoundingBox bounds, int fisc, Type type, Property property, int rating) {
		this(name, buildDescription(type,property,null,null),
		     bounds, fisc, type, property, rating, new String[0], new double[0][]);
	}
	
	protected Projection(
			String name, BoundingBox bounds, int fisc, Type type, Property property,
			int rating, String adjective) {
		this(name, buildDescription(type,property,adjective,null), bounds,
		     fisc, type, property, rating, new String[0], new double[0][]);
	}
	
	protected Projection(
			String name, BoundingBox bounds, int fisc, Type type, Property property,
			int rating, String adjective, String addendum) {
		this(name, buildDescription(type,property,adjective,addendum), bounds,
		     fisc, type, property, rating, new String[0], new double[0][]);
	}
	
	protected Projection(
			String name, String description, BoundingBox bounds, int fisc,
			Type type, Property property, int rating) {
		this(name, description, bounds, fisc, type, property, rating,
		     new String[0], new double[0][]);
	}
	
	protected Projection(
			String name, String description, BoundingBox bounds, int fisc, Type type,
			Property property, int rating, String[] paramNames, double[][] paramValues) {
		this(name, description, bounds, fisc, type, property, rating,
		     paramNames, paramValues, true);
	}
	
	protected Projection(
			String name, String description, BoundingBox bounds, int fisc, Type type,
			Property property, int rating, String[] paramNames, double[][] paramValues,
			boolean hasAspect) {
		this(name, description, bounds,
		     (fisc&0b1000) > 0, (fisc&0b0100) > 0, (fisc&0b0010) > 0, (fisc&0b0001) > 0,
		     type, property, rating, paramNames, paramValues, hasAspect);
	}
	
	protected Projection (
			String name, String description, BoundingBox bounds,
			boolean finite, boolean invertable, boolean solveable, boolean continuous, Type type, Property property, int rating,
			String[] paramNames, double[][] paramValues, boolean hasAspect) {
		this.name = name;
		this.description = description;
		this.paramNames = paramNames;
		this.paramValues = paramValues;
		this.hasAspect = hasAspect;
		this.bounds = bounds;
		this.finite = finite;
		this.invertable = invertable;
		this.solveable = solveable;
		this.continuous = continuous;
		this.type = type;
		this.property = property;
		this.rating = rating;
	}
	
	protected Projection(String name, Projection base) {
		this(name, base.description, base.bounds, base.finite, base.invertable,
		     base.solveable, base.continuous, base.type, base.property, base.rating,
		     base.paramNames, base.paramValues, base.hasAspect);
	}
	
	private static String buildDescription(Type type, Property property, String adjective, String addendum) { //these should all be lowercase
		String description = property+" "+type+" projection";
		if (adjective != null)
			description = adjective+" "+description;
		if (addendum != null)
			description += " "+addendum;
		if (description.charAt(0) == 'a' || description.charAt(0) == 'e' || description.charAt(0) == 'i' || description.charAt(0) == 'o' || description.charAt(0) == 'u')
			return "An "+description+".";
		else
			return "A "+description+".";
	}


	/**
	 * convert a location on the globe to a location on the map plane
	 * @param lat the latitude in radians
	 * @param lon the longitude in radians
	 * @return the x value and y value in the same units as this.bounds
	 */
	public abstract double[] project(double lat, double lon);

	/**
	 * convert a location on the map plane to a location on the globe
	 * @param x the x value in the same units as this.bounds
	 * @param y the y value in the same units as this.bounds
	 * @return the latitude and longitude in radians
	 */
	public abstract double[] inverse(double x, double y);


	/**
	 * load in any user-specified parameters and perform any tasks that must be done before calling project() and
	 * inverse(). this method only gets called if and when the user selects this Projection, so it’s a good place to put
	 * computationally expensive setup stuff.
	 * @param params a list of numbers that may relate to the map projection. each parameterized Projection will have
	 *               its own interpretation for these. most Projections will take an empty array for this.
	 * @throws IllegalArgumentException if the number of parameters is not what the Projection was expecting
	 */
	public void initialize(double... params) throws IllegalArgumentException {
	}


	/**
	 * convert a location on the globe to a location on the map plane
	 * @param coords the latitude and longitude in radians
	 * @return the x value and y value in the same units as this.bounds
	 */
	public double[] project(double[] coords) {
		return project(coords[0], coords[1]);
	}

	/**
	 * convert a location on the rotated globe to a location on the map plane
	 * @param coords the absolute latitude and longitude in radians
	 * @param pole the desired aspect: the latitude and longitude of the location that should appear as the North Pole
	 *             and the angle to rotate the globe about its new axis
	 * @return the x value and y value in the same units as this.bounds
	 */
	public double[] project(double[] coords, double[] pole) {
		return project(coords[0], coords[1], pole);
	}

	/**
	 * convert a location on the rotated globe to a location on the map plane
	 * @param lat the absolute latitude in radians
	 * @param lon the absolute longitude in radians
	 * @param pole the desired aspect: the latitude and longitude of the location that should appear as the North Pole
	 *             and the angle to rotate the globe about its new axis
	 * @return the x value and y value in the same units as this.bounds
	 */
	public double[] project(double lat, double lon, double[] pole) {
		return project(transformFromOblique(lat, lon, hasAspect ? pole : null));
	}

	/**
	 * convert a location on the map plane to a location on the globe
	 * @param coords the x and y value, in the same units as this.bounds
	 * @return the latitude and longitude in radians
	 */
	public double[] inverse(double[] coords) {
		return inverse(coords[0], coords[1]);
	}

	/**
	 * convert a location on the map plane to a location on the rotated globe
	 * @param coords the x and y value, in the same units as this.bounds
	 * @param pole the desired aspect: the latitude and longitude of the location that should appear as the North Pole
	 *             and the angle to rotate the globe about its new axis
	 * @return the latitude and longitude in radians
	 */
	public double[] inverse(double[] coords, double[] pole) {
		return inverse(coords[0], coords[1], pole);
	}

	/**
	 * convert a location on the map plane to a location on the rotated globe
	 * @param x the x value in the same units as this.bounds
	 * @param y the y value in the same units as this.bounds
	 * @param pole the desired aspect: the latitude and longitude of the location that should appear as the North Pole
	 *             and the angle to rotate the globe about its new axis
	 * @return the latitude and longitude in radians
	 */
	public double[] inverse(double x, double y, double[] pole) {
		return inverse(x, y, pole, false);
	}

	/**
	 * convert a location on the map plane to a location on the rotated globe
	 * @param x the x value in the same units as this.bounds
	 * @param y the y value in the same units as this.bounds
	 * @param pole the desired aspect: the latitude and longitude of the location that should appear as the North Pole
	 *             and the angle to rotate the globe about its new axis
	 * @param cropAtPi whether to forbid longitudes outside of [-π, π], returning NaN for any points that can only be
	 *                 reached by extrapolating the projection to such extreme longitudes. if false, these points will
	 *                 be considered valid and those extreme longitudes will be returned like normal. for funnily-shaped
	 *                 interrupted projections, this option will typically restrict the result so that only one x and y
	 *                 returns each latitude and longitude, whereas setting it false may cause some globe locations to
	 *                 be duplicated.
	 * @return the latitude and longitude in radians
	 */
	public double[] inverse(double x, double y, double[] pole, boolean cropAtPi) {
		final double[] relCoords = inverse(x, y);
		if (relCoords == null || (cropAtPi && abs(relCoords[1]) > PI))
			return null; //cropAtPi removes all points with longitudes outside +- PI
		else
			return transformToOblique(relCoords, hasAspect ? pole : null);
	}


	/**
	 * perform the inverse projection on a 2D array of points
	 * @param size the maximum linear dimension of the 2D array
	 * @return an array of latitude-longitude pairs, where each row is at a particular x value,
	 *         and each column is at a particular y value.  y decreases with increasing row index.
	 *         elements corresponding to points not on the map will be set to zero.
	 */
	public double[][][] map(int size) {
		return map(size, false);
	}

	/**
	 * perform the inverse projection on a 2D array of points
	 * @param size the maximum linear dimension of the 2D array
	 * @param cropAtPi whether to insist that each point on the globe is included exactly once
	 * @return an array of latitude-longitude pairs, where each row is at a particular x value,
	 *         and each column is at a particular y value.  y decreases with increasing row index.
	 *         elements corresponding to points not on the map will be set to zero.
	 */
	public double[][][] map(int size, boolean cropAtPi) {
		return map(size, null, cropAtPi);
	}

	/**
	 * perform the inverse projection on a 2D array of points
	 * @param size the maximum linear dimension of the 2D array
	 * @param pole the aspect to use for the projection
	 * @param cropAtPi whether to insist that each point on the globe is included exactly once
	 * @return an array of latitude-longitude pairs, where each row is at a particular x value,
	 *         and each column is at a particular y value.  y decreases with increasing row index.
	 *         elements corresponding to points not on the map will be set to zero.
	 */
	public double[][][] map(int size, double[] pole, boolean cropAtPi) {
		double width = bounds.xMin - bounds.xMax;
		double height = bounds.yMin - bounds.yMax;
		if (width >= height)
			return map(size, max(round(size*height/width),1), pole, cropAtPi, null);
		else
			return map(max(round(size*width/height),1), size, pole, cropAtPi, null);
	}

	/**
	 * perform the inverse projection on a 2D array of points
	 * @param w the number of columns in the 2D array
	 * @param h the number of rows in the 2D array
	 * @param pole the aspect to use for the projection
	 * @param cropAtPi whether to insist that each point on the globe is included exactly once
	 * @return an array of latitude-longitude pairs, where each row is at a particular x value,
	 *         and each column is at a particular y value.  y decreases with increasing row index.
	 *         elements corresponding to points not on the map will be set to zero.
	 */
	public double[][][] map(double w, double h, double[] pole, boolean cropAtPi,
			DoubleConsumer tracker) { //generate a matrix of coordinates based on a map projection
		final double[][][] output = new double[(int) h][(int) w][2];
		for (int y = 0; y < h; y ++) {
			for (int x = 0; x < w; x ++)
				output[y][x] = inverse(
						bounds.xMin + (x + 0.5)/w*(bounds.xMax - bounds.xMin),
						bounds.yMax - (y + 0.5)/h*(bounds.yMax - bounds.yMin), pole, cropAtPi);
			if (tracker != null)
				tracker.accept((double)y / (int)h);
		}
		return output;
	}
	
	
	/**
	 * Create a series of paths that draw a graticule mesh
	 * @param spacing The number of radians between each parallel or meridian
	 * @param precision The maximum allowable distance from the true path
	 * @param maxLat The maximum absolute value of latitude for any graticule curve
	 * @param maxLon The maximum absolute value of longitude for any graticule curve
	 * @param outW The maximum width of this graticule
	 * @param outH The maximum height of this graticule
	 * @param pole The aspect of this graticule
	 * @return list of curves where each curve is a list of {x,y} arrays
	 */
	public Path drawGraticule(double spacing, double precision, double outW, double outH,
			double maxLat, double maxLon, double[] pole) {
		Path output = new Path();
		
		for (int y = 0; y < (int)(maxLat/spacing); y ++) {
			output.addAll(drawLoxodrome( //northern parallel
					 y*spacing,-maxLon, y*spacing, maxLon, precision, outW, outH, pole));
			if (y == 0) 	continue;
			output.addAll(drawLoxodrome( //southern parallel
					-y*spacing,-maxLon,-y*spacing, maxLon, precision, outW, outH, pole));
		}
		maxLat -= .0001; //don't draw on the poles; it makes things easier
		for (int x = 0; x <= (int)(maxLon/spacing); x ++) {
			output.addAll(drawLoxodrome( //western meridian
					-maxLat,-x*spacing, maxLat,-x*spacing, precision, outW, outH, pole));
			if (x == 0 || x == (int)(maxLon/spacing)) 	continue;
			output.addAll(drawLoxodrome( //eastern meridian
					-maxLat, x*spacing, maxLat, x*spacing, precision, outW, outH, pole));
		}
		
		return output;
	}
	
	
	private Path drawLoxodrome(double lat0, double lon0, double lat1, double lon1,
			double precision, double outW, double outH, double[] pole) {
		final double[][] baseRange = {{bounds.xMin, bounds.yMax}, {bounds.xMax, bounds.yMin}};
		final double[][] imgRange = {{0, 0}, {outW, outH}}; //define some constants for changing coordinates
		
		double[] endPt0 = new double[] {lat0, lon0};
		double[] endPt1 = new double[] {lat1, lon1};
		List<double[]> spherical = new ArrayList<double[]>(); //the spherical coordinates of the vertices
		for (double a = 0; a <= 1; a += 1/32.) //populated with vertices along the loxodrome
			spherical.add(new double[] {endPt0[0]*a+endPt1[0]*(1-a), endPt0[1]*a+endPt1[1]*(1-a)});
		Path planar = new Path(); //the planar coordinates of the vertices
		for (int i = 0; i < spherical.size(); i ++) {
			double[] si = spherical.get(i); //populated with projections of spherical, in image coordinates
			double[] pi = linInterp(this.project(si, pole), baseRange, imgRange);
			char type = (i == 0) ? 'M' : 'L';
			planar.add(new Command(type, pi));
		}
		
		Queue<double[]> queue = new LinkedList<double[]>(spherical.subList(0, spherical.size()-1));
		double[] s0;
		while ((s0 = queue.poll()) != null) { //now iteratively flesh out the rest
			int i = spherical.indexOf(s0); //s0 is the first spherical endpoint
			double[] s1 = spherical.get(i+1); //second spherical endpoint
			double[] sm = new double[] {(s0[0]+s1[0])/2, (s0[1]+s1[1])/2}; //spherical (loxodromic) midpoint
			double[] p0 = planar.get(i).args; //first planar endpoint
			double[] p1 = planar.get(i+1).args; //second planar endpoint
			if (outOfBoundsInSameDirection(imgRange, p0, p1)) // if we're talking about things entirely off the map
				continue; // just forget about it
			double[] pm = linInterp(this.project(sm, pole), baseRange, imgRange); //planar (loxodromic) midpoint
			
			double error = hypot(pm[0] - (p0[0] + p1[0])/2, pm[1] - (p0[1] + p1[1])/2); // midpoint error
			if (error > precision) { //if the calculated midpoint is too far from what we expect
				if ((i-1 < 0 || hypot(planar.get(i-1).args, p0) <= precision) &&
						(i+2 >= planar.size() || hypot(planar.get(i+2).args, p1) <= precision)) { // check if it's getting real close on each side
					planar.set(i+1, new Command('M', p1)); // if so, it's probably an interruption. Change the second one to 'M'.
					continue;
				}
				else if (hypot(s1[0] - s0[0], s1[1] - s0[1]) < 1e-4) { // alternatively, if we're getting to arcsecond scale,
					planar.set(i+1, new Command('M', p1)); // it's just not worth it
					continue;
				}
				spherical.add(i+1, sm); //if there's still work to do, add the midpoint to the curve
				planar.add(i+1, new Command('L', pm));
				queue.add(s0); //and see if you need to recurse this at all
				queue.add(sm);
			}
		}
		
		return planar;
	}
	
	
	public static double[][][] globe(double dt) { //generate a matrix of coordinates based on the sphere
		List<double[]> points = new ArrayList<double[]>();
		for (double phi = -PI/2+dt/2; phi < PI/2; phi += dt) { // make sure phi is never exactly +-tau/4
			for (double lam = -PI+dt/cos(phi)/2; lam < PI; lam += dt/cos(phi)) {
				points.add(new double[] {phi, lam});
			}
		}
		return new double[][][] {points.toArray(new double[0][])};
	}
	
	
	public static double[][][] hemisphere(double dt) { //like globe(), but for the eastern hemisphere. Good for doing projections that are symmetrical in longitude (i.e. pretty much all of them)
		List<double[]> points = new ArrayList<double[]>();
		for (double phi = -PI/2+dt/2; phi < PI/2; phi += dt) { // make sure phi is never exactly +-tau/4
			for (double lam = dt/cos(phi)/2; lam < PI; lam += dt/cos(phi)) {
				points.add(new double[] {phi, lam});
			}
		}
		return new double[][][] {points.toArray(new double[0][])};
	}
	
	
	public double[] avgDistortion(double[][][] points, double[] params) {
		this.initialize(params);
		return avgDistortion(points);
	}
	
	public double[] avgDistortion(double[][][] points) {
		final double[][][] distDist = calculateDistortion(points);
		return new double[] {stdDev(distDist[0]), rms(distDist[1])};
	}
	
	
	public double[][][] calculateDistortion(double[][][] points) {
		return calculateDistortion(points, () -> false, (d) -> {});
	}
	
	public double[][][] calculateDistortion(double[][][] points,
			BooleanSupplier cancelation, DoubleConsumer progressTracker) { //calculate both kinds of distortion over the given region
		double[][][] output = new double[2][points.length][points[0].length]; //the distortion matrix
		
		for (int y = 0; y < points.length; y ++) {
			if (cancelation.getAsBoolean()) 	return null;
			progressTracker.accept((double)y/points.length);
			for (int x = 0; x < points[y].length; x ++) {
				if (points[y][x] != null) {
					final double[] dists = getDistortionAt(points[y][x]);
					output[0][y][x] = dists[0]; //the output matrix has two layers:
					output[1][y][x] = dists[1]; //area and angular distortion
				}
				else {
					output[0][y][x] = NaN;
					output[1][y][x] = NaN; //NaN means no map here
				}
			}
		}
		
		final double avgArea = mean(output[0]); //don't forget to normalize output[0] so the average is zero
		for (int y = 0; y < output[0].length; y ++)
			for (int x = 0; x < output[0][y].length; x ++)
				output[0][y][x] -= avgArea;
		
		return output;
	}
	
	
	public double[] getDistortionAt(double[] s0) { //calculate both kinds of distortion at the given point
		final double[] output = new double[2];
		final double dx = 1e-8;
		
		final double[] sC = { s0[0]+dx, s0[1] }; //first, step to the side a bit to help us avoid interruptions
		final double[] sE = { sC[0], sC[1]+dx/cos(sC[0]) }; //consider a point slightly to the east
		final double[] sN = { sC[0]+dx, sC[1] }; //and slightly to the north
		final double[] pC = project(sC);
		final double[] pE = project(sE);
		final double[] pN = project(sN);
		
		final double dA = 
				(pE[0]-pC[0])*(pN[1]-pC[1]) - (pE[1]-pC[1])*(pN[0]-pC[0]);
		output[0] = log(abs(dA/(dx*dx))); //the zeroth output is the size (area) distortion
		if (abs(output[0]) > 25)
			output[0] = NaN; //discard outliers
		
		final double s1ps2 = hypot((pE[0]-pC[0])+(pN[1]-pC[1]), (pE[1]-pC[1])-(pN[0]-pC[0]));
		final double s1ms2 = hypot((pE[0]-pC[0])-(pN[1]-pC[1]), (pE[1]-pC[1])+(pN[0]-pC[0]));
		output[1] = abs(log(abs((s1ps2-s1ms2)/(s1ps2+s1ms2)))); //the first output is the shape (angle) distortion
		if (output[1] > 25)
			output[1] = NaN; //discard outliers
		
		return output;
	}
	
	
	/**
	 * Calculate relative latitude and longitude for an oblique pole
	 * @param latF the absolute latitude in radians
	 * @param lonF the absolute longitude in radians
	 * @param pole the pole location
	 * @return { latr, lonr }, or coords if pole is null
	 */
	protected static double[] transformFromOblique(double latF, double lonF, double[] pole) {
		if (pole == null || Arrays.equals(pole, NORTH_POLE)) // null pole indicates that this procedure should be bypassed
			return new double[] {latF, lonF};
		
		final double lat0 = pole[0];
		final double lon0 = pole[1];
		final double tht0 = pole[2];
		
		double lat1;
		if (lat0 == PI/2)
			lat1 = latF;
		else
			lat1 = asin(sin(lat0)*sin(latF) + cos(lat0)*cos(latF)*cos(lon0-lonF)); // relative latitude
		
		double lon1;
		if (lat0 == PI/2) // accounts for all the 0/0 errors at the poles
			lon1 = lonF - lon0;
		else if (lat0 == -PI/2)
			lon1 = lon0 - lonF - PI;
		else {
			lon1 = acos((cos(lat0)*sin(latF) - sin(lat0)*cos(latF)*cos(lon0-lonF))/cos(lat1))-PI; // relative longitude
			if (isNaN(lon1)) {
				if ((cos(lon0-lonF) >= 0 && latF < lat0) || (cos(lon0-lonF) < 0 && latF < -lat0))
					lon1 = 0;
				else
					lon1 = -PI;
			}
			else if (sin(lonF - lon0) > 0) // it's a plus-or-minus arccos.
				lon1 = -lon1;
		}
		lon1 = lon1-tht0;
		if (abs(lon1) > PI) //put all longitudes in [-pi,pi], for convenience
			lon1 = coerceAngle(lon1);
		if (lon1 >= PI - 1e-7) // finally, kill any roundoff error on the edge
			lon1 = -PI;
		
		return new double[] {lat1, lon1};
	}
	
	
	/**
	 * Calculate absolute latitude and longitude for an oblique pole
	 * @param coords the relative coordinates
	 * @param pole the pole location
	 * @return { LAT, LON }, or coords if pole is null
	 */
	protected static double[] transformToOblique(double[] coords, double[] pole) {
		if (pole == null) //this indicates that you just shouldn't do this calculation
			return coords;
		
		double lat1 = coords[0], lon1 = coords[1];
		final double lat0 = pole[0], lon0 = pole[1], tht0 = pole[2];
		
		lon1 += tht0;
		double latf = asin(sin(lat0)*sin(lat1) - cos(lat0)*cos(lon1)*cos(lat1));
		double lonf;
		double innerFunc = sin(lat1)/cos(lat0)/cos(latf) - tan(lat0)*tan(latf);
		if (lat0 == PI/2) // accounts for special case when lat0 = pi/2
			lonf = lon1+lon0;
		else if (lat0 == -PI/2) // accounts for special case when lat0 = -pi/2
			lonf = -lon1+lon0 + PI;
		else if (abs(innerFunc) > 1) { // accounts for special case when cos(lat1) -> 0
			if ((lon1 == 0 && lat1 < -lat0) || (lon1 != 0 && lat1 < lat0))
				lonf = lon0 + PI;
			else
				lonf = lon0;
		}
		else if (sin(lon1) > 0)
			lonf = lon0 + acos(innerFunc);
		else
			lonf = lon0 - acos(innerFunc);
		
		if (abs(lonf) > PI)
			lonf = coerceAngle(lonf);
		
		double thtf = pole[2];
		
		return new double[] {latf, lonf, thtf};
	}
	
	
	@Override
	public String toString() {
		return this.getName();
	}
	
	
	public final Projection transverse() {
		return transverse(getName());
	}
	
	public final Projection transverse(String name) {
		return new Oblique(this, name, 0, 0, 0);
	}
	
	public final Projection withAspect(String name ,double... aspect) {
		return new Oblique(this, name, aspect);
	}
	
	
	public final String getName() {
		return this.name;
	}
	
	public final String getDescription() {
		return this.description;
	}
	
	public final boolean isParametrized() {
		return this.paramNames.length > 0;
	}
	
	public final int getNumParameters() {
		return this.paramNames.length;
	}
	
	public final String[] getParameterNames() {
		return this.paramNames;
	}
	
	public final double[] getDefaultParameters() {
		final double[] params = new double[this.getNumParameters()];
		for (int i = 0; i < this.getNumParameters(); i ++)
			params[i] = this.paramValues[i][2];
		return params;
	}
	
	public final double[][] getParameterValues() {
		return this.paramValues;
	}
	
	public final boolean hasAspect() {
		return this.hasAspect;
	}
	
	public final boolean isFinite() {
		return this.finite;
	}
	
	public final boolean isInvertable() {
		return this.invertable;
	}
	
	public final boolean isSolveable() {
		return this.solveable;
	}
	
	public final boolean isContinuous() {
		return this.continuous;
	}
	
	public final Type getType() {
		return this.type;
	}
	
	public final Property getProperty() {
		return this.property;
	}
	
	/**
	 * @return the rating:
	 * 0 if I hate it with a burning passion;
	 * 1 if it is bad and you shouldn't use it;
	 * 2 if it has its use cases but isn't very good outside of them;
	 * 3 if it is a solid, defensible choice; and
	 * 4 if I love it with a burning passion.
	 */
	public final int getRating() {
		return this.rating;
	}

	public final BoundingBox getBounds() {
		return this.bounds;
	}
	
	public final double getAspectRatio() {
		return (bounds.xMax - bounds.xMin)/(bounds.yMax - bounds.yMin);
	}

	public static final Projection NULL_PROJECTION = //this exists solely for the purpose of a "More..." option at the end of menus
			new Projection("More...", null, null, 0, null, null, 0) {
		
		public double[] project(double lat, double lon) {
			return null;
		}
		
		public double[] inverse(double x, double y) {
			return null;
		}
		
	};
	
	
	
	/**
	 * The most common geometric configurations of projections
	 * @author jkunimune
	 */
	public enum Type {
		CYLINDRICAL("Cylindrical"), CONIC("Conic"), AZIMUTHAL("Azimuthal"),
		PSEUDOCYLINDRICAL("Pseudocylindrical"), PSEUDOCONIC("Pseudoconic"),
		PSEUDOAZIMUTHAL("Pseudoazimuthal"),
		TETRAHEDRAL("Tetrahedral"), OCTOHEDRAL("Octohedral"),
		TETRADECAHEDRAL("Truncated Octohedral"), ICOSOHEDRAL("Icosohedral"),
		POLYNOMIAL("Polynomial"), STREBE("Strebe blend"), PLANAR("Planar"), OTHER("Other");
		
		private final String name;
		
		Type(String name) {
			this.name = name;
		}
		
		public String toString() {
			return this.name.toLowerCase();
		}
		
		public String getName() {
			return this.name;
		}
	}
	
	
	/**
	 * The useful quantities that projections can preserve
	 * @author jkunimune
	 */
	public enum Property {
		CONFORMAL("Conformal"), EQUIDISTANT("Equidistant"), EQUAL_AREA("Equal-area"),
		PERSPECTIVE("Perspective"), GNOMONIC("Gnomonic"), RETROAZIMUTHAL("Retroazimuthal"),
		COMPROMISE("Compromise"), POINTLESS("Pointless"), TRUE("True");
		
		private final String name;
		
		Property(String name) {
			this.name = name;
		}
		
		public String toString() {
			return this.name.toLowerCase();
		}
		
		public String getName() {
			return this.name;
		}
	}
}
