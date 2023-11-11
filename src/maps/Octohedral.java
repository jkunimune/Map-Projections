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

import de.jtem.mfc.field.Complex;
import maps.Projection.Property;
import utils.BoundingBox;

import static java.lang.Double.NaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;
import static utils.Math2.floorMod;
import static utils.Math2.min;

/**
 * A class of maps that use octohedral octants. Very similar to Polyhedral, but much faster since
 * it takes advantage of the fact that everything is orthogonal.
 * 
 * @author jkunimune
 */
public class Octohedral {
	
	public static final Projection WATERMAN = new OctohedralProjection(
			"Waterman Butterfly", "A simple Cahill-esque octohedral map arrangement, with Antarctica left on.",
			2*sqrt(3), (sqrt(3)-1)/2, 0b1010, Property.COMPROMISE, 3,
			Configuration.BUTTERFLY) {
		
		protected double[] faceProject(double lat, double lon) {
			return Waterman.faceProject(lat, lon);
		}
		
		protected double[] faceInverse(double x, double y) {
			return Waterman.faceInverse(x, y);
		}
	};


	public static final Projection KEYES_BASIC_M = new OctohedralProjection(
			"Cahill\u2013Keyes (simplified)", "A simple M-shaped octohedral projection, with Antarctica broken into three pieces.",
			CahillKeyes.lMG, CahillKeyes.lMA, 0b1010, Property.COMPROMISE, 3,
			Configuration.M_PROFILE) {
		
		protected double[] faceProject(double lat, double lon) {
			return CahillKeyes.faceProjectD(toDegrees(lat), toDegrees(lon));
		}
		
		protected double[] faceInverse(double x, double y) {
			double[] coords = CahillKeyes.faceInverseD(x, y);
			return (coords == null) ? null :
				new double[] {toRadians(coords[0]), toRadians(coords[1])};
		}
	};


	public static final Projection KEYES_STANDARD = new OctohedralProjection(
			"Cahill\u2013Keyes", "An M-shaped octohedral projection with Antarctica assembled in the center.",
			CahillKeyes.lMG, CahillKeyes.lMA, 0b1010, Property.COMPROMISE, 4,
			Configuration.M_W_S_POLE) {

		public double[] project(double lat, double lon) {
			return super.project(lat, lon + PI/9); // apply the central meridian manually
		}

		protected double[] faceProject(double lat, double lon) {
			return CahillKeyes.faceProjectD(toDegrees(lat), toDegrees(lon));
		}

		public double[] inverse(double x, double y) {
			double[] coords = super.inverse(x, y);
			if (coords == null)	return null;
			coords[1] = floorMod(coords[1] - PI/9 + PI, 2*PI) - PI; // apply the central meridian manually
			return coords;
		}

		protected double[] faceInverse(double x, double y) {
			double[] coords = CahillKeyes.faceInverseD(x, y);
			return (coords == null) ? null :
					new double[] {toRadians(coords[0]), toRadians(coords[1])};
		}
	};


	public static final Projection KEYES_OCTANT = new OctohedralProjection(
			"Cahill\u2013Keyes (single octant)", "A single octant of the Cahill\u2013Keyes projection (for memory economization in the case of very large maps).",
			CahillKeyes.lMG, CahillKeyes.lMA, 0b1010, Property.COMPROMISE, 3,
			Configuration.SINGLE_OCTANT) {

		protected double[] faceProject(double lat, double lon) {
			return CahillKeyes.faceProjectD(toDegrees(lat), toDegrees(lon));
		}

		protected double[] faceInverse(double x, double y) {
			double[] coords = CahillKeyes.faceInverseD(x, y);
			return (coords == null) ? null :
					new double[] {toRadians(coords[0]), toRadians(coords[1])};
		}
	};


	public static final OctohedralProjection CONFORMAL_CAHILL = new OctohedralProjection(
			"Cahill Conformal", "The conformal and only reproducable variant of Cahill's original map.",
			sqrt(3)/2, 0, 0b1000, Property.CONFORMAL, 3, Configuration.BUTTERFLY) {

		private final double HEXAGON_SCALE = 1.112913; //this is 2^(2/3)/6*\int_0^\pi sin^(-1/3) x dx
		private final double TOLERANCE = 1e-3;
		private final double[] VERTEX = {0, PI/4, -3*PI/4}; // TODO this needs to be tilted a bit
		
		protected double[] faceProject(double lat, double lon) {
			double[] poleCoords = {lat, lon};
			double[] vertCoords = transformFromOblique(lat, lon, VERTEX); //look at an oblique aspect from the nearest vertex
			if (poleCoords[0] > vertCoords[0]) { //if this point is closer to the pole
				Complex w = Complex.fromPolar(pow(tan(PI/4-lat/2), 2/3.), lon*2/3.);
				Complex z = polynomial(w); //project it as normal
				return new double[] {z.getRe(), z.getIm()};
			}
			else { //if it is closer to the vertex
				Complex w = Complex.fromPolar(
						pow(tan(PI/4-vertCoords[0]/2), 2/3.), vertCoords[1]*2/3.);
				Complex zSkew = polynomial(w); //use the maclaurin series centred there
				return new double[] {
						-1/2.*zSkew.getRe() + sqrt(3)/2*zSkew.getIm() + sqrt(3)/2,
						-sqrt(3)/2*zSkew.getRe() - 1/2.*zSkew.getIm() + 1/2. };
			}
		}
		
		protected double[] faceInverse(double x, double y) {
			Complex z;
			if (x < (1-y)/sqrt(3)) //do the Newton Raphson from whichever vertex to which it is closest
				z = new Complex(x, y);
			else
				z = new Complex(-1/2.*(x-sqrt(3)/2) - sqrt(3)/2*(y-1/2.),
						sqrt(3)/2*(x-sqrt(3)/2) - 1/2.*(y-1/2.));
			Complex w = z.divide(HEXAGON_SCALE);
			Complex error = polynomial(w).minus(z);
			for (int i = 0; i < 8 && error.abs() > TOLERANCE; i ++) {
				Complex dzdw = derivative(w);
				w = w.minus(error.divide(dzdw));
				error = polynomial(w).minus(z);
			}
			double[] latLon = { PI/2 - 2*atan(pow(w.abs(), 3/2.)), w.arg()*3/2. }; //inverse conic it back to spherical coordinates
			if (x < (1-y)/sqrt(3)) //if it was closest to that vertex, the result is easy
				return latLon;
			else //if it was closer to the other vertex, do some obliquifying
				return transformToOblique(latLon, VERTEX);
		}
		
		private Complex polynomial(Complex w) { //an approximation of the true conformal mapping function
			w = w.times(Complex.fromPolar(1, -PI/6));
			Complex z = w.plus(w.pow(7).divide(21))
					.plus(w.pow(11).divide(99)).plus(w.pow(13).divide(1287/16.));
			return z.divide(Complex.fromPolar(HEXAGON_SCALE, -PI/6));
		}
		
		private Complex derivative(Complex w) { //the derivative of polynomial()
			w = w.times(Complex.fromPolar(1, -PI/6));
			Complex z = new Complex(1).plus(w.pow(6).divide(3))
					.plus(w.pow(10).divide(9)).plus(w.pow(12).divide(99/16.));
			return z.divide(Complex.fromPolar(HEXAGON_SCALE, -PI/6));
		}
	};
	
	
	public static final Projection CAHILL_CONCIALDI = new OctohedralProjection(
			"Cahill\u2013Concialdi", "A conformal octohedral projection with no extra cuts and a unique arrangement.",
			sqrt(3)/2, 0, 0b1000, Property.CONFORMAL, 4, Configuration.BAT_SHAPE) {
		
		private final double lon0 = toRadians(20);
		
		public double[] project(double lat, double lon) {
			lon = floorMod(lon - lon0 + PI, 2*PI) - PI; // change the central meridian
			return super.project(lat, lon);
		}
		
		protected double[] faceProject(double lat, double lon) {
			return CONFORMAL_CAHILL.faceProject(lat, lon);
		}
		
		public double[] inverse(double x, double y) {
			double[] coords = super.inverse(x, y);
			if (coords != null)
				coords[1] = floorMod(coords[1] + lon0 + PI, 2*PI) - PI; // change the central meridian
			return coords;
		}
		
		protected double[] faceInverse(double x, double y) {
			return CONFORMAL_CAHILL.faceInverse(x, y);
		}
	};
	
	
	
	private static abstract class OctohedralProjection extends Projection {
		
		protected final double size;
		private final Configuration config;
		
		
		public OctohedralProjection(String name, String desc, double altitude, double cutSize,
				int fisc, Property property, int rating, Configuration config) {
			super(name, desc,
			      new BoundingBox(config.fullWidth*altitude - config.cutWidth*cutSize,
			                      config.fullHeight*altitude - config.cutHeight*cutSize),
			      fisc, (cutSize == 0) ? Type.OCTOHEDRAL : Type.TETRADECAHEDRAL, property, rating,
                  new String[] {}, new double[][] {}, config.hasAspect);
			this.size = altitude;
			this.config = config;
		}
		
		
		protected abstract double[] faceProject(double lat, double lon);
		
		protected abstract double[] faceInverse(double x, double y);
		
		
		public double[] project(double lat, double lon) {
			for (Octant octant: config.octants) { // try each octant
				double lonr = floorMod(lon - octant.centralLongitude + PI, 2*PI) - PI; // relative longitude
				if (abs(lonr) > PI/4 || lat < octant.minLatitude || lat > octant.maxLatitude) // if it doesn't fit...
					continue; // check the next one
				if (lonr < octant.minLongitude || lonr > octant.maxLongitude) // also check the longitude restriction
					continue;
				double xP = octant.x*size, yP = octant.y*size; // vertex coordinates
				double th = octant.planeRotation; // rotation angle
				
				double[] coords = this.faceProject(abs(lat), abs(lonr));
				double xMj = coords[0], yMj = coords[1]; //relative octant coordinates (Mj stands for "Mary Jo Graca")
				
				if (lat < 0) //reflect the southern hemisphere over the equator
					xMj = 2*size - xMj;
				if (lonr < 0)
					yMj = -yMj;
				
				return new double[] {
						xP + sin(th)*xMj + cos(th)*yMj,
						yP - cos(th)*xMj + sin(th)*yMj + config.fullHeight*size/2 };
			}
			return new double[] {NaN, NaN}; // if none of the octants fit, return null
		}
		
		
		public double[] inverse(double x, double y) {
			y -= config.fullHeight*size/2;
			
			for (Octant octant: config.octants) { // try each octant
				double xV = octant.x*size, yV = octant.y*size; // vertex coordinates
				double th = octant.planeRotation; // rotation angle
				
				double xMj = sin(th)*(x-xV) - cos(th)*(y-yV); // do the coordinate change
				double yMj = cos(th)*(x-xV) + sin(th)*(y-yV);
				if (sqrt(3)*abs(yMj) > min(xMj, 2*size-xMj) + 1e-12) // if the angle is wrong,
					continue; // check the next one
				
				double[] coords = this.faceInverse(min(xMj, 2*size-xMj), abs(yMj));
				if (coords == null)
					continue; // if you got nothing, keep looking
				double lat = coords[0], lon = coords[1]; // project
				
				lat *= signum(size-xMj); // undo the reflections
				lon *= signum(yMj);
				if (lat < octant.minLatitude - 1e-6 || lat > octant.maxLatitude + 1e-6) // if the resulting coordinates are wrong
					continue; // move on
				if (lon < octant.minLongitude - 1e-6 || lon > octant.maxLongitude + 1e-6) // also check the longitude restriction
					continue;
				
				lon = floorMod(lon + octant.centralLongitude + PI, 2*PI)- PI;
				
				return new double[] { lat, lon };
			}
			return null;
		}
		
	}
	
	
	
	private enum Configuration {

		/** the classic four quadrants splayed out in a nice butterfly shape, with Antarctica divided and attached */
		BUTTERFLY(4, 2, 4/sqrt(3), sqrt(3), true, new Octant[] {
			new Octant(0, -1/sqrt(3), -PI/2, -PI/2, PI/2, -3*PI/4),
			new Octant(0, -1/sqrt(3), -PI/6, -PI/2, PI/2,   -PI/4),
			new Octant(0, -1/sqrt(3),  PI/6, -PI/2, PI/2,    PI/4),
			new Octant(0, -1/sqrt(3),  PI/2, -PI/2, PI/2,  3*PI/4),
		}),
		/** The more compact zigzag configuration with Antarctica divided and attached */
		M_PROFILE(4, 0, sqrt(3), sqrt(3), true, new Octant[] {
			new Octant(-1, 0, -PI/6, -PI/2, PI/2, -3*PI/4),
			new Octant(-1, 0,  PI/6, -PI/2, PI/2,   -PI/4),
			new Octant( 1, 0, -PI/6, -PI/2, PI/2,    PI/4),
			new Octant( 1, 0,  PI/6, -PI/2, PI/2,  3*PI/4),
		}),
		/** Gene Keyes's current configuration, with Antarctica reassembled in the center */
		M_W_S_POLE(4, 0, 2.008, sqrt(3), false, new Octant[] {
			new Octant(-1.    ,   0.    ,   -PI/6, toRadians(-64),           PI/2, -3*PI/4),
			new Octant(-1.    ,   0.    ,    PI/6,          -PI/2,           PI/2,   -PI/4),
			new Octant( 1.    ,   0.    ,   -PI/6, toRadians(-64),           PI/2,    PI/4),
			new Octant( 1.    ,   0.    ,    PI/6, toRadians(-64),           PI/2,  3*PI/4),
			new Octant(-1.6976,  -2.6036,  2*PI/3,          -PI/2, toRadians(-64), -3*PI/4),
			new Octant( 1.6036,  -0.6976,   -PI/3,          -PI/2, toRadians(-64),    PI/4),
			new Octant( 0.9060,  -3.3013, -5*PI/6,          -PI/2, toRadians(-64),  3*PI/4),
		}),
		/** Luca Concialdi's "Bat" arrangement */
		BAT_SHAPE(3.47, 0, 2.03, 0, false, rotate(0, -1.55, toRadians(5), new Octant[] {
			new Octant(         0, -0.5, -2*PI/3,   0  ,  PI/2, -PI  , toRadians(-9),  1),
			new Octant(-3/sqrt(3),  0.5,       0, -PI/2,   0  , -PI  , toRadians( 9),  1),
			new Octant(         0, -0.5,   -PI/3, -PI/2,  PI/2, -PI/2),
			new Octant(         0, -0.5,       0, -PI/4,  PI/2,   0  ),
			new Octant(         0, -2.5, -2*PI/3, -PI/2, -PI/4,   0  , -1,  0),
			new Octant(         0, -2.5,  2*PI/3, -PI/2, -PI/4,   0  ,  0,  1),
			new Octant(         0, -0.5,    PI/3, -PI/2,  PI/2,  PI/2),
			new Octant(         0, -0.5,  2*PI/3,   0  ,  PI/2,  PI  , -1, toRadians(-9)),
			new Octant( 3/sqrt(3),  0.5,       0, -PI/2,   0  ,  PI  , -1, toRadians( 9)),
		})),
		/** an octohedron that actually only covers the positive octant, in case you want to do each octant separately */
        SINGLE_OCTANT(1, 0, 2/sqrt(3), sqrt(3), true, new Octant[] {
			new Octant(-0.5, 0, PI/6, 0, PI/2, PI/4),
		});
		
		public final double fullWidth, cutWidth, fullHeight, cutHeight;
		public final boolean hasAspect;
		/**
		 * array of {
		 *  x, y, rotation, λ_0, ф_min, ф_max[, λ_min, λ_max]
		 * } for each face of the octohedron
		 */
		public final Octant[] octants;

		/**
		 * @param fullWidth the size of the configuration in altitudes ignoring cuts
		 * @param cutWidth the change in width due to the cut-in lengths
		 * @param fullHeight the heit of the configuration in altitudes ignoring cuts
		 * @param cutHeight the change in heit due to the cut-in lengths
		 * @param hasAspect whether it would make any sense to change the aspect
		 * @param octants the array of octant specifications.  each octant must
		 *                specify the x and y of the pole, the central meridian,
		 *                and bounding latitudes and longitudes of the face.
		 */
		Configuration(double fullWidth, double cutWidth,
				double fullHeight, double cutHeight, boolean hasAspect,
				Octant[] octants) {
			this.fullWidth = fullWidth;
			this.cutWidth = cutWidth;
			this.fullHeight = fullHeight;
			this.cutHeight = cutHeight;
			this.hasAspect = hasAspect;
			this.octants = octants;
		}
		
		private static Octant[] rotate(double xP, double yP, double th, Octant[] in) { // apply that rotation to the bat
			Octant[] out = new Octant[in.length];
			for (int i = 0; i < in.length; i ++) {
				out[i] = new Octant(
						cos(th)*(in[i].x - xP) - sin(th)*(in[i].y - yP) + xP,
						sin(th)*(in[i].x - xP) + cos(th)*(in[i].y - yP) + yP,
						in[i].planeRotation + th,
						in[i].minLatitude, in[i].maxLatitude,
						in[i].centralLongitude,
						in[i].minLongitude, in[i].maxLongitude);
			}
			return out;
		}
	}



	private static class Octant {
		public final double x;
		public final double y;
		public final double planeRotation;
		public final double minLatitude;
		public final double maxLatitude;
		public final double centralLongitude;
		public final double minLongitude;
		public final double maxLongitude;

		Octant(double x, double y, double planeRotation, double minLatitude, double maxLatitude, double centralLongitude) {
			this(x, y, planeRotation, minLatitude, maxLatitude, centralLongitude, -PI, PI);
		}

		Octant(double x, double y, double planeRotation, double minLatitude, double maxLatitude,
		       double centralLongitude, double minLongitude, double maxLongitude) {
			this.x = x;
			this.y = y;
			this.planeRotation = planeRotation;
			this.minLatitude = minLatitude;
			this.maxLatitude = maxLatitude;
			this.centralLongitude = centralLongitude;
			this.minLongitude = minLongitude;
			this.maxLongitude = maxLongitude;
		}
	}
	
}
