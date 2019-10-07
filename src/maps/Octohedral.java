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
import utils.Math2;

/**
 * A class of maps that use octohedral octants. Very similar to Polyhedral, but much faster since
 * it takes advantage of the fact that everything is orthogonal.
 * 
 * @author jkunimune
 */
public class Octohedral {
	
	public static final Projection WATERMAN = new OctohedralProjection(
			"Waterman Butterfly", "A simple Cahill-esque octohedral map arrangement, with Antarctica left on.",
			2*Math.sqrt(3), (Math.sqrt(3)-1)/2, 0b1010, Property.COMPROMISE, 3,
			Configuration.BUTTERFLY) {
		
		protected double[] faceProject(double lat, double lon) {
			return Waterman.faceProject(lat, lon);
		}
		
		protected double[] faceInverse(double x, double y) {
			return Waterman.faceInverse(x, y);
		}
	};
	
	
	public static final Projection KEYES_BUTTERFLY = new OctohedralProjection(
			"Cahill\u2013Keyes Butterfly", "A simple Cahill-esque octohedral map arrangement, with Antarctica left on.",
			CahillKeyes.lMG, CahillKeyes.lMA, 0b1010, Property.COMPROMISE, 4,
			Configuration.BUTTERFLY) {
		
		protected double[] faceProject(double lat, double lon) {
			return CahillKeyes.faceProjectD(Math.toDegrees(lat), Math.toDegrees(lon));
		}
		
		protected double[] faceInverse(double x, double y) {
			double[] coords = CahillKeyes.faceInverseD(x, y);
			return (coords == null) ? null :
				new double[] {Math.toRadians(coords[0]), Math.toRadians(coords[1])};
		}
	};
	
	
	public static final Projection KEYES_BASIC_M = new OctohedralProjection(
			"Cahill\u2013Keyes Basic", "A simple M-shaped octohedral projection, with Antarctica broken into three pieces.",
			CahillKeyes.lMG, CahillKeyes.lMA, 0b1010, Property.COMPROMISE, 3,
			Configuration.M_PROFILE) {
		
		protected double[] faceProject(double lat, double lon) {
			return CahillKeyes.faceProjectD(Math.toDegrees(lat), Math.toDegrees(lon));
		}
		
		protected double[] faceInverse(double x, double y) {
			double[] coords = CahillKeyes.faceInverseD(x, y);
			return (coords == null) ? null :
				new double[] {Math.toRadians(coords[0]), Math.toRadians(coords[1])};
		}
	};
	
	
	public static final Projection CAHILL_KEYES = new OctohedralProjection(
			"Cahill\u2013Keyes", "An M-shaped octohedral projection with Antarctica assembled in the center.",
			CahillKeyes.lMG, CahillKeyes.lMA, 0b1010, Property.COMPROMISE, 4,
			Configuration.M_W_S_POLE) {
		
		public double[] project(double lat, double lon) {
			return super.project(lat, lon + Math.PI/9); // apply the central meridian manually
		}
		
		protected double[] faceProject(double lat, double lon) {
			return CahillKeyes.faceProjectD(Math.toDegrees(lat), Math.toDegrees(lon));
		}
		
		public double[] inverse(double x, double y) {
			double[] coords = super.inverse(x, y);
			if (coords == null)	return null;
			coords[1] = Math2.floorMod(coords[1] - Math.PI/9 + Math.PI, 2*Math.PI) - Math.PI; // apply the central meridian manually
			return coords;
		}
		
		protected double[] faceInverse(double x, double y) {
			double[] coords = CahillKeyes.faceInverseD(x, y);
			return (coords == null) ? null :
				new double[] {Math.toRadians(coords[0]), Math.toRadians(coords[1])};
		}
	};
	
	
	public static final OctohedralProjection CONFORMAL_CAHILL = new OctohedralProjection(
			"Cahill Conformal", "The conformal and only reproducable variant of Cahill's original map.",
			Math.sqrt(3)/2, 0, 0b1000, Property.CONFORMAL, 3, Configuration.BUTTERFLY) {

		private final double HEXAGON_SCALE = 1.112913; //this is 2^(2/3)/6*\int_0^\pi sin^(-1/3) x dx
		private final double TOLERANCE = 1e-3;
		private final double[] VERTEX = {0, Math.PI/4, -3*Math.PI/4}; // TODO this needs to be tilted a bit
		
		protected double[] faceProject(double lat, double lon) {
			double[] poleCoords = {lat, lon};
			double[] vertCoords = obliquifySphc(lat, lon, VERTEX); //look at an oblique aspect from the nearest vertex
			if (poleCoords[0] > vertCoords[0]) { //if this point is closer to the pole
				Complex w = Complex.fromPolar(Math.pow(Math.tan(Math.PI/4-lat/2), 2/3.), lon*2/3.);
				Complex z = polynomial(w); //project it as normal
				return new double[] {z.getRe(), z.getIm()};
			}
			else { //if it is closer to the vertex
				Complex w = Complex.fromPolar(
						Math.pow(Math.tan(Math.PI/4-vertCoords[0]/2), 2/3.), vertCoords[1]*2/3.);
				Complex zSkew = polynomial(w); //use the maclaurin series centred there
				return new double[] {
						-1/2.*zSkew.getRe() + Math.sqrt(3)/2*zSkew.getIm() + Math.sqrt(3)/2,
						-Math.sqrt(3)/2*zSkew.getRe() - 1/2.*zSkew.getIm() + 1/2. };
			}
		}
		
		protected double[] faceInverse(double x, double y) {
			Complex z;
			if (x < (1-y)/Math.sqrt(3)) //do the Newton Raphson from whichever vertex to which it is closest
				z = new Complex(x, y);
			else
				z = new Complex(-1/2.*(x-Math.sqrt(3)/2) - Math.sqrt(3)/2*(y-1/2.),
						Math.sqrt(3)/2*(x-Math.sqrt(3)/2) - 1/2.*(y-1/2.));
			Complex w = z.divide(HEXAGON_SCALE);
			Complex error = polynomial(w).minus(z);
			for (int i = 0; i < 8 && error.abs() > TOLERANCE; i ++) {
				Complex dzdw = derivative(w);
				w = w.minus(error.divide(dzdw));
				error = polynomial(w).minus(z);
			}
			double[] latLon = { Math.PI/2 - 2*Math.atan(Math.pow(w.abs(), 3/2.)), w.arg()*3/2. }; //inverse conic it back to spherical coordinates
			if (x < (1-y)/Math.sqrt(3)) //if it was closest to that vertex, the result is easy
				return latLon;
			else //if it was closer to the other vertex, do some obliquifying
				return obliquifyPlnr(latLon, VERTEX);
		}
		
		private Complex polynomial(Complex w) { //an approximation of the true conformal mapping function
			w = w.times(Complex.fromPolar(1, -Math.PI/6));
			Complex z = w.plus(w.pow(7).divide(21))
					.plus(w.pow(11).divide(99)).plus(w.pow(13).divide(1287/16.));
			return z.divide(Complex.fromPolar(HEXAGON_SCALE, -Math.PI/6));
		}
		
		private Complex derivative(Complex w) { //the derivative of polynomial()
			w = w.times(Complex.fromPolar(1, -Math.PI/6));
			Complex z = new Complex(1).plus(w.pow(6).divide(3))
					.plus(w.pow(10).divide(9)).plus(w.pow(12).divide(99/16.));
			return z.divide(Complex.fromPolar(HEXAGON_SCALE, -Math.PI/6));
		}
	};
	
	
	public static final Projection CAHILL_CONCIALDI = new OctohedralProjection(
			"Cahill\u2013Concialdi Bat", "A conformal octohedral projection with no cuts and a unique arrangement.",
			Math.sqrt(3)/2, 0, 0b1000, Property.CONFORMAL, 4, Configuration.BAT_SHAPE) {
		
		private final double lon0 = Math.toRadians(20);
		private final double correction = Math.toRadians(9);
		private final double tilt = Math.toRadians(4.7);
		private final double trueWidth = width*Math.cos(tilt);
		
		public double[] project(double lat, double lon) {
			lon = Math2.floorMod(lon - lon0 + Math.PI, 2*Math.PI) - Math.PI; // first, change the central meridian
			
			double[] coords = super.project(lat, lon);
			
			if (lat >= 0 && lon > Math.PI-correction) { // next, fix the Bering strait
				double xi = coords[0], yi = coords[1] - size/2;
				double xf = Math.cos(2*Math.PI/3)*xi - Math.sin(2*Math.PI/3)*yi;
				double yf = Math.sin(2*Math.PI/3)*xi + Math.cos(2*Math.PI/3)*yi;
				coords = new double[] { xf, yf + size/2};
			}
			
			double xi = coords[0], yi = coords[1] - size; // and rotate 5.5 degrees
			double xf = Math.cos(tilt)*xi - Math.sin(tilt)*yi;
			double yf = Math.sin(tilt)*xi + Math.cos(tilt)*yi;
			coords = new double[] { xf, yf + size};
			
			coords[0] -= 1.5*size*Math.sin(tilt);
			if (coords[0] < -trueWidth/2) { // finally, move the cropped right stuff to the left
				coords[0] += trueWidth;
				coords[1] += trueWidth*Math.tan(tilt);
			}
			
			return coords;
		}
		
		protected double[] faceProject(double lat, double lon) {
			return CONFORMAL_CAHILL.faceProject(lat, lon);
		}
		
		public double[] inverse(double lat, double lon) {
			double[] symmetric = super.inverse(lat, lon);
			
			if (symmetric != null)
				symmetric[1] += lon0;
			return symmetric;
		}
		
		protected double[] faceInverse(double x, double y) {
			return CONFORMAL_CAHILL.faceInverse(x, y);
		}
	};
	
	
	
	private static abstract class OctohedralProjection extends Projection {
		
		protected final double size;
		private Configuration config;
		
		
		public OctohedralProjection(String name, String desc, double altitude, double cutSize,
				int fisc, Property property, int rating, Configuration config) {
			super(name, desc,
					config.fullWidth*altitude-config.cutWidth*cutSize,
					config.fullHeight*altitude-config.cutHeight*cutSize, fisc,
					(cutSize == 0) ? Type.OCTOHEDRAL : Type.TETRADECAHEDRAL, property, rating,
					new String[] {}, new double[][] {}, config.hasAspect);
			this.size = altitude;
			this.config = config;
			this.config.setCutRatio(cutSize/altitude);
		}
		
		
		protected abstract double[] faceProject(double lat, double lon);
		
		protected abstract double[] faceInverse(double x, double y);
		
		
		public double[] project(double lat, double lon) {
			for (double[] octant: config.octants) { // try each octant
				double lonr = lon - octant[3]; // relative longitude
				if (Math.abs(lonr) > Math.PI/4 || lat < octant[5] || lat > octant[6]) // if it doesn't fit...
					continue; // check the next one
				if (octant[4] != 0 && (lonr < 0) != (octant[4] < 0)) // also check the sign restriction
					continue;
				double xP = octant[0]*size, yP = octant[1]*size; // vertex coordinates
				double th = octant[2]; // rotation angle
				
				double[] coords = this.faceProject(Math.abs(lat), Math.abs(lonr));
				double xMj = coords[0], yMj = coords[1]; //relative octant coordinates (Mj stands for "Mary Jo Graca")
				
				if (lat < 0) //reflect the southern hemisphere over the equator
					xMj = 2*size - xMj;
				if (lonr < 0)
					yMj = -yMj;
				
				return new double[] {
						xP + Math.sin(th)*xMj + Math.cos(th)*yMj,
						yP - Math.cos(th)*xMj + Math.sin(th)*yMj + config.fullHeight*size/2 };
			}
			return new double[] {Double.NaN, Double.NaN}; // if none of the octants fit, return null
		}
		
		
		public double[] inverse(double x, double y) {
			y -= config.fullHeight*size/2;
			
			for (double[] octant: config.octants) { // try each octant
				double xV = octant[0]*size, yV = octant[1]*size; // vertex coordinates
				double th = octant[2]; // rotation angle
				
				double xMj = Math.sin(th)*(x-xV) - Math.cos(th)*(y-yV); // do the coordinate change
				double yMj = Math.cos(th)*(x-xV) + Math.sin(th)*(y-yV);
				if (Math.sqrt(3)*Math.abs(yMj) > Math.min(xMj, 2*size-xMj)) // if the angle is wrong,
					continue; // check the next one
				if (octant[4] != 0 && (yMj < 0) != (octant[4] < 0)) // also check the sign restriction
					continue;
				
				double[] coords = this.faceInverse(Math.min(xMj, 2*size-xMj), Math.abs(yMj));
				if (coords == null)
					continue; // if you got nothing, keep looking
				double lat = coords[0], lon = coords[1]; // project
				
				lat *= Math.signum(size-xMj); // undo the reflections
				if (lat < octant[5] || lat > octant[6]) // if the resulting coordinates are wrong
					continue; // move on
				
				lon = Math.signum(yMj)*lon + octant[3];
				
				return new double[] { lat, lon };
			}
			return null;
		}
		
	}
	
	
	
	private enum Configuration {
		
		BUTTERFLY(4, 2, 4/Math.sqrt(3), Math.sqrt(3), true, new double[][] { //the classic four quadrants splayed out in a nice butterfly shape, with Antarctica divided and attached
			{  0, -1/Math.sqrt(3), -Math.PI/2, -3*Math.PI/4, 0, -Math.PI/2,  Math.PI/2 },
			{  0, -1/Math.sqrt(3), -Math.PI/6,   -Math.PI/4, 0, -Math.PI/2,  Math.PI/2 },
			{  0, -1/Math.sqrt(3),  Math.PI/6,    Math.PI/4, 0, -Math.PI/2,  Math.PI/2 },
			{  0, -1/Math.sqrt(3),  Math.PI/2,  3*Math.PI/4, 0, -Math.PI/2,  Math.PI/2 },
		}),
		
		M_PROFILE(4, 0, Math.sqrt(3), Math.sqrt(3), true, new double[][] { //The more compact zigzag configuration with Antarctica divided and attached
			{ -1,               0, -Math.PI/6, -3*Math.PI/4, 0, -Math.PI/2,  Math.PI/2 },
			{ -1,               0,  Math.PI/6,   -Math.PI/4, 0, -Math.PI/2,  Math.PI/2 },
			{  1,               0, -Math.PI/6,    Math.PI/4, 0, -Math.PI/2,  Math.PI/2 },
			{  1,               0,  Math.PI/6,  3*Math.PI/4, 0, -Math.PI/2,  Math.PI/2 },
		}),
		
		M_W_S_POLE(4, 0, 2.008, Math.sqrt(3), false, new double[][] { //Keyes's current configuration, with Antarctica reassembled in the center
			{ -1,               0, -Math.PI/6, -3*Math.PI/4, 0, Math.toRadians(-64),  Math.PI/2 },
			{ -1,               0,  Math.PI/6,   -Math.PI/4, 0,          -Math.PI/2,  Math.PI/2 },
			{  1,               0, -Math.PI/6,    Math.PI/4, 0, Math.toRadians(-64),  Math.PI/2 },
			{  1,               0,  Math.PI/6,  3*Math.PI/4, 0, Math.toRadians(-64),  Math.PI/2 },
			{ -1.6976,  -2.6036,  2*Math.PI/3, -3*Math.PI/4, 0, -Math.PI/2, Math.toRadians(-64) },
			{  1.6036,  -0.6976,   -Math.PI/3,    Math.PI/4, 0, -Math.PI/2, Math.toRadians(-64) },
			{  0.9060,  -3.3013, -5*Math.PI/6,  3*Math.PI/4, 0, -Math.PI/2, Math.toRadians(-64) },
		}),
		
		BAT_SHAPE(2*Math.sqrt(3), 0, 2, 0, false, new double[][] { //Luca Concialdi's obscure "Bat" arrangement that I liked.
			{               0, -0.5, -2*Math.PI/3, -Math.PI  ,  0,          0,  Math.PI/2 },
			{ -3/Math.sqrt(3),  0.5,            0, -Math.PI  ,  0, -Math.PI/2,          0 },
			{               0, -0.5,   -Math.PI/3, -Math.PI/2,  0, -Math.PI/2,  Math.PI/2 },
			{               0, -0.5,            0,          0,  0, -Math.PI/4,  Math.PI/2 },
			{               0, -2.5, -2*Math.PI/3,          0, -1, -Math.PI/2, -Math.PI/4 },
			{               0, -2.5,  2*Math.PI/3,          0,  1, -Math.PI/2, -Math.PI/4 },
			{               0, -0.5,    Math.PI/3,  Math.PI/2,  0, -Math.PI/2,  Math.PI/2 },
			{               0, -0.5,  2*Math.PI/3,  Math.PI  ,  0,          0,  Math.PI/2 },
			{  3/Math.sqrt(3),  0.5,            0,  Math.PI  ,  0, -Math.PI/2,          0 },
		});
		
		public final double fullWidth, cutWidth, fullHeight, cutHeight;
		public final boolean hasAspect;
		public final double[][] octants; // array of {x, y (from top), rotation, \lambda_0, sign, \phi_min, \phi_max}
//		public double cutRatio; //this variable should be set by the map projection so that the configuration knows how big the cuts actually are
		
		private Configuration(double fullWidth, double cutWidth,
				double fullHeight, double cutHeight, boolean hasAspect,
				double[][] octants) {
			this.fullWidth = fullWidth; //the size of the configuration in altitudes ignoring cuts
			this.cutWidth = cutWidth; //the change in width due to the cut in cut lengths
			this.fullHeight = fullHeight; //and the cut sizes
			this.cutHeight = cutHeight;
			this.hasAspect = hasAspect;
			this.octants = octants;
		}
		
		public void setCutRatio(double cutRatio) {
//			this.cutRatio = cutRatio; //the length of the cut in altitudes
		}
	}
	
}
