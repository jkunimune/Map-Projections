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

import maps.Projection.Property;
import maps.Projection.Type;
import utils.BoundingBox;
import utils.Dixon;
import utils.NumericalAnalysis;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.asin;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.hypot;
import static java.lang.Math.pow;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.toRadians;
import static utils.Math2.coerceAngle;
import static utils.Math2.floorMod;

/**
 * Projections created by projecting onto and then unfolding some kind of polyhedron
 * 
 * @author jkunimune
 */
public class Polyhedral {
	
	public static final PolyhedralProjection LEE_TETRAHEDRAL_RECTANGULAR = new PolyhedralProjection(
			"Lee Tetrahedral", 0b1001, Configuration.TETRAHEDRON_WIDE_FACE, Property.CONFORMAL,
			4, null, "that really deserves more attention") {
		
		public double[] faceProject(double lat, double lon) {
			final de.jtem.mfc.field.Complex z = de.jtem.mfc.field.Complex.fromPolar(
					pow(2, 5/6.)*tan(PI/4-lat/2), lon);
			final de.jtem.mfc.field.Complex w = Dixon.invFunc(z);
			return new double[] { w.abs()*2/Dixon.PERIOD_THIRD, w.arg() }; //I don't understand Dixon functions well enough to say whence the 1.132 comes
		}
		
		public double[] faceInverse(double r, double tht) {
			final de.jtem.mfc.field.Complex w = de.jtem.mfc.field.Complex.fromPolar(
					r*Dixon.PERIOD_THIRD/2, tht);
			final de.jtem.mfc.field.Complex ans = Dixon.leeFunc(w).times(pow(2, -5/6.));
			return new double[] {
					PI/2 - 2*atan(ans.abs()),
					ans.arg() };
		}
	};
	
	
	public static final PolyhedralProjection LEE_TETRAHEDRAL_TRIANGULAR = new PolyhedralProjection(
			"Lee Tetrahedral (triangular)", 0b1001, Configuration.TRIANGLE_FACE, Property.CONFORMAL,
			2, null,
			"in a triangle, because this is the form in which it was published, even though the rectangle is clearly better") {
		
		public double[] faceProject(double lat, double lon) {
			return LEE_TETRAHEDRAL_RECTANGULAR.faceProject(lat, lon);
		}
		
		public double[] faceInverse(double r, double tht) {
			return LEE_TETRAHEDRAL_RECTANGULAR.faceInverse(r, tht);
		}
	};
	
	
	public static final PolyhedralProjection TETRAGRAPH = new PolyhedralProjection(
			"TetraGraph", 0b1111, Configuration.TETRAHEDRON_WIDE_FACE, Property.EQUIDISTANT,
			2, null, "that I invented") {
		
		public double[] faceProject(double lat, double lon) {
			return new double[] {
					atan(1/tan(lat)*cos(lon))/cos(lon) /atan(sqrt(2)),
					lon };
		}
		
		public double[] faceInverse(double r, double tht) {
			return new double[] {
					PI/2 - atan(tan(r*cos(tht)*atan(sqrt(2)))/cos(tht)),
					tht };
		}
	};
	
	
	public static final PolyhedralProjection AUTHAGRAPH = new PolyhedralProjection(
			"IMAGO (AuthaGraph)", "Authagraph is a hip new Japanese map that would be super great if "
					+ "they actually published their equations. This is technically just an approximation, also known as the Infinitessimal "
					+ "Mutated AuthaGraph Offspring.",
			0b1011, Configuration.AUTHAGRAPH, Property.COMPROMISE, 3,
			new String[] {"Power"}, new double[][] {{.5,1,.68}}) {
		
		private final double[] POLE = {toRadians(77), toRadians(143), toRadians(17)};
		private double k;
		
		public void initialize(double... params) {
				this.k = params[0];
		}
		
		@Override
		public double[] project(double lat, double lon) { //apply a pole shift to AuthaGraph
			double[] relCoords = transformFromOblique(lat, lon, POLE);
			return super.project(relCoords[0], relCoords[1]);
		}
		
		@Override
		public double[] inverse(double x, double y) { //because AuthaGraph needs its obliquity, and I didn't want to program that into the Configuration
			return transformToOblique(super.inverse(x, y), POLE);
		}
		
		
		public double[] faceProject(double lat, double lon) {
			final double tht = atan((lon - asin(sin(lon)/sqrt(3)))/PI*sqrt(12));
			final double p = (PI/2 - lat) / atan(sqrt(2)/cos(lon));
			return new double[] { pow(p,k)*sqrt(3)/cos(tht), tht };
		}
		
		protected double[] faceInverse(double r, double th) {
			final double lon = NumericalAnalysis.newtonRaphsonApproximation(th, th*2,
					(l) -> atan((l - asin(sin(l)/sqrt(3)))/PI*sqrt(12)),
					(l) -> (1-1/sqrt(1+2*pow(cos(l),-2)))/sqrt(pow(PI,2)/12+pow(l-asin(sin(l)/sqrt(3)),2)),
					.001);
			final double R = r / (sqrt(3)/cos(th));
			return new double[] {
					PI/2 - pow(R,1/k)*atan(sqrt(2)/cos(lon)), lon };
		}
	};
	
	
	public static final PolyhedralProjection AUTHAPOWER = new PolyhedralProjection(
			"TetraPower", "A parametrised, simplified version of my AuthaGraph approximation.",
			0b1011, Configuration.TETRAHEDRON_WIDE_VERTEX, Property.COMPROMISE, 4,
			new String[] {"Power"}, new double[][] {{.5,1,.6}}) {
		
		private double k;
		
		public void initialize(double... params) {
				this.k = params[0];
		}
		
		public double[] faceProject(double lat, double lon) {
			final double tht = atan((lon - asin(sin(lon)/sqrt(3)))/PI*sqrt(12));
			final double p = (PI/2 - lat) / atan(sqrt(2)/cos(lon));
			return new double[] { pow(p,k)*sqrt(3)/cos(tht), tht };
		}
		
		protected double[] faceInverse(double r, double th) {
			final double lon = NumericalAnalysis.newtonRaphsonApproximation(th, th*2,
					(l) -> atan((l - asin(sin(l)/sqrt(3)))/PI*sqrt(12)),
					(l) -> (1-1/sqrt(1+2*pow(cos(l),-2)))/sqrt(pow(PI,2)/12+pow(l-asin(sin(l)/sqrt(3)),2)),
					.001);
			final double R = r / (sqrt(3)/cos(th));
			return new double[] {
					PI/2 - pow(R,1/k)*atan(sqrt(2)/cos(lon)), lon };
		}
	};
	
	
	public static final PolyhedralProjection ACTUAUTHAGRAPH = new PolyhedralProjection(
			"EquaHedral", "An interrupted authalic tetrahedral projection.",
			0b1010, Configuration.TETRAHEDRON_WIDE_VERTEX, Property.EQUAL_AREA, 3,
			new String[] {"Sinus length"}, new double[][] {{0, 60, 20}}) {
		
		private double sig, a0, scale;
		
		public void initialize(double... params) {
			this.sig = params[0]/60;
			this.a0 = 3 - 1.5*sig*sig;
			this.scale = sqrt(3)*a0/PI;
		}
		
		public double[] faceProject(double lat, double lon) {
			double bet = atan((lon-asin(sin(lon)/sqrt(3)))/(a0/2)*scale);
			double f = (1-sin(lat))/(1-1/sqrt(1+2/pow(cos(lon), 2)));
			if (f < (3*sig*sig)/(2*a0)) { //sinus zone
				double alf = atan(2*tan(abs(bet)) - 1/sqrt(3));
				double rA = sqrt(2*f*a0/pow(cos(alf), 2))/2;
				return toPolar(rA, alf, signum(lon));
			}
			else { //primary zone
				double rB = sqrt(a0*f - a0 + 3)/cos(bet);
				return new double[] {rB, bet};
			}
		}
		
		public double[] faceInverse(double r, double th) {
			double bet;
			double f;
			if (r < sig*sqrt(3)/2/cos(abs(th)-PI/3)) //empty
				return null;
			else if (r*cos(th) < sig*sqrt(3)) { //sinus zone
				double[] relCoords = fromPolar(r, th);
				double rA = relCoords[0];
				double alf = relCoords[1];
				bet = signum(th)*atan(tan(alf)/2+1/sqrt(12));
				f = pow(rA, 2) * 2*pow(cos(alf), 2)/a0;
			}
			else { //primary zone
				bet = th;
				f = (r*r*pow(cos(bet), 2) - 1.5*sig*sig)/a0;
			}
			double lon = NumericalAnalysis.newtonRaphsonApproximation(
					a0/2*tan(bet)/scale, bet*2,
					(l) -> l - asin(sin(l)/sqrt(3)),
					(l) -> 1 - 1/sqrt(1 + 2/pow(cos(l), 2)), 1e-4);
			double lat = asin(1 - f*(1 - 1/sqrt(1+2/pow(cos(lon), 2))));
			return new double[] {lat, lon};
		}
		
		private double[] toPolar(double rA, double alf, double s) {
			double x = rA*cos(alf) + sig*sqrt(3)/2;
			double y = rA*sin(alf) + sig/2;
			return new double[] {hypot(x, y), s*atan2(y, x)};
		}
		
		private double[] fromPolar(double rB, double bet) {
			double x = rB*cos(bet) - sig*sqrt(3)/2;
			double y = abs(rB*sin(bet)) - sig/2;
			return new double[] {hypot(x, y), atan2(y, x)};
		}
	};
	
	
	public static final Projection VAN_LEEUWEN = new PolyhedralProjection(
			"Van Leeuwen", "An uninterrupted equal-area tetrahedral projection. It's more accurately known as \"the Vertex-oriented great circle projection applied to a tetrahedron\", but the guy who copublished it with Leeuwen calls it \"the van Leeuwen projection\" on his website, so I think this is fine.",
			0b1011, Configuration.TETRAHEDRON_WIDE_VERTEX, Property.EQUAL_AREA, 3) {
		
		public double[] faceProject(double lat, double lon) {
			ACTUAUTHAGRAPH.initialize(0);
			return ACTUAUTHAGRAPH.faceProject(lat, lon);
		}
		
		public double[] faceInverse(double r, double th) {
			ACTUAUTHAGRAPH.initialize(0);
			return ACTUAUTHAGRAPH.faceInverse(r, th);
		}
	};
	
	
	public static final Projection DYMAXION = new PolyhedralProjection(
			"Dymaxion", "A polyhedral projection that slices up the oceans as much as possible without slicing up any landmasses.",
			0b1110, Configuration.DYMAXION, Property.COMPROMISE, 3) {
		
		private final double[] POLE = {0.040158, -0.091549,-2.015269}; //I derived these numbers from [Robert Gray](http://www.rwgrayprojects.com/rbfnotes/maps/graymap4.html)
		private final double X_0 = 0.75;
		private final double Y_0 = -sqrt(3)/4;
		
		private final double sin36 = sqrt(10-2*sqrt(5))/4;
		private final double cos36 = (1+sqrt(5))/4;
		
		@Override
		public double[] project(double lat, double lon) { //apply a pole shift and Cartesian shift to Dymaxion
			double[] coords = transformFromOblique(lat, lon, POLE);
			coords = super.project(coords[0], coords[1]);
			return new double[] {coords[0] + X_0, coords[1] + Y_0};
		}
		
		@Override
		public double[] inverse(double x, double y) { //because Dymaxion needs its obliquity, and I didn't want to program that into the Configuration
			double[] coords = super.inverse(x - X_0, y - Y_0);
			if (coords == null) 	return null;
			return transformToOblique(coords, POLE);
		}
		
		public double[] faceProject(double lat, double lon) {
			double xG = cos(lon)/tan(lat)/cos36; //normalised gnomonic coordinates
			double yG = sin(lon)/tan(lat)/sin36;
			double a = asin((xG+yG)/(2*sqrt(1+xG*xG))) + atan(xG); //angular distance up each side of the triangle
			double b = asin((xG-yG)/(2*sqrt(1+xG*xG))) + atan(xG);
			double x = (a + b)/(2*sqrt(3)); //final cartesian coordinates in radians
			double y = (a - b)/2;
			return new double[] {hypot(x,y)/atan(2), atan2(y,x)}; //scale to fit to layout, where side length is 1
		}
		
		public double[] faceInverse(double r, double th) {
			if (abs(th) > PI/6+1e-15) 	throw new IllegalArgumentException("Wait, what?"+th);
			double x = r*cos(th)*atan(2); //cartesian coordinates in radians
			double y = r*sin(th)*atan(2);
			double a = sqrt(3)*x + y; //angular distance up each side of the triangle
			double b = sqrt(3)*x - y;
			double xG = cos36*(sin(a) + sin(b))/(1 + cos(a) + cos(b)); //unnormalised gnomonic coordinates
			double yG = sin36*
					(sin(a) - sin(b) + 2*sin(a-b))/(1 + cos(a) + cos(b));
			return new double[] {atan(1/hypot(xG, yG)), atan2(yG, xG)}; //inverse gnomonic projection
		}
	};
	
	
	
	/**
	 * A base for polyhedral Projections
	 * 
	 * @author jkunimune
	 */
	private static abstract class PolyhedralProjection extends Projection {
		
		private final Configuration configuration;
		
		
		public PolyhedralProjection(
				String name, int fisc, Configuration config, Property property, int rating,
				String adjective, String addendum) {
			super(name, config.bounds, fisc, config.type, property, rating,
					adjective, addendum);
			this.configuration = config;
		}
		
		public PolyhedralProjection(
				String name, String description, int fisc, Configuration config, Property property,
				int rating) {
			super(name, description, config.bounds, fisc, config.type, property,
					rating);
			this.configuration = config;
		}
		
		public PolyhedralProjection(
				String name, String description, int fisc, Configuration config, Property property,
				int rating, String[] paramNames, double[][] paramValues) {
			super(name, description, config.bounds, fisc, config.type, property,
					rating, paramNames, paramValues);
			this.configuration = config;
		}
		
		
		protected abstract double[] faceProject(double lat, double lon); //the projection from spherical to polar within a face
		
		protected abstract double[] faceInverse(double x, double y); //I think you can guess
		
		
		public double[] project(double lat, double lon) {
			final int numSym = configuration.sphereSym; //we're about to be using this variable a lot
			double latR = NEGATIVE_INFINITY;
			double lonR = NEGATIVE_INFINITY;
			double[] centrum = null;
			for (double[] testCentrum: configuration.centrumSet) { //iterate through the centrums to see which goes here
				final double[] relCoords = transformFromOblique(lat, lon, testCentrum);
				if (testCentrum.length > 6) { //if the centrum is long, then it contains longitude bounds
					double minL = testCentrum[6]*PI/numSym;
					double maxL = testCentrum[7]*PI/numSym;
					relCoords[1] = floorMod(relCoords[1]-minL, 2*PI) + minL;
					if (relCoords[1] < minL || relCoords[1] > maxL)
						continue; //ignore any longitudes not in the bounds described in [6:7]
				}
				
				if (relCoords[0] > latR) { //pick the centrum that maxes out latitude
					latR = relCoords[0];
					lonR = relCoords[1];
					centrum = testCentrum;
				}
			}
			
			final double lonR0 = floor((lonR+PI/numSym)/(2*PI/numSym))
					*(2*PI/numSym); //because most face projections are periodic
			
			final double[] rth = faceProject(latR, lonR - lonR0); //apply the projection to the relative coordinates
			final double r = rth[0];
			final double th = rth[1] + centrum[3] + lonR0*numSym/configuration.planarSym;
			final double x0 = centrum[4];
			final double y0 = centrum[5];
			
			double[] output = { r*cos(th) + x0, r*sin(th) + y0 };
			if (output[0] < bounds.xMin || output[0] > bounds.xMax || output[1] < bounds.yMin || output[1] > bounds.yMax) { //rotate OOB bits around nearest singularity
				output = configuration.rotateOOB(output[0], output[1], x0, y0);
			}
			return output;
		}
		
		
		public double[] inverse(double x, double y) {
			if (!configuration.inBounds(x, y)) 	return null;
			
			final int numSym = configuration.planarSym; //we'll be using this variable a lot soon
			
			double rM = POSITIVE_INFINITY;
			double[] centrum = null; //iterate to see which centrum we get
			for (double[] testCentrum: configuration.centrumSet) {
				final double rR = hypot(x-testCentrum[4], y-testCentrum[5]);
				if (rR < rM) { //pick the centrum that minimises r
					rM = rR;
					centrum = testCentrum;
				}
			}
			
			final double th0 = centrum[3];
			final double x0 = centrum[4];
			final double y0 = centrum[5];
			final double r = hypot(x - x0, y - y0);
			final double th = coerceAngle(atan2(y - y0, x - x0) - th0);
			
			if (centrum.length > 6) { //if the centrum has extra values, they are angle bounds
				if (th < centrum[6]*PI/numSym || th > centrum[7]*PI/numSym)
					return null; //ignore any angles not in the bounds described in [6:7]
			}
			
			final double thBase = floor((th+PI/numSym)/(2*PI/numSym))
					*(2*PI/numSym); //because most face projections are periodic
			
			double[] relCoords = faceInverse(r, th - thBase);
			
			if (relCoords == null)
				return null;
			
			relCoords[1] = thBase*numSym/configuration.sphereSym + relCoords[1];
			double[] absCoords = transformToOblique(relCoords, centrum);
			if (abs(absCoords[1]) > PI)
				absCoords[1] = coerceAngle(absCoords[1]);
			return absCoords;
		}
	}
	
	
	
	/**
	 * A set of objects that determine the layouts of tetrahedral projections
	 * 
	 * @author jkunimune
	 */
	private enum Configuration {
		
		/*		  LATITUDE,		 LONGITUDE,		 CTR_MERID,		 PLANE_ROT,		 X,	 Y */
		TETRAHEDRON_WIDE_FACE(3, 3, new BoundingBox(6., 2*sqrt(3)), new double[][] { // [<|>] arrangement, face-centred
				{ asin(1/3.),	 PI,	-2*PI/3,	-2*PI/3,	 2,	 sqrt(3),-1,2 },
				{ asin(1/3.),	 PI,	 2*PI/3,	-PI/3,		-2,	 sqrt(3) },
				{-PI/2,	 0,			 2*PI/3,	 2*PI/3,	 2,	-sqrt(3),-2,1 },
				{-PI/2,	 0,			-2*PI/3,	 PI/3,		-2,	-sqrt(3) },
				{ asin(1/3.),	 PI/3,	-2*PI/3,	 PI,		 1,	 0 },
				{ asin(1/3.),	-PI/3,	 2*PI/3,	 0,				-1,	 0 }}),
		
		TRIANGLE_FACE(3, 3, new BoundingBox(4*sqrt(3), 6.), new double[][] { // \delta arrangement, like they are often published
				{ asin(1/3.),  PI/3,	 0,				-5*PI/6,	 sqrt(3),	2 },
				{ asin(1/3.), -PI/3,	 0,				-PI/6,		-sqrt(3),	2 },
				{ asin(1/3.),  PI,	 0,				 PI/2,		 0,			-1 },
				{-PI/2,  0,			 0,				-PI/2,		 0,			 1 }}) {
			@Override public boolean inBounds(double x, double y) {
				return y > sqrt(3)*abs(x) - 3;
			}
		},
		TETRAHEDRON_WIDE_VERTEX(3, 6, new BoundingBox(6., 2*sqrt(3)), new double[][] { // [<|>] arrangement, vertex-centred
				{ PI/2,	 0,				 0,				-PI/2,		 0,	 sqrt(3) },
				{-asin(1/3.),	 0,				 PI,		 PI/2,		 0,	-sqrt(3) },
				{-asin(1/3.),	 2*PI/3,	 PI,		 5*PI/6,	 3,	 0 },
				{-asin(1/3.),	-2*PI/3,	 PI,		 PI/6,		-3,	 0 }}) {
			@Override public double[] rotateOOB(double x, double y, double xCen, double yCen) {
				if (x < bounds.xMin || x > bounds.xMax)
					return new double[] {2*xCen - x, -y};
				else
					return new double[] {-x, 2*bounds.yMax*signum(y) - y};
			}
		},
		AUTHAGRAPH(3, 6, new BoundingBox(4*sqrt(3), 3), new double[][] { // |\/\/`| arrangement, vertex-centred
				{-asin(1/3.),	 PI,		 PI,	 0,	-2*sqrt(3)-.6096,	 1.5 },
				{-asin(1/3.),	-PI/3,		 PI/3,	 0,	-sqrt(3)-.6096,	-1.5 },
				{ PI/2,	 0,				 PI,	 0,	 0-.6096,				 1.5 },
				{-asin(1/3.),	 PI/3,		-PI/3,	 0,	 sqrt(3)-.6096,	-1.5 },
				{-asin(1/3.),	 PI,		 PI,	 0,	 2*sqrt(3)-.6096,	 1.5 },
				{-asin(1/3.),	-PI/3,		 PI/3,	 0,	 3*sqrt(3)-.6096,	-1.5 }}) {
			@Override public double[] rotateOOB(double x, double y, double xCen, double yCen) {
				if (y < bounds.yMin || y > bounds.yMax) {
					x = 2*xCen - x;
					y = 2*yCen - y;
				}
				if (x < bounds.xMin || x > bounds.xMax)
					x = floorMod(x + bounds.xMax,2*bounds.xMax) - bounds.xMax;
				return new double[] {x, y};
			}
		},
		DYMAXION(5, 6, new BoundingBox(5.5, 1.5*sqrt(3)), new double[][] { // I can't draw this in ASCII. You know what "Dymaxion" means
				{ PI/2,    0.0,		-3*PI/5,-PI/2,  -1.5, sqrt(3),  -3,3 }, //West Africa
				{ PI/2,    0.0,		 PI/5,	 -PI/2,   0.5, sqrt(3),  -1,1 }, //Brazil
				{ PI/2,    0.0,		 3*PI/5,-PI/2,   1.5, sqrt(3),  -1,1 }, //South Atlantic O.
				{ atan(1/2.),-4*PI/5, 2*PI/5,-PI/6,  -2.0, sqrt(3)/2,-5,5 }, //Arabia
				{ atan(1/2.),-2*PI/5,-2*PI/5,-5*PI/6,-1.0, sqrt(3)/2,-5,5 }, //Scandanavia
				{ atan(1/2.), 0.0,		 0.0,		 -PI/2,   0.0, sqrt(3)/2,-3,5 }, //Caribbean
				{ atan(1/2.), 0.0,		-4*PI/5,-5*PI/6,-0.5, sqrt(3),  -1,1 }, //North Atlantic O.
				{ atan(1/2.), 2*PI/5, 0.0,		 -PI/2,   1.0, sqrt(3)/2,-5,5 }, //Patagonia
				{ atan(1/2.), 4*PI/5,-2*PI/5,-5*PI/6, 2.0, sqrt(3)/2,-3,2 }, //East Antarctica
				{ atan(1/2.), 4*PI/5, 0.0,		 -PI/6,  -3.5, 0.0,		    0,1 }, //South Indian O.
				{ atan(1/2.), 4*PI/5, 2*PI/5,-PI/6,  -3.0, sqrt(3)/2,-1,1 }, //North Indian O.
				{ atan(1/2.), 4*PI/5, 4*PI/5,-PI/6,  -2.5, sqrt(3),  -1,1 }, //South Africa
				{-atan(1/2.),-PI,	 PI/5,  -PI/6,  -2.5, 0.0,		   -5,5 }, //Australia
				{-atan(1/2.),-3*PI/5, PI,	  PI/2,  -1.5, 0.0,		   -6,4 }, //China
				{-atan(1/2.),-PI/5,	 PI,	  PI/2,  -0.5, 0.0,		   -5,5 }, //North America
				{-atan(1/2.), PI/5,	-3*PI/5, 5*PI/6, 0.5, 0.0,		   -5,5 }, //East Pacific O.
				{-atan(1/2.), 3*PI/5, PI,	  PI/2,   1.5, 0.0,		   -3,3 }, //West Antarctica
				{-atan(1/2.), 3*PI/5,-PI/5,   5*PI/6, 1.0,-sqrt(3)/2,-1,1 }, //South Pacific O.
				{-atan(1/2.), 3*PI/5, PI/5,   PI/6,  -3.0,-sqrt(3)/2,-1,1 }, //New Zealand
				{-PI/2,    0.0,		-PI,	  PI/2,   0.0,-sqrt(3)/2,-3,1 }, //Hawai`i
				{-PI/2,    0.0,		-3*PI/5, PI/2,  -1.0,-sqrt(3)/2,-1,2 }, //West Pacific O.
				{-PI/2,    0.0,		-PI/5,	  PI/2,  -2.0,-sqrt(3)/2, 0,3 }}); //Melanesia
		/*		  LATITUDE,	    LONGITUDE,	 CTR_MERID,	  PLANE_ROT,   X,   Y               RANGE */
		public final int sphereSym, planarSym; //the numbers of symmetries in the two coordinate systems
		public final BoundingBox bounds; //the width and height of a map with this configuration
		public final double[][] centrumSet; //the mathematical information about this configuration
		public final Type type; //holds the number of faces
		
		private Configuration(int sphereSym, int planarSym, BoundingBox bounds, double[][] centrumSet) {
			this.bounds = bounds;
			this.sphereSym = sphereSym;
			this.planarSym = planarSym;
			this.centrumSet = centrumSet;
			if (sphereSym == 3)
				this.type = Type.TETRAHEDRAL;
			else
				this.type = Type.ICOSOHEDRAL;
		}
		
		public double[] rotateOOB(double x, double y, double xCen, double yCen) { //move points that are out of bounds for project()
			return new double[] {x, y}; //this method should be overridden by projections with weird geometry
		}
		
		public boolean inBounds(double x, double y) {return true;} //determine whether a point is in bounds for inverse()
	}
}
