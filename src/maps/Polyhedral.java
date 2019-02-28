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
import utils.Dixon;
import utils.Math2;
import utils.NumericalAnalysis;

/**
 * Projections created by projecting onto and then unfolding some kind of polyhedron
 * 
 * @author jkunimune
 */
public class Polyhedral {
	
	private static final double ASIN_ONE_THD = Math.asin(1/3.); //the complement of the angular radius of a tetrahedron face
	private static final double ATAN_ONE_HLF = Math.atan(1/2.); //the complement of the angular length of an icosahedron edge
	
	
	public static final PolyhedralProjection LEE_TETRAHEDRAL_RECTANGULAR = new PolyhedralProjection(
			"Lee Tetrahedral", 0b1001, Configuration.TETRAHEDRON_WIDE_FACE, Property.CONFORMAL,
			4, null, "that really deserves more attention") {
		
		public double[] faceProject(double lat, double lon) {
			final de.jtem.mfc.field.Complex z = de.jtem.mfc.field.Complex.fromPolar(
					Math.pow(2, 5/6.)*Math.tan(Math.PI/4-lat/2), lon);
			final de.jtem.mfc.field.Complex w = Dixon.invFunc(z);
			return new double[] { w.abs()*2/Dixon.PERIOD_THIRD, w.arg() }; //I don't understand Dixon functions well enough to say whence the 1.132 comes
		}
		
		public double[] faceInverse(double r, double tht) {
			final de.jtem.mfc.field.Complex w = de.jtem.mfc.field.Complex.fromPolar(
					r*Dixon.PERIOD_THIRD/2, tht);
			final de.jtem.mfc.field.Complex ans = Dixon.leeFunc(w).times(Math.pow(2, -5/6.));
			return new double[] {
					Math.PI/2 - 2*Math.atan(ans.abs()),
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
					Math.atan(1/Math.tan(lat)*Math.cos(lon))/Math.cos(lon) /Math.atan(Math.sqrt(2)),
					lon };
		}
		
		public double[] faceInverse(double r, double tht) {
			return new double[] {
					Math.PI/2 - Math.atan(Math.tan(r*Math.cos(tht)*Math.atan(Math.sqrt(2)))/Math.cos(tht)),
					tht };
		}
	};
	
	
	public static final PolyhedralProjection AUTHAGRAPH = new PolyhedralProjection(
			"IMAGO (AuthaGraph)", "Authagraph is a hip new Japanese map that is almost equal-area and would be super great if "
					+ "they actually published their equations. This is technically just an approximation, also known as the Infinitessimal "
					+ "Mutated AuthaGraph Offspring.",
			0b1011, Configuration.AUTHAGRAPH, Property.COMPROMISE, 3,
			new String[] {"Power"}, new double[][] {{.5,1,.68}}) {
		
		private final double[] POLE = {Math.toRadians(77), Math.toRadians(143), Math.toRadians(17)};
		private double k;
		
		public void setParameters(double... params) {
				this.k = params[0];
		}
		
		@Override
		public double[] project(double lat, double lon) { //apply a pole shift to AuthaGraph
			double[] relCoords = obliquifySphc(lat, lon, POLE);
			return super.project(relCoords[0], relCoords[1]);
		}
		
		@Override
		public double[] inverse(double x, double y) { //because AuthaGraph needs its obliquity, and I didn't want to program that into the Configuration
			return obliquifyPlnr(super.inverse(x, y), POLE);
		}
		
		
		public double[] faceProject(double lat, double lon) {
			final double tht = Math.atan((lon - Math.asin(Math.sin(lon)/Math.sqrt(3)))/Math.PI*Math.sqrt(12));
			final double p = (Math.PI/2 - lat) / Math.atan(Math.sqrt(2)/Math.cos(lon));
			return new double[] { Math.pow(p,k)*Math.sqrt(3)/Math.cos(tht), tht };
		}
		
		protected double[] faceInverse(double r, double th) {
			final double lon = NumericalAnalysis.newtonRaphsonApproximation(th, th*2,
					(l) -> Math.atan((l - Math.asin(Math.sin(l)/Math.sqrt(3)))/Math.PI*Math.sqrt(12)),
					(l) -> (1-1/Math.sqrt(1+2*Math.pow(Math.cos(l),-2)))/Math.sqrt(Math.pow(Math.PI,2)/12+Math.pow(l-Math.asin(Math.sin(l)/Math.sqrt(3)),2)),
					.001);
			final double R = r / (Math.sqrt(3)/Math.cos(th));
			return new double[] {
					Math.PI/2 - Math.pow(R,1/k)*Math.atan(Math.sqrt(2)/Math.cos(lon)), lon };
		}
	};
	
	
	public static final PolyhedralProjection AUTHAPOWER = new PolyhedralProjection(
			"TetraPower", "A parametrised, simplified version of my AuthaGraph approximation.",
			0b1011, Configuration.TETRAHEDRON_WIDE_VERTEX, Property.COMPROMISE, 4,
			new String[] {"Power"}, new double[][] {{.5,1,.6}}) {
		
		private double k;
		
		public void setParameters(double... params) {
				this.k = params[0];
		}
		
		public double[] faceProject(double lat, double lon) {
			final double tht = Math.atan((lon - Math.asin(Math.sin(lon)/Math.sqrt(3)))/Math.PI*Math.sqrt(12));
			final double p = (Math.PI/2 - lat) / Math.atan(Math.sqrt(2)/Math.cos(lon));
			return new double[] { Math.pow(p,k)*Math.sqrt(3)/Math.cos(tht), tht };
		}
		
		protected double[] faceInverse(double r, double th) {
			final double lon = NumericalAnalysis.newtonRaphsonApproximation(th, th*2,
					(l) -> Math.atan((l - Math.asin(Math.sin(l)/Math.sqrt(3)))/Math.PI*Math.sqrt(12)),
					(l) -> (1-1/Math.sqrt(1+2*Math.pow(Math.cos(l),-2)))/Math.sqrt(Math.pow(Math.PI,2)/12+Math.pow(l-Math.asin(Math.sin(l)/Math.sqrt(3)),2)),
					.001);
			final double R = r / (Math.sqrt(3)/Math.cos(th));
			return new double[] {
					Math.PI/2 - Math.pow(R,1/k)*Math.atan(Math.sqrt(2)/Math.cos(lon)), lon };
		}
	};
	
	
	public static final PolyhedralProjection ACTUAUTHAGRAPH = new PolyhedralProjection(
			"EquaHedral", "An interrupted authalic tetrahedral projection.",
			0b1010, Configuration.TETRAHEDRON_WIDE_VERTEX, Property.EQUAL_AREA, 3,
			new String[] {"Sinus length"}, new double[][] {{0, 60, 20}}) {
		
		private double sig, a0, scale;
		
		public void setParameters(double... params) {
			this.sig = params[0]/60;
			this.a0 = 3 - 1.5*sig*sig;
			this.scale = Math.sqrt(3)*a0/Math.PI;
		}
		
		public double[] faceProject(double lat, double lon) {
			double bet = Math.atan((lon-Math.asin(Math.sin(lon)/Math.sqrt(3)))/(a0/2)*scale);
			double f = (1-Math.sin(lat))/(1-1/Math.sqrt(1+2/Math2.cos2(lon)));
			if (f < (3*sig*sig)/(2*a0)) { //sinus zone
				double alf = Math.atan(2*Math.tan(Math.abs(bet)) - 1/Math.sqrt(3));
				double rA = Math.sqrt(2*f*a0/Math2.cos2(alf))/2;
				return toPolar(rA, alf, Math.signum(lon));
			}
			else { //primary zone
				double rB = Math.sqrt(a0*f - a0 + 3)/Math.cos(bet);
				return new double[] {rB, bet};
			}
		}
		
		public double[] faceInverse(double r, double th) {
			double bet;
			double f;
			if (r < sig*Math.sqrt(3)/2/Math.cos(Math.abs(th)-Math.PI/3)) //empty
				return null;
			else if (r*Math.cos(th) < sig*Math.sqrt(3)) { //sinus zone
				double[] relCoords = fromPolar(r, th);
				double rA = relCoords[0];
				double alf = relCoords[1];
				bet = Math.signum(th)*Math.atan(Math.tan(alf)/2+1/Math.sqrt(12));
				f = Math.pow(rA, 2) * 2*Math2.cos2(alf)/a0;
			}
			else { //primary zone
				bet = th;
				f = (r*r*Math2.cos2(bet) - 1.5*sig*sig)/a0;
			}
			double lon = NumericalAnalysis.newtonRaphsonApproximation(
					a0/2*Math.tan(bet)/scale, bet*2,
					(l) -> l - Math.asin(Math.sin(l)/Math.sqrt(3)),
					(l) -> 1 - 1/Math.sqrt(1 + 2/Math2.cos2(l)), 1e-4);
			double lat = Math.asin(1 - f*(1 - 1/Math.sqrt(1+2/Math2.cos2(lon))));
			return new double[] {lat, lon};
		}
		
		private double[] toPolar(double rA, double alf, double s) {
			double x = rA*Math.cos(alf) + sig*Math.sqrt(3)/2;
			double y = rA*Math.sin(alf) + sig/2;
			return new double[] {Math.hypot(x, y), s*Math.atan2(y, x)};
		}
		
		private double[] fromPolar(double rB, double bet) {
			double x = rB*Math.cos(bet) - sig*Math.sqrt(3)/2;
			double y = Math.abs(rB*Math.sin(bet)) - sig/2;
			return new double[] {Math.hypot(x, y), Math.atan2(y, x)};
		}
	};
	
	
	public static final Projection VAN_LEEUWEN = new PolyhedralProjection(
			"van Leeuwen", "An uninterrupted equal-area tetrahedral projection. It's more accurately known as \"the Vertex-oriented great circle projection applied to a tetrahedron\", but the guy who copublished it with Leeuwen calls it \"the van Leeuwen projection\" on his website, so I think this is fine.",
			0b1011, Configuration.TETRAHEDRON_WIDE_VERTEX, Property.EQUAL_AREA, 3) {
		
		public double[] faceProject(double lat, double lon) {
			ACTUAUTHAGRAPH.setParameters(0);
			return ACTUAUTHAGRAPH.faceProject(lat, lon);
		}
		
		public double[] faceInverse(double r, double th) {
			ACTUAUTHAGRAPH.setParameters(0);
			return ACTUAUTHAGRAPH.faceInverse(r, th);
		}
	};
	
	
	public static final Projection DYMAXION = new PolyhedralProjection(
			"Dymaxion", "A polyhedral projection that slices up the oceans as much as possible without slicing up any landmasses.",
			0b1110, Configuration.DYMAXION, Property.COMPROMISE, 3) {
		
		private final double[] POLE = {0.040158, -0.091549,-2.015269}; //I derived these numbers from [Robert Gray](http://www.rwgrayprojects.com/rbfnotes/maps/graymap4.html)
		private final double X_0 = 0.75;
		private final double Y_0 = -Math.sqrt(3)/4;
		
		private final double sin36 = Math.sqrt(10-2*Math.sqrt(5))/4;
		private final double cos36 = (1+Math.sqrt(5))/4;
		
		@Override
		public double[] project(double lat, double lon) { //apply a pole shift and Cartesian shift to Dymaxion
			double[] coords = obliquifySphc(lat, lon, POLE);
			coords = super.project(coords[0], coords[1]);
			return new double[] {coords[0] + X_0, coords[1] + Y_0};
		}
		
		@Override
		public double[] inverse(double x, double y) { //because Dymaxion needs its obliquity, and I didn't want to program that into the Configuration
			double[] coords = super.inverse(x - X_0, y - Y_0);
			if (coords == null) 	return null;
			return obliquifyPlnr(coords, POLE);
		}
		
		public double[] faceProject(double lat, double lon) {
			double xG = Math.cos(lon)/Math.tan(lat)/cos36; //normalised gnomonic coordinates
			double yG = Math.sin(lon)/Math.tan(lat)/sin36;
			double a = Math.asin((xG+yG)/(2*Math.sqrt(1+xG*xG))) + Math.atan(xG); //angular distance up each side of the triangle
			double b = Math.asin((xG-yG)/(2*Math.sqrt(1+xG*xG))) + Math.atan(xG);
			double x = (a + b)/(2*Math.sqrt(3)); //final cartesian coordinates in radians
			double y = (a - b)/2;
			return new double[] {Math.hypot(x,y)/Math.atan(2), Math.atan2(y,x)}; //scale to fit to layout, where side length is 1
		}
		
		public double[] faceInverse(double r, double th) {
			if (Math.abs(th) > Math.PI/6) 	throw new IllegalArgumentException("Wait, what?");
			double x = r*Math.cos(th)*Math.atan(2); //cartesian coordinates in radians
			double y = r*Math.sin(th)*Math.atan(2);
			double a = Math.sqrt(3)*x + y; //angular distance up each side of the triangle
			double b = Math.sqrt(3)*x - y;
			double xG = cos36*(Math.sin(a) + Math.sin(b))/(1 + Math.cos(a) + Math.cos(b)); //unnormalised gnomonic coordinates
			double yG = sin36*
					(Math.sin(a) - Math.sin(b) + 2*Math.sin(a-b))/(1 + Math.cos(a) + Math.cos(b));
			return new double[] {Math.atan(1/Math.hypot(xG, yG)), Math.atan2(yG, xG)}; //inverse gnomonic projection
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
			super(name, config.width, config.height, fisc, config.type, property, rating,
					adjective, addendum);
			this.configuration = config;
		}
		
		public PolyhedralProjection(
				String name, String description, int fisc, Configuration config, Property property,
				int rating) {
			super(name, description, config.width, config.height, fisc, config.type, property,
					rating);
			this.configuration = config;
		}
		
		public PolyhedralProjection(
				String name, String description, int fisc, Configuration config, Property property,
				int rating, String[] paramNames, double[][] paramValues) {
			super(name, description, config.width, config.height, fisc, config.type, property,
					rating, paramNames, paramValues);
			this.configuration = config;
		}
		
		
		protected abstract double[] faceProject(double lat, double lon); //the projection from spherical to polar within a face
		
		protected abstract double[] faceInverse(double x, double y); //I think you can guess
		
		
		public double[] project(double lat, double lon) {
			final int numSym = configuration.sphereSym; //we're about to be using this variable a lot
			double latR = Double.NEGATIVE_INFINITY;
			double lonR = Double.NEGATIVE_INFINITY;
			double[] centrum = null;
			for (double[] testCentrum: configuration.centrumSet) { //iterate through the centrums to see which goes here
				final double[] relCoords = obliquifySphc(lat, lon, testCentrum);
				if (testCentrum.length > 6) { //if the centrum is long, then it contains longitude bounds
					double minL = testCentrum[6]*Math.PI/numSym;
					double maxL = testCentrum[7]*Math.PI/numSym;
					relCoords[1] = Math2.floorMod(relCoords[1]-minL, 2*Math.PI) + minL;
					if (relCoords[1] < minL || relCoords[1] > maxL)
						continue; //ignore any longitudes not in the bounds described in [6:7]
				}
				
				if (relCoords[0] > latR) { //pick the centrum that maxes out latitude
					latR = relCoords[0];
					lonR = relCoords[1];
					centrum = testCentrum;
				}
			}
			
			final double lonR0 = Math.floor((lonR+Math.PI/numSym)/(2*Math.PI/numSym))
					*(2*Math.PI/numSym); //because most face projections are periodic
			
			final double[] rth = faceProject(latR, lonR - lonR0); //apply the projection to the relative coordinates
			final double r = rth[0];
			final double th = rth[1] + centrum[3] + lonR0*numSym/configuration.planarSym;
			final double x0 = centrum[4];
			final double y0 = centrum[5];
			
			double[] output = { r*Math.cos(th) + x0, r*Math.sin(th) + y0 };
			if (Math.abs(output[0]) > width/2 || Math.abs(output[1]) > height/2) { //rotate OOB bits around nearest singularity
				output = configuration.rotateOOB(output[0], output[1], x0, y0);
			}
			return output;
		}
		
		
		public double[] inverse(double x, double y) {
			if (!configuration.inBounds(x, y)) 	return null;
			
			final int numSym = configuration.planarSym; //we'll be using this variable a lot soon
			
			double rM = Double.POSITIVE_INFINITY;
			double[] centrum = null; //iterate to see which centrum we get
			for (double[] testCentrum: configuration.centrumSet) {
				final double rR = Math.hypot(x-testCentrum[4], y-testCentrum[5]);
				if (rR < rM) { //pick the centrum that minimises r
					rM = rR;
					centrum = testCentrum;
				}
			}
			
			final double th0 = centrum[3];
			final double x0 = centrum[4];
			final double y0 = centrum[5];
			final double r = Math.hypot(x - x0, y - y0);
			final double th = Math2.coerceAngle(Math.atan2(y - y0, x - x0) - th0);
			
			if (centrum.length > 6) { //if the centrum has extra values, they are angle bounds
				if (th < centrum[6]*Math.PI/numSym || th > centrum[7]*Math.PI/numSym)
					return null; //ignore any angles not in the bounds described in [6:7]
			}
			
			final double thBase = Math.floor((th+Math.PI/numSym)/(2*Math.PI/numSym))
					*(2*Math.PI/numSym); //because most face projections are periodic
			
			double[] relCoords = faceInverse(r, th - thBase);
			
			if (relCoords == null)
				return null;
			
			relCoords[1] = thBase*numSym/configuration.sphereSym + relCoords[1];
			double[] absCoords = obliquifyPlnr(relCoords, centrum);
			if (Math.abs(absCoords[1]) > Math.PI)
				absCoords[1] = Math2.coerceAngle(absCoords[1]);
			return absCoords;
		}
	}
	
	
	
	/**
	 * A set of objects that determine the layouts of tetrahedral projections
	 * 
	 * @author jkunimune
	 */
	private static enum Configuration {
		
		/*		  LATITUDE,		 LONGITUDE,		 CTR_MERID,		 PLANE_ROT,		 X,	 Y */
		TETRAHEDRON_WIDE_FACE(3, 3, 6., 2*Math.sqrt(3), new double[][] { // [<|>] arrangement, face-centred
				{ ASIN_ONE_THD,	 Math.PI,	-2*Math.PI/3,	-2*Math.PI/3,	 2,	 Math.sqrt(3),-1,2 },
				{ ASIN_ONE_THD,	 Math.PI,	 2*Math.PI/3,	-Math.PI/3,		-2,	 Math.sqrt(3) },
				{-Math.PI/2,	 0,			 2*Math.PI/3,	 2*Math.PI/3,	 2,	-Math.sqrt(3),-2,1 },
				{-Math.PI/2,	 0,			-2*Math.PI/3,	 Math.PI/3,		-2,	-Math.sqrt(3) },
				{ ASIN_ONE_THD,	 Math.PI/3,	-2*Math.PI/3,	 Math.PI,		 1,	 0 },
				{ ASIN_ONE_THD,	-Math.PI/3,	 2*Math.PI/3,	 0,				-1,	 0 }}),
		
		TRIANGLE_FACE(3, 3, 4*Math.sqrt(3), 6., new double[][] { // \delta arrangement, like they are often published
				{ ASIN_ONE_THD,	 Math.PI/3,	 0,				-5*Math.PI/6,	 Math.sqrt(3),	2 },
				{ ASIN_ONE_THD,	-Math.PI/3,	 0,				-Math.PI/6,		-Math.sqrt(3),	2 },
				{ ASIN_ONE_THD,	 Math.PI,	 0,				 Math.PI/2,		 0,			-1 },
				{-Math.PI/2,	 0,			 0,				-Math.PI/2,		 0,			 1 }}) {
			@Override public boolean inBounds(double x, double y) {
				return y > Math.sqrt(3)*Math.abs(x) - 3;
			}
		},
		TETRAHEDRON_WIDE_VERTEX(3, 6, 6., 2*Math.sqrt(3), new double[][] { // [<|>] arrangement, vertex-centred
				{ Math.PI/2,	 0,				 0,				-Math.PI/2,		 0,	 Math.sqrt(3) },
				{-ASIN_ONE_THD,	 0,				 Math.PI,		 Math.PI/2,		 0,	-Math.sqrt(3) },
				{-ASIN_ONE_THD,	 2*Math.PI/3,	 Math.PI,		 5*Math.PI/6,	 3,	 0 },
				{-ASIN_ONE_THD,	-2*Math.PI/3,	 Math.PI,		 Math.PI/6,		-3,	 0 }}) {
			@Override public double[] rotateOOB(double x, double y, double xCen, double yCen) {
				if (Math.abs(x) > width/2)
					return new double[] {2*xCen - x, -y};
				else
					return new double[] {-x, height*Math.signum(y) - y};
			}
		},
		AUTHAGRAPH(3, 6, 4*Math.sqrt(3), 3, new double[][] { // |\/\/`| arrangement, vertex-centred
				{-ASIN_ONE_THD,	 Math.PI,		 Math.PI,	 0,	-2*Math.sqrt(3)-.6096,	 1.5 },
				{-ASIN_ONE_THD,	-Math.PI/3,		 Math.PI/3,	 0,	-Math.sqrt(3)-.6096,	-1.5 },
				{ Math.PI/2,	 0,				 Math.PI,	 0,	 0-.6096,				 1.5 },
				{-ASIN_ONE_THD,	 Math.PI/3,		-Math.PI/3,	 0,	 Math.sqrt(3)-.6096,	-1.5 },
				{-ASIN_ONE_THD,	 Math.PI,		 Math.PI,	 0,	 2*Math.sqrt(3)-.6096,	 1.5 },
				{-ASIN_ONE_THD,	-Math.PI/3,		 Math.PI/3,	 0,	 3*Math.sqrt(3)-.6096,	-1.5 }}) {
			@Override public double[] rotateOOB(double x, double y, double xCen, double yCen) {
				if (Math.abs(y) > height/2) {
					x = 2*xCen - x;
					y = 2*yCen - y;
				}
				if (Math.abs(x) > width/2)
					x = Math2.floorMod(x+width/2,width)-width/2;
				return new double[] {x, y};
			}
		},
		DYMAXION(5, 6, 5.5, 1.5*Math.sqrt(3), new double[][] { // I can't draw this in ASCII. You know what "Dymaxion" means
				{ Math.PI/2,    0.0,		-3*Math.PI/5,-Math.PI/2,  -1.5, Math.sqrt(3),  -3,3 }, //West Africa
				{ Math.PI/2,    0.0,		 Math.PI/5,	 -Math.PI/2,   0.5, Math.sqrt(3),  -1,1 }, //Brazil
				{ Math.PI/2,    0.0,		 3*Math.PI/5,-Math.PI/2,   1.5, Math.sqrt(3),  -1,1 }, //South Atlantic O.
				{ ATAN_ONE_HLF,-4*Math.PI/5, 2*Math.PI/5,-Math.PI/6,  -2.0, Math.sqrt(3)/2,-5,5 }, //Arabia
				{ ATAN_ONE_HLF,-2*Math.PI/5,-2*Math.PI/5,-5*Math.PI/6,-1.0, Math.sqrt(3)/2,-5,5 }, //Scandanavia
				{ ATAN_ONE_HLF, 0.0,		 0.0,		 -Math.PI/2,   0.0, Math.sqrt(3)/2,-3,5 }, //Caribbean
				{ ATAN_ONE_HLF, 0.0,		-4*Math.PI/5,-5*Math.PI/6,-0.5, Math.sqrt(3),  -1,1 }, //North Atlantic O.
				{ ATAN_ONE_HLF, 2*Math.PI/5, 0.0,		 -Math.PI/2,   1.0, Math.sqrt(3)/2,-5,5 }, //Patagonia
				{ ATAN_ONE_HLF, 4*Math.PI/5,-2*Math.PI/5,-5*Math.PI/6, 2.0, Math.sqrt(3)/2,-3,2 }, //East Antarctica
				{ ATAN_ONE_HLF, 4*Math.PI/5, 0.0,		 -Math.PI/6,  -3.5, 0.0,		    0,1 }, //South Indian O.
				{ ATAN_ONE_HLF, 4*Math.PI/5, 2*Math.PI/5,-Math.PI/6,  -3.0, Math.sqrt(3)/2,-1,1 }, //North Indian O.
				{ ATAN_ONE_HLF, 4*Math.PI/5, 4*Math.PI/5,-Math.PI/6,  -2.5, Math.sqrt(3),  -1,1 }, //South Africa
				{-ATAN_ONE_HLF,-Math.PI,	 Math.PI/5,  -Math.PI/6,  -2.5, 0.0,		   -5,5 }, //Australia
				{-ATAN_ONE_HLF,-3*Math.PI/5, Math.PI,	  Math.PI/2,  -1.5, 0.0,		   -6,4 }, //China
				{-ATAN_ONE_HLF,-Math.PI/5,	 Math.PI,	  Math.PI/2,  -0.5, 0.0,		   -5,5 }, //North America
				{-ATAN_ONE_HLF, Math.PI/5,	-3*Math.PI/5, 5*Math.PI/6, 0.5, 0.0,		   -5,5 }, //East Pacific O.
				{-ATAN_ONE_HLF, 3*Math.PI/5, Math.PI,	  Math.PI/2,   1.5, 0.0,		   -3,3 }, //West Antarctica
				{-ATAN_ONE_HLF, 3*Math.PI/5,-Math.PI/5,   5*Math.PI/6, 1.0,-Math.sqrt(3)/2,-1,1 }, //South Pacific O.
				{-ATAN_ONE_HLF, 3*Math.PI/5, Math.PI/5,   Math.PI/6,  -3.0,-Math.sqrt(3)/2,-1,1 }, //New Zealand
				{-Math.PI/2,    0.0,		-Math.PI,	  Math.PI/2,   0.0,-Math.sqrt(3)/2,-3,1 }, //Hawai`i
				{-Math.PI/2,    0.0,		-3*Math.PI/5, Math.PI/2,  -1.0,-Math.sqrt(3)/2,-1,2 }, //West Pacific O.
				{-Math.PI/2,    0.0,		-Math.PI/5,	  Math.PI/2,  -2.0,-Math.sqrt(3)/2, 0,3 }}); //Melanesia
		/*		  LATITUDE,	    LONGITUDE,	 CTR_MERID,	  PLANE_ROT,   X,   Y               RANGE */
		public final int sphereSym, planarSym; //the numbers of symmetries in the two coordinate systems
		public final double width, height; //the width and height of a map with this configuration
		public final double[][] centrumSet; //the mathematical information about this configuration
		public final Type type; //holds the number of faces
		
		private Configuration(int sphereSym, int planarSym, double width, double height, double[][] centrumSet) {
			this.width = width;
			this.height = height;
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
