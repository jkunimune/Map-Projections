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
import utils.NumericalAnalysis;
import utils.Shape;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.NaN;
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
 * Projections created by projecting onto and then unfolding some kind of polyhedron – specifically a tetrahedron or an
 * icosohedron (octahedrons are in a different file)
 * 
 * @author jkunimune
 */
public class Polyhedral {
	
	public static final PolyhedralProjection LEE_TETRAHEDRAL_RECTANGULAR = new PolyhedralProjection(
			"Lee Tetrahedral", "A conformal tetrahedral projection that really deserves more attention",
			true, false, false, Polyhedron.TETRAHEDRON_WIDE_FACE, Property.CONFORMAL, 4) {
		
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
			"Lee Tetrahedral (triangular)", "A conformal tetrahedral projection in a triangle, because this is the form in which it was published, even though the rectangle is clearly better",
			true, false, false, Polyhedron.TRIANGLE_FACE, Property.CONFORMAL, 2) {
		
		public double[] faceProject(double lat, double lon) {
			return LEE_TETRAHEDRAL_RECTANGULAR.faceProject(lat, lon);
		}
		
		public double[] faceInverse(double r, double tht) {
			return LEE_TETRAHEDRAL_RECTANGULAR.faceInverse(r, tht);
		}
	};
	
	
	public static final PolyhedralProjection TETRAGRAPH = new PolyhedralProjection(
			"TetraGraph", "An equidistant tetrahedral projection that I invented",
			true, true, true, Polyhedron.TETRAHEDRON_WIDE_FACE, Property.EQUIDISTANT, 2) {
		
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
			"IMAGO (AuthaGraph)",
			"AuthaGraph is a hip new Japanese map that would be super great if they actually " +
			"published their equations. This is technically just an approximation, also known as " +
			"the Infinitesimal Mutated AuthaGraph Offspring.",
			true, true, false, Polyhedron.AUTHAGRAPH, Property.COMPROMISE, 3,
			new String[] {"Power"}, new double[][] {{.5,1,.68}}) {
		
		private final double[] POLE = {toRadians(77), toRadians(143), toRadians(17)};
		
		@Override
		public void initialize(double... params) {
			AUTHAGRAPH_ALT.initialize(params[0]);
		}
		
		@Override
		public double[] project(double lat, double lon) { //apply a pole shift to AuthaGraph
			double[] relCoords = transformFromOblique(lat, lon, POLE);
			return super.project(relCoords[0], relCoords[1]);
		}
		
		@Override
		public double[] inverse(double x, double y) { //because AuthaGraph needs its obliquity, and I didn't want to program that into the Polyhedron
			return transformToOblique(super.inverse(x, y), POLE);
		}
		
		public double[] faceProject(double lat, double lon) {
			return AUTHAGRAPH_ALT.faceProject(lat, lon);
		}
		
		public double[] faceInverse(double r, double th) {
			return AUTHAGRAPH_ALT.faceInverse(r, th);
		}
	};
	
	
	public static final PolyhedralProjection AUTHAGRAPH_ALT = new PolyhedralProjection(
			"IMAGO (simplified)", "An approximation of the AuthaGraph projection, rearranged to an alternate layout that's more north-up",
			true, true, false, Polyhedron.TETRAHEDRON_WIDE_VERTEX, Property.COMPROMISE, 4,
			new String[] {"Power"}, new double[][] {{.5,1,.68}}) {
		
		private double k;
		
		public void initialize(double... params) {
				this.k = params[0];
		}
		
		public double[] faceProject(double lat, double lon) {
			final double tht = atan((lon - asin(sin(lon)/sqrt(3)))/PI*sqrt(12));
			final double p = (PI/2 - lat) / atan(sqrt(2)/cos(lon));
			return new double[] { pow(p,k)*sqrt(3)/cos(tht), tht };
		}
		
		public double[] faceInverse(double r, double th) {
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
			"EquaHedral", "An interrupted authalic tetrahedral projection",
			true, true, false, Polyhedron.TETRAHEDRON_WIDE_VERTEX, Property.EQUAL_AREA, 3,
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
			"Van Leeuwen",
			"An uninterrupted equal-area tetrahedral projection. It's more accurately known as " +
			"\"the Vertex-oriented great circle projection applied to a tetrahedron\", but the " +
			"guy who copublished it with Leeuwen calls it \"the van Leeuwen projection\" on his " +
			"website, so I think this is fine.",
			true, true, false, Polyhedron.TETRAHEDRON_WIDE_VERTEX, Property.EQUAL_AREA, 2) {
		
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
			"Dymaxion", "A polyhedral projection that slices up the oceans as much as possible without slicing up any landmasses",
			true, true, true, Polyhedron.DYMAXION, Property.COMPROMISE, 3) {
		
		private final double[] POLE = {0.040158, -0.091549,-2.015269}; //I derived these numbers from [Robert Gray](http://www.rwgrayprojects.com/rbfnotes/maps/graymap4.html)

		private final double sin36 = sqrt(10-2*sqrt(5))/4;
		private final double cos36 = (1+sqrt(5))/4;
		
		@Override
		public double[] project(double lat, double lon) { //apply a pole shift to Dymaxion
			double[] coords = transformFromOblique(lat, lon, POLE);
			return super.project(coords[0], coords[1]);
		}
		
		@Override
		public double[] inverse(double x, double y) { //because Dymaxion needs its obliquity, and I didn't want to program that into the Polyhedron
			double[] coords = super.inverse(x, y);
			if (coords == null)
				return null;
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
			if (abs(th) > PI/6+1e-15)
				throw new IllegalArgumentException("Wait, what?"+th);
			double x = r*cos(th)*atan(2); //cartesian coordinates in radians
			double y = r*sin(th)*atan(2);
			double a = sqrt(3)*x + y; //angular distance up each side of the triangle
			double b = sqrt(3)*x - y;
			double correction = 1 + cos(a) + cos(b);
			double xG = cos36*(sin(a) + sin(b))/correction; //unnormalised gnomonic coordinates
			double yG = sin36*(sin(a) - sin(b) + 2*sin(a-b))/correction;
			return new double[] {atan(1/hypot(xG, yG)), atan2(yG, xG)}; //inverse gnomonic projection
		}
	};
	
	
	
	/**
	 * A base for polyhedral Projections
	 * 
	 * @author jkunimune
	 */
	private static abstract class PolyhedralProjection extends Projection {
		
		private final Polyhedron configuration;
		
		
		public PolyhedralProjection(
				String name, String description, boolean comprehensive, boolean solvable,
				boolean invertible, Polyhedron config, Property property, int rating) {
			super(name, description, config.shape, false, comprehensive, solvable, invertible,
			      config.type, property, rating);
			this.configuration = config;
		}
		
		public PolyhedralProjection(
				String name, String description, boolean comprehensive, boolean solvable,
				boolean invertible, Polyhedron config, Property property, int rating,
				String[] paramNames, double[][] paramValues) {
			super(name, description, config.shape, false, comprehensive, solvable, invertible,
			      config.type, property, rating, paramNames, paramValues);
			this.configuration = config;
		}
		
		
		protected abstract double[] faceProject(double lat, double lon); //the projection from spherical to polar within a face
		
		protected abstract double[] faceInverse(double x, double y); //I think you can guess
		
		
		public double[] project(double lat, double lon) {
			final int numSym = configuration.sphereSym; //we're about to be using this variable a lot
			double latR = NEGATIVE_INFINITY;
			double lonR = NEGATIVE_INFINITY;
			Facet facet = null;
			for (Facet testFacet: configuration.facets) { //iterate through the facets to see which applies here
				final double[] relCoords = transformFromOblique(lat, lon, new double[] {
						testFacet.latitude, testFacet.longitude, testFacet.centralMeridian});
				relCoords[1] = floorMod(relCoords[1] - testFacet.minAngle, 2*PI) + testFacet.minAngle;
				if (relCoords[1] < testFacet.minAngle || relCoords[1] > testFacet.maxAngle)
					continue; //ignore any longitudes not in the bounds

				if (relCoords[0] > latR) { //pick the facet that maximizes latitude
					latR = relCoords[0];
					lonR = relCoords[1];
					facet = testFacet;
				}
			}
			if (facet == null)
				return new double[] {NaN, NaN}; // return NaN if no facet appears to contain this point (should never happen)
			
			final double lonR0 = floor((lonR+PI/numSym)/(2*PI/numSym))
					*(2*PI/numSym); //because most face projections are periodic
			
			final double[] rth = faceProject(latR, lonR - lonR0); //apply the projection to the relative coordinates
			final double r = rth[0];
			final double th = rth[1] + facet.planeRotation + lonR0*numSym/configuration.planarSym;
			final double x0 = facet.x;
			final double y0 = facet.y;
			
			double[] output = { r*cos(th) + x0, r*sin(th) + y0 };
			if (output[0] < shape.xMin || output[0] > shape.xMax || output[1] < shape.yMin || output[1] > shape.yMax) { //rotate OOB bits around nearest singularity
				output = configuration.rotateOOB(output[0], output[1], x0, y0);
			}
			return output;
		}
		
		
		public double[] inverse(double x, double y) {
			if (!configuration.inBounds(x, y))
				return null;
			
			final int numSym = configuration.planarSym; //we'll be using this variable a lot soon
			
			double rM = POSITIVE_INFINITY;
			Facet facet = null; //iterate to see which facet we get
			for (Facet testFacet: configuration.facets) {
				final double rR = hypot(x-testFacet.x, y-testFacet.y);
				if (rR < rM) { //pick the facet that minimises r
					rM = rR;
					facet = testFacet;
				}
			}
			if (facet == null)
				return null; // return null if no facet appears to contain this point (should never happen)
			
			final double th0 = facet.planeRotation;
			final double x0 = facet.x;
			final double y0 = facet.y;
			final double r = hypot(x - x0, y - y0);
			final double th = coerceAngle(atan2(y - y0, x - x0) - th0);
			
			if (th < facet.minAngle*configuration.sphereSym/configuration.planarSym ||
			    th > facet.maxAngle*configuration.sphereSym/configuration.planarSym)
				return null; //ignore any angles not in the bounds

			final double thBase = floor((th+PI/numSym)/(2*PI/numSym))
					*(2*PI/numSym); //because most face projections are periodic

			double[] relCoords = faceInverse(r, th - thBase);

			if (relCoords == null)
				return null;

			relCoords[1] = thBase*numSym/configuration.sphereSym + relCoords[1];
			double[] absCoords = transformToOblique(relCoords, new double[] {
					facet.latitude, facet.longitude, facet.centralMeridian});
			if (abs(absCoords[1]) > PI)
				absCoords[1] = coerceAngle(absCoords[1]);
			return absCoords;
		}
	}


	/**
	 * an object containing all of the information needed to map a polyhedron to a plane, including the nature of the
	 * polyhedron and the arrangement of its facets in the plane.  not included is how to map the globe to the polyhedron.
	 */
	private enum Polyhedron {
		
		/** [<|>] arrangement, face-centerd */
		TETRAHEDRON_WIDE_FACE(3, 3, Shape.rectangle(6., 2*sqrt(3)), new Facet[] {
				new Facet( 2,  sqrt(3), -2*PI/3, asin(1/3.),  PI  , -2*PI/3,   -PI/3, 2*PI/3),
				new Facet(-2,  sqrt(3),   -PI/3, asin(1/3.),  PI  ,  2*PI/3),
				new Facet( 2, -sqrt(3),  2*PI/3,      -PI/2,   0  ,  2*PI/3, -2*PI/3,   PI/3),
				new Facet( 2, -sqrt(3),  2*PI/3,      -PI/2,   0  ,  2*PI/3, -2*PI/3,   PI/3),
				new Facet(-2, -sqrt(3),    PI/3,      -PI/2,   0  , -2*PI/3),
				new Facet( 1,        0,    PI  , asin(1/3.),  PI/3, -2*PI/3),
				new Facet(-1,        0,     0  , asin(1/3.), -PI/3,  2*PI/3),
				}),
		/** \nabla arrangement, like they are often published */
		TRIANGLE_FACE(3, 3, Shape.polygon(new double[][] {
				{0, -4}, {2*sqrt(3), 2.}, {-2*sqrt(3), 2.}
		}), new Facet[] {
				new Facet( sqrt(3),  1, -5*PI/6, asin(1/3.),  PI/3, 0),
				new Facet(-sqrt(3),  1,   -PI/6, asin(1/3.), -PI/3, 0),
				new Facet(       0, -2,    PI/2, asin(1/3.),  PI  , 0),
				new Facet(       0,  0,   -PI/2,      -PI/2,   0  , 0),
		}) {
			@Override public boolean inBounds(double x, double y) {
				return y > sqrt(3)*abs(x) - 3;
			}
		},
		/** [<|>] arrangement, vertex-centred */
		TETRAHEDRON_WIDE_VERTEX(3, 6, Shape.rectangle(6., 2*sqrt(3)), new Facet[] {
				new Facet( 0,  sqrt(3),  -PI/2,        PI/2,     0  ,  0),
				new Facet( 0, -sqrt(3),   PI/2, -asin(1/3.),     0  , PI),
				new Facet( 3,        0, 5*PI/6, -asin(1/3.),  2*PI/3, PI),
				new Facet(-3,        0,   PI/6, -asin(1/3.), -2*PI/3, PI),
				}) {
			@Override public double[] rotateOOB(double x, double y, double xCen, double yCen) {
				if (x < shape.xMin || x > shape.xMax)
					return new double[] {2*xCen - x, -y};
				else
					return new double[] {-x, 2*shape.yMax*signum(y) - y};
			}
		},
		/** |\/\/`| arrangement, vertex-centred */
		AUTHAGRAPH(3, 6, Shape.rectangle(4*sqrt(3), 3), new Facet[] {
				new Facet(-2*sqrt(3)-.6096,  1.5, 0, -asin(1/3.),  PI  ,  PI  ),
				new Facet(  -sqrt(3)-.6096, -1.5, 0, -asin(1/3.), -PI/3,  PI/3),
				new Facet(          -.6096,  1.5, 0,        PI/2,   0  ,  PI  ),
				new Facet(   sqrt(3)-.6096, -1.5, 0, -asin(1/3.),  PI/3, -PI/3),
				new Facet( 2*sqrt(3)-.6096,  1.5, 0, -asin(1/3.),  PI  ,  PI  ),
				new Facet( 3*sqrt(3)-.6096, -1.5, 0, -asin(1/3.), -PI/3,  PI/3),
				}) {
			@Override public double[] rotateOOB(double x, double y, double xCen, double yCen) {
				if (y < shape.yMin || y > shape.yMax) {
					x = 2*xCen - x;
					y = 2*yCen - y;
				}
				if (x < shape.xMin || x > shape.xMax)
					x = floorMod(x + shape.xMax,2*shape.xMax) - shape.xMax;
				return new double[] {x, y};
			}
		},
		/** I can't draw this in ASCII. You know what "Dymaxion" means */
		DYMAXION(5, 6, Shape.polygon(new double[][] {
				{ 0.0,  0.0        }, { 0.5, -0.5*sqrt(3)}, { 1.5, -0.5*sqrt(3)}, { 1.0,  0.0        },
				{ 2.5,  0.0        }, { 2.5,  0.5*sqrt(3)}, { 2.0,  1.0*sqrt(3)}, { 1.5,  0.5*sqrt(3)},
				{ 1.0,  1.0*sqrt(3)}, { 0.5,  0.5*sqrt(3)}, {-0.5,  0.5*sqrt(3)}, { 0.0,  1.0*sqrt(3)},
				{-2.0,  1.0*sqrt(3)}, {-1.5,  0.5*sqrt(3)}, {-2.5,  0.5*sqrt(3)}, {-2.0,  0.0        },
				{-3.0,  0.0        }, {-2.25,-.25*sqrt(3)}, {-2.5, -0.5*sqrt(3)}, {-1.5, -0.5*sqrt(3)},
				{-1.5,-1/6.*sqrt(3)}, {-1.0,  0.0        }, {-1.0,-1/3.*sqrt(3)}, {-0.5, -0.5*sqrt(3)},
		}), new Facet[] {
				new Facet(-1.0,  1.0*sqrt(3),   -PI/2,        PI/2,     0  , -3*PI/5, -3*PI/5, 3*PI/5), //West Africa
				new Facet( 1.0,  1.0*sqrt(3),   -PI/2,        PI/2,     0  ,    PI/5,   -PI/5,   PI/5), //Brazil
				new Facet( 2.0,  1.0*sqrt(3),   -PI/2,        PI/2,     0  ,  3*PI/5,   -PI/5,   PI/5), //South Atlantic O.
				new Facet(-1.5,  0.5*sqrt(3),   -PI/6,  atan(1/2.), -4*PI/5,  2*PI/5,   -PI  ,   PI  ), //Arabia
				new Facet(-0.5,  0.5*sqrt(3), -5*PI/6,  atan(1/2.), -2*PI/5, -2*PI/5,   -PI  ,   PI  ), //Scandanavia
				new Facet( 0.5,  0.5*sqrt(3),   -PI/2,  atan(1/2.),     0  ,     0  , -3*PI/5,   PI  ), //Caribbean
				new Facet( 0.0,  1.0*sqrt(3), -5*PI/6,  atan(1/2.),     0  , -4*PI/5,   -PI/5,   PI/5), //North Atlantic O.
				new Facet( 1.5,  0.5*sqrt(3),   -PI/2,  atan(1/2.),  2*PI/5,     0  ,   -PI  ,   PI  ), //Patagonia
				new Facet( 2.5,  0.5*sqrt(3), -5*PI/6,  atan(1/2.),  4*PI/5, -2*PI/5, -3*PI/5, 2*PI/5), //East Antarctica
				new Facet(-3.0,  0.0        ,   -PI/6,  atan(1/2.),  4*PI/5,     0  ,     0  ,   PI/5), //South Indian O.
				new Facet(-2.5,  0.5*sqrt(3),   -PI/6,  atan(1/2.),  4*PI/5,  2*PI/5,   -PI/5,   PI/5), //North Indian O.
				new Facet(-2.0,  1.0*sqrt(3),   -PI/6,  atan(1/2.),  4*PI/5,  4*PI/5,   -PI/5,   PI/5), //South Africa
				new Facet(-2.0,  0.0        ,   -PI/6, -atan(1/2.),   -PI  ,    PI/5,   -PI  ,   PI  ), //Australia
				new Facet(-1.0,  0.0        ,    PI/2, -atan(1/2.), -3*PI/5,    PI  , -6*PI/5, 4*PI/5), //China
				new Facet( 0.0,  0.0        ,    PI/2, -atan(1/2.),   -PI/5,    PI  ,   -PI  ,   PI  ), //North America
				new Facet( 1.0,  0.0        ,  5*PI/6, -atan(1/2.),    PI/5, -3*PI/5,   -PI  ,   PI  ), //East Pacific O.
				new Facet( 2.0,  0.0        ,    PI/2, -atan(1/2.),  3*PI/5,    PI  , -3*PI/5, 3*PI/5), //West Antarctica
				new Facet( 1.5, -0.5*sqrt(3),  5*PI/6, -atan(1/2.),  3*PI/5,   -PI/5,   -PI/5,   PI/5), //South Pacific O.
				new Facet(-2.5, -0.5*sqrt(3),    PI/6, -atan(1/2.),  3*PI/5,    PI/5,   -PI/5,   PI/5), //New Zealand
				new Facet( 0.5, -0.5*sqrt(3),    PI/2,       -PI/2,     0  ,   -PI  , -3*PI/5,   PI/5), //Hawai`i
				new Facet(-0.5, -0.5*sqrt(3),    PI/2,       -PI/2,     0  , -3*PI/5,   -PI/5, 2*PI/5), //West Pacific O.
				new Facet(-1.5, -0.5*sqrt(3),    PI/2,       -PI/2,     0  ,   -PI/5,     0  , 3*PI/5), //Melanesia
		});
		/** the number of the rotational symmetry in the spherical coordinate system */
		public final int sphereSym;
		/** the number of the rotational symmetry in the planar coordinate system */
		public final int planarSym;
		/** the bounding shape of the planar map */
		public final Shape shape;
		/** the location and orientation of each facet in both coordinate systems */
		public final Facet[] facets;
		/** the number of faces */
		public final Type type;

		Polyhedron(int sphereSym, int planarSym, Shape shape, Facet[] facets) {
			this.shape = shape;
			this.sphereSym = sphereSym;
			this.planarSym = planarSym;
			this.facets = facets;
			if (sphereSym == 3)
				this.type = Type.TETRAHEDRAL;
			else
				this.type = Type.ICOSAHEDRAL;
		}
		
		public double[] rotateOOB(double x, double y, double xCen, double yCen) { //move points that are out of bounds for project()
			return new double[] {x, y}; //this method should be overridden by projections with weird geometry
		}
		
		public boolean inBounds(double x, double y) {return true;} //determine whether a point is in bounds for inverse()
	}



	/**
	 * the information regarding a single facet of a polyhedron, including its location and orientation in 3D space and
	 * its location and orientation in the post-projection map plane
	 */
	private static class Facet {
		/** the x position of this Facet's main vertex in the plane */
		public final double x;
		/** the y position of this Facet's main vertex in the plane */
		public final double y;
		/** the orientation of this Facet with respect to its main vertex in the plane */
		public final double planeRotation;
		/** the latitude of this Facet's main vertex in 3D space */
		public final double latitude;
		/** the longitude of this Facet's main vertex in 3D space */
		public final double longitude;
		/** the orientation of this Facet with respect to its main vertex in 3D space */
		public final double centralMeridian;
		/** the orientation of the left bound of this Facet with respect to its main vertex in 3D space */
		public final double minAngle;
		/** the orientation of the right bound of this Facet with respect to its main vertex in 3D space */
		public final double maxAngle;

		public Facet(double x, double y, double planeRotation, double latitude, double longitude, double centralMeridian) {
			this(x, y, planeRotation, latitude, longitude, centralMeridian, -2*PI, 2*PI);
		}

		public Facet(double x, double y, double planeRotation, double latitude, double longitude, double centralMeridian, double minAngle, double maxAngle) {
			this.x = x;
			this.y = y;
			this.planeRotation = planeRotation;
			this.latitude = latitude;
			this.longitude = longitude;
			this.centralMeridian = centralMeridian;
			this.minAngle = minAngle;
			this.maxAngle = maxAngle;
		}
	}
}
