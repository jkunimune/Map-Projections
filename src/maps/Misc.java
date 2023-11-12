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

import de.jtem.ellipticFunctions.Jacobi;
import image.SVGMap.Command;
import maps.Projection.Property;
import maps.Projection.Type;
import org.apache.commons.math3.complex.Complex;
import utils.BoundingBox;
import utils.Elliptic;
import utils.NumericalAnalysis;

import java.util.List;

import static java.lang.Double.NaN;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
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
import static utils.Math2.max;
import static utils.Math2.sigone;

/**
 * All the projections that don't fit into any of the other categories.
 * 
 * @author jkunimune
 */
public class Misc {
	
	public static final Projection PEIRCE_QUINCUNCIAL =
			new Projection(
					"Peirce Quincuncial", "A conformal projection that uses complex elliptic functions.",
					new BoundingBox(2, 2), 0b1001, Type.OTHER, Property.CONFORMAL, 3) {
		
		private static final double K_RT_HALF = 1.854074677; //this is approx K(sqrt(1/2))
		
		public double[] project(double lat, double lon) {
			double quadNum = floor((lon-PI/4)/(PI/2));
			double wArg = lon - quadNum*PI/2;
			double wAbs = tan(PI/4-abs(lat)/2);
			Complex w = new Complex(wAbs*sin(wArg), -wAbs*cos(wArg)); //this Complex comes from Apache Commons
			Complex k = new Complex(sqrt(0.5));
			Complex z = Elliptic.F(w.acos(),k).divide(K_RT_HALF).subtract(Complex.ONE);
			z = z.multiply(Complex.I.pow(quadNum-2));
			double x = z.getReal(), y = z.getImaginary();
			if (lat < 0)
				return new double[] {sigone(x)*(1 - abs(y)),
						sigone(y)*(1 - abs(x))}; //reflect over equator if necessary
			else
				return new double[] {x, y};
		}
		
		public double[] inverse(double x, double y) {
			de.jtem.mfc.field.Complex u = new de.jtem.mfc.field.Complex(K_RT_HALF*(x+1), K_RT_HALF*y);
			de.jtem.mfc.field.Complex k = new de.jtem.mfc.field.Complex(sqrt(0.5)); //This Complex comes from a German university
			de.jtem.mfc.field.Complex ans = Jacobi.cn(u, k);
			double p = 2 * atan(ans.abs());
			double theta = atan2(-ans.getRe(), ans.getIm());
			double lambda = PI/2 - p;
			return new double[] {lambda, theta};
		}
	};
	
	
	public static final Projection GUYOU =
			new Projection(
					"Guyou", "Peirce Quincuncial, rearranged a bit.",
					new BoundingBox(2, 1), 0b1001,
					Type.OTHER, Property.CONFORMAL, 3) {
		
		private static final double K_RT_HALF = 1.854074677; //this is approx K(sqrt(1/2))
		private final double[] POLE = {0, -PI/2, PI/4};
		
		public double[] project(double lat, double lon) {
			final double[] coords = transformFromOblique(lat, lon, POLE);
			double quadNum = floor((coords[1]-PI/4)/(PI/2));
			double wArg = coords[1] - quadNum*PI/2;
			double wAbs = tan(PI/4-abs(coords[0])/2);
			Complex w = new Complex(wAbs*sin(wArg), -wAbs*cos(wArg)); //this Complex comes from Apache Commons
			Complex k = new Complex(sqrt(0.5));
			Complex z = Elliptic.F(w.acos(),k).divide(K_RT_HALF).subtract(Complex.ONE);
			z = z.multiply(Complex.I.pow(quadNum-1.5));
			double x = z.getReal()/sqrt(2), y = z.getImaginary()/sqrt(2);
			if (coords[0] < 0)
				return new double[] {.5 - x, y}; //reflect over equator if necessary
			else
				return new double[] {x - .5, y};
		}
		
		public double[] inverse(double x, double y) {
			de.jtem.mfc.field.Complex u = new de.jtem.mfc.field.Complex(x-y-.5, x+y+.5).times(K_RT_HALF);
			de.jtem.mfc.field.Complex k = new de.jtem.mfc.field.Complex(sqrt(0.5)); //just some fancy complex calculus stuff
			de.jtem.mfc.field.Complex ans = Jacobi.cn(u, k);
			double p = 2 * atan(ans.abs());
			double theta = ans.arg();
			double lambda = PI/2 - p;
			return transformToOblique(new double[] {lambda, theta}, POLE);
		}
	};
	
	
	public static final Projection HAMMER_RETROAZIMUTHAL =
			new Projection(
					"Hammer Retroazimuthal", "The full version of a map where bearing and distance to a reference point is preserved.",
					new BoundingBox(2*PI, 2*PI), 0b1110, Type.PSEUDOCONIC, Property.RETROAZIMUTHAL, 2,
					new String[] {"Latitude","Longitude"},
					new double[][] {{-89,89,21.4}, {-180,180,39.8}}, false) {
		
		private double phi0, lam0;
		
		public void initialize(double... params) {
			this.phi0 = toRadians(params[0]);
			this.lam0 = toRadians(params[1]);
		}
		
		public double[] project(double lat, double lon) {
			final double z = acos(sin(phi0)*sin(lat) +
					cos(phi0)*cos(lat)*cos(lon-lam0));
			final double K = z/sin(z);
			final double x = K*cos(phi0)*sin(lon-lam0);
			final double y = -K*(sin(phi0)*cos(lat) -
					cos(phi0)*sin(lat)*cos(lon-lam0));
			if (cos(lon-lam0) < 0)
				return new double[] {-x, -y};
			else
				return new double[] {x, y};
		}
		
		public double[] inverse(double x, double y) {
			double phi1 = PI/2 - hypot(x, y);
			if (phi1 < -PI/2) 	return null;
			double lam1 = atan2(x, -y);
			double phiP = asin(sin(phi0)/hypot(sin(phi1),cos(phi1)*cos(lam1))) - atan2(cos(phi1)*cos(lam1),sin(phi1));
			if (abs(phiP) > PI/2)
				phiP = signum(phiP)*PI - phiP;
			double delL = acos(sin(phi1)/cos(phiP)/cos(phi0) - tan(phiP)*tan(phi0));
			double lamP = lam0 + signum(x)*delL;
			if (isNaN(phiP) || isNaN(lamP)) 	return null;
			if (lamP > PI) 	lamP -= 2*PI;
			if (lamP < -PI) 	lamP += 2*PI;
			return new double[] {phiP, lamP};
		}
		
		@Override
		public double[] project(double lat, double lon, double[] pole) {
			return super.project(lat, lon, null);
		}
		
		@Override
		public double[] inverse(double x, double y, double[] pole, boolean crop) {
			return super.inverse(x, y, null, crop);
		}
	};
	
	
	public static final Projection TWO_POINT_EQUIDISTANT =
			new Projection(
					"Two-point Equidistant", "A map that preserves distances, but not azimuths, to two arbitrary points.",
					null, 0b1111, Type.OTHER, Property.EQUIDISTANT, 3,
					new String[] {"Latitude 1","Longitude 1","Latitude 2","Longitude 2"},
					new double[][] {{-90,90,41.9},{-180,180,12.5},{-90,90,34.7},{-180,180,112.4}},
					false) {
		
		private double lat1, lon1, lat2, lon2, D, a, c;
		
		public void initialize(double... params) {
			this.lat1 = toRadians(params[0]); //coordinates of first reference
			this.lon1 = toRadians(params[1]);
			this.lat2 = toRadians(params[2]); //coordinates of second reference
			this.lon2 = toRadians(params[3]);
			this.D = dist(lat1,lon1, lat2,lon2); //distance between references
			this.a = PI - D/2; //semimajor axis
			this.c = D/2; //focal distance
			double b = sqrt(pow(a, 2) - pow(c, 2)); //semiminor axis
			this.bounds = new BoundingBox(2*a, 2*b);
		}
		
		public double[] project(double lat0, double lon0) {
			final double d1 = dist(lat0,lon0, lat1,lon1);
			final double d2 = dist(lat0,lon0, lat2,lon2);
			final double s =signum(
					tan(lat0)*sin(lon2-lon1) +
					tan(lat1)*sin(lon0-lon2) +
					tan(lat2)*sin(lon1-lon0));
			return new double[] {
					(d1*d1-d2*d2)/(2*D),
					s*sqrt(d1*d1 - pow((d1*d1-d2*d2+D*D)/(2*D), 2)) };
		}
		
		public double[] inverse(double x, double y) {
			final double d1 = hypot(x + c, y);
			final double d2 = hypot(x - c, y);
			if (d1 + d2 > 2*a) 	return null; //TODO find out why it throws a hissy fit when y=0
			final double t1 = -(cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon1-lon2))/sin(D);
			double t2 = signum(lon1-lon2)*(cos(d1)*cos(D) - cos(d2))/(sin(d1)*sin(D));
			if (abs(t2) > 1) 	t2 = signum(t2);
			final double s0 = signum(lon1-lon2)*signum(y);
			final double casab = t1*t2 +s0* sqrt((t1*t1 - 1)*(t2*t2 - 1));
			final double s1 = signum(sin(acos(t1)-s0*acos(t2)));
			final double PHI = asin(sin(lat1)*cos(d1) - cos(lat1)*sin(d1)*casab);
			final double LAM = lon1 +s1* acos((cos(d1) - sin(lat1)*sin(PHI))/(cos(lat1)*cos(PHI)));
			return new double[] { PHI, coerceAngle(LAM) };
		}
		
		private double dist(double lat1, double lon1, double lat2, double lon2) {
			return acos(
					sin(lat1)*sin(lat2) +
					cos(lat1)*cos(lat2)*cos(lon1-lon2));
		}
		
		@Override
		public double[] project(double lat, double lon, double[] pole) {
			return super.project(lat, lon, null);
		}
		
		@Override
		public double[] inverse(double x, double y, double[] pole, boolean crop) {
			return super.inverse(x, y, null, crop);
		}
	};
	
	
	public static final Projection BRAUN_CONIC =
			new Projection(
					"Braun conic", "A particular perspective conic that is tangent at 30\u00B0.",
					new BoundingBox(-2*sqrt(3), 2*sqrt(3), -2*sqrt(3), 0), 0b1111,
					Type.CONIC, Property.PERSPECTIVE, 3) {
		
		public double[] project(double lat, double lon) {
			double r = 1.5*(sqrt(3) - tan((lat + PI/6)/2));
			double th = lon/2;
			double x = r*sin(th);
			double y = -r*cos(th);
			return new double[] {x, y};
		}
		
		public double[] inverse(double x, double y) {
			double r = hypot(x, y);
			double th = atan2(x, -y);
			if (r > bounds.xMax)
				return null;
			double lat = 2*atan(sqrt(3) - 2/3.*r) - PI/6;
			double lon = th*2;
			return new double[] {lat, lon};
		}
	};
	
	
	public static final Projection BONNE =
			new Projection(
					"Bonne", "A traditional pseudoconic projection, also known as the Sylvanus projection.",
					null, 0b1111, Type.PSEUDOCONIC, Property.EQUAL_AREA, 1,
					new String[] {"Std. Parallel"}, new double[][] {{-90, 90, 45}}) {
		
		private double r0;
		private boolean reversed; // if the standard parallel is southern
		
		public void initialize(double... params) {
			double lat0 = toRadians(params[0]);
			this.reversed = (lat0 < 0);
			if (reversed)
				lat0 = -lat0;
			this.r0 = 1/tan(lat0) + lat0;
			if (isInfinite(r0)) {
				this.bounds = Pseudocylindrical.SINUSOIDAL.bounds;
			}
			else { // for such a simple map projection...
				double argmaxX = NumericalAnalysis.bisectionFind(
						(p) -> (PI*(1/(r0-p)*cos(p) - sin(p))*cos(PI/(r0-p)*cos(p)) - sin(PI/(r0-p)*cos(p))),
						-PI/2, 0, 1e-3); // it sure is complicated to find its dimensions!
				double maxX = (r0 - argmaxX)*sin(PI/(r0 - argmaxX)*cos(argmaxX));
				double argmaxY;
				try {
					argmaxY = NumericalAnalysis.bisectionFind(
							(p) -> (PI*(1/(r0-p)*cos(p) - sin(p))*sin(PI/(r0-p)*cos(p)) + cos(PI/(r0-p)*cos(p))),
							0, PI/4, 1e-3);
				} catch (IllegalArgumentException e) {
					argmaxY = PI/2;
				}
				double maxY = max(-r0 + PI/2,
						-(r0 - argmaxY)*cos(PI/(r0 - argmaxY)*cos(argmaxY)));
				this.bounds = new BoundingBox(-maxX, maxX, -r0 - PI/2, maxY);
				if (reversed)
					this.bounds = new BoundingBox(-bounds.xMax, -bounds.xMin, -bounds.yMax, -bounds.yMin);
			}
		}
		
		public double[] project(double lat, double lon) {
			if (isInfinite(r0))
				return Pseudocylindrical.SINUSOIDAL.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			
			double r = r0 - lat;
			double th = lon*cos(lat)/r;
			double x = r*sin(th);
			double y = -r*cos(th);
			
			if (reversed)
				return new double[] {-x,-y};
			else
				return new double[] { x, y};
		}
		
		public double[] inverse(double x, double y) {
			if (isInfinite(r0))
				return Pseudocylindrical.SINUSOIDAL.inverse(x, y);
			if (reversed) {
				x = -x;
				y = -y;
			}
			
			double r = hypot(x, y);
			if (r < r0 - PI/2 || r > r0 + PI/2)
				return null;
			double th = atan2(x, -y);
			double lat = r0 - r;
			double lon = th*r/cos(lat);
			
			if (reversed)
				return new double[] {-lat,-lon};
			else
				return new double[] { lat, lon};
		}
	};
	
	
	public static final Projection T_SHIRT =
			new Projection(
					"T-Shirt", "A conformal projection onto a torso.",
					new BoundingBox(10, 6), 0b1001, Type.OTHER, Property.CONFORMAL, 3) {

		private final double[] X = {0, .507, .753, 1};
		private final double[] A = {.128, .084, .852, -.500};
		private final de.jtem.mfc.field.Complex K = new de.jtem.mfc.field.Complex(1, 0);
		
		public double[] project(double lat, double lon) {
			double wAbs = tan(PI/4-abs(lat)/2);
			de.jtem.mfc.field.Complex w = new de.jtem.mfc.field.Complex(wAbs*sin(lon), -wAbs*cos(lon));
			w = (w.plus(1)).divide(w.minus(1)).neg().timesI();
			de.jtem.mfc.field.Complex z = NumericalAnalysis.simpsonIntegrate(
					new de.jtem.mfc.field.Complex(0,1), w, this::integrand, 1e-2);
			double x = z.getIm(), y = -z.getRe();
			if (lat >= 0)
				return new double[] {x-2.5, y+1}; //move the back of the shirt over
			else
				return new double[] {2.5-x, y+1};
		}
		
		public double[] inverse(double x, double y) {
			return null;
		}
		
		private de.jtem.mfc.field.Complex integrand(de.jtem.mfc.field.Complex z) {
			de.jtem.mfc.field.Complex w = K;
			for (int i = X.length-1; i > 0; i --)
				w = w.times(z.plus(X[i]).pow(-A[i]));
			for (int i = 0; i < X.length; i ++)
				w = w.times(z.minus(X[i]).pow(-A[i]));
			return w;
		}
	};
	
	
	public static final Projection FLAT_EARTH =
			new Projection(
					"Flat Earth", "The one true map.", new BoundingBox(2, 2), 0b1111,
					Type.PLANAR, Property.TRUE, 5, new String[0], new double[0][], false) {
		
		private final double[] CORE_LONGITUDES = {
				toRadians(-60), toRadians(20), toRadians(135) };
		
		public double[] project(double lat, double lon) {
			double a;
			if (lat >= 0 || lat <= -PI/3)
				a = 1;
			else
				a = 1+(4-sqrt(6)+sqrt(2))/4*3/PI*lat;
			double lonR = POSITIVE_INFINITY; //the round earth model was computed using three central longitudes
			for (double lon0: CORE_LONGITUDES)
				if (abs(coerceAngle(lon-lon0)) < abs(lonR))
					lonR = coerceAngle(lon-lon0); //pick the one to which you are closest
			
			double r = tan((PI/2-lat)/4); //compute the distance from the center of the Universe
			double th = lonR*a + (lon-lonR); //adjust longitude to remove RET-induced distortion
			return new double[] {r*sin(th), -r*cos(th)};
		}
		
		public double[] inverse(double x, double y) {
			double lat = PI/2 - 4*atan(hypot(x, y)); //simulate the fake concept of latitude
			double a;
			if (lat >= 0 || lat <= -PI/3)
				a = 1;
			else
				a = 1+(4-sqrt(6)+sqrt(2))/4*3/PI*lat;
			
			if (lat < -PI/2) 	lat = -PI/2; //the ice wall extends as far as man knows
			double th = atan2(x, -y);
			
			double th0 = POSITIVE_INFINITY;
			for (double lon: CORE_LONGITUDES)
				if (isNaN(th0) ||
						abs(coerceAngle(th-lon)) < abs(th-th0))
					th0 = lon; //find the central longitude to which you are closest
			double thR = coerceAngle(th - th0);
			
			double th1 = NaN; //to fix the empty space
			for (double lon: CORE_LONGITUDES) {
				if (isNaN(th1) ||
						floorMod((lon-th)*signum(thR), 2*PI) <
						floorMod((th1-th)*signum(thR), 2*PI))
					th1 = lon; //find the central longitude on your other side
			}
			double thM = coerceAngle((th0 + th1)/2); //the line between the two
			if (abs(thR)/a > abs(coerceAngle(thM-th0))) { //if you are farther than the midpoint is, once adjusted
				if (thR>0 == coerceAngle(thM-th0)>0)
					return new double[] {lat, thM}; //then this point does not exist in the sphere earth model; fill it with guess
				else
					return new double[] {lat, coerceAngle(thM+PI)};
			}
			
			return new double[] {lat, coerceAngle(thR/a + th0)};
		}
		
		@Override
		public double[] getDistortionAt(double[] s0) {
			return new double[] {0,0};
		}
		
		@Override
		public List<Command> drawGraticule(double spacing, double precision, double outW, double outH,
		                                          double maxLat, double maxLon, double[] pole) {
			return Azimuthal.POLAR.drawGraticule(spacing, precision, outW, outH, maxLat, maxLon, null);
		}
	};
	
}
