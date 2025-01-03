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
import image.Path;
import maps.Projection.Property;
import maps.Projection.Type;
import org.apache.commons.math3.complex.Complex;
import utils.Elliptic;
import utils.NumericalAnalysis;
import utils.Shape;

import java.util.ArrayList;
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
import static utils.Math2.sigone;

/**
 * All the projections that don't fit into any of the other categories.
 * 
 * @author jkunimune
 */
public class Misc {
	
	public static final Projection PEIRCE_QUINCUNCIAL =
			new Projection(
					"Peirce Quincuncial", "Charles S. Pierce", "A conformal projection that uses complex elliptic functions",
					Shape.rectangle(2, 2), false, true, false, false, Type.OTHER, Property.CONFORMAL, 3) {

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
					"Guyou", "Emile Guyou", "Peirce Quincuncial, rearranged a bit",
					Shape.rectangle(2, 1), false, true, false, false,
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
	
	
	public static final Projection TWO_POINT_EQUIDISTANT =
			new Projection(
					"Two-point Equidistant", "Hans Maurer", "A map that preserves distances, but not azimuths, to two arbitrary points",
					null, true, true, true, true, Type.OTHER, Property.EQUIDISTANT, 3,
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
			this.shape = Shape.ellipse(a, b);
		}
		
		public double[] project(double lat0, double lon0) {
			final double d1 = dist(lat0,lon0, lat1,lon1);
			final double d2 = dist(lat0,lon0, lat2,lon2);
			final double s = (
					tan(lat0)*sin(lon2-lon1) +
					tan(lat1)*sin(lon0-lon2) +
					tan(lat2)*sin(lon1-lon0) > 0) ? 1 : -1;
			return new double[] {
					(d1*d1-d2*d2)/(2*D),
					s*sqrt(d1*d1 - pow((d1*d1-d2*d2+D*D)/(2*D), 2)) };
		}
		
		public double[] inverse(double x, double y) {
			final double d1 = hypot(x + c, y);
			final double d2 = hypot(x - c, y);
			if (d1 + d2 > 2*a) 	return null; //TODO find out why it throws a hissy fit when y=0
			final double t1 = -(cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon1-lon2))/sin(D);
			double t2 = (lon1 > lon2 ? 1 : -1)*(cos(d1)*cos(D) - cos(d2))/(sin(d1)*sin(D));
			t2 = Math.max(-1., Math.min(1., t2));
			final double s0 = ((lon1 > lon2) == (y > 0)) ? 1 : -1;
			final double casab = t1*t2 +s0* sqrt((t1*t1 - 1)*(t2*t2 - 1));
			final double s1 = coerceAngle(acos(t1)-s0*acos(t2)) > 0 ? 1 : -1;
			final double PHI = asin(sin(lat1)*cos(d1) - cos(lat1)*sin(d1)*casab);
			final double LAM = lon1 + s1*acos((cos(d1) - sin(lat1)*sin(PHI))/(cos(lat1)*cos(PHI)));
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
					"Braun conic", "Carl Braun", "A particular perspective conic that is tangent at 30\u00B0",
					Shape.annularSector(0, 2*sqrt(3), PI, false), true, true, true, true,
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
			if (r > shape.xMax)
				return null;
			double lat = 2*atan(sqrt(3) - 2/3.*r) - PI/6;
			double lon = th*2;
			return new double[] {lat, lon};
		}
	};
	
	
	public static final Projection BONNE =
			new Projection(
					"Bonne", "Bernardo Sylvano", "A traditional pseudoconic projection, also known as the Sylvanus projection",
					null, true, true, true, true, Type.PSEUDOCONIC, Property.EQUAL_AREA, 1,
					new String[] {"Std. Parallel"}, new double[][] {{-90, 90, 45}}) {
		
		private double r0;
		private boolean reversed; // if the standard parallel is southern
		
		public void initialize(double... params) {
			double lat0 = toRadians(params[0]);
			this.reversed = (lat0 < 0);
			if (reversed)
				lat0 = -lat0;
			this.r0 = 1/tan(lat0) + lat0;
			this.shape = Shape.meridianEnvelope(this);
		}
		
		public double[] project(double lat, double lon) {
			if (isInfinite(r0))
				return Pseudocylindrical.SINUSOIDAL.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			
			double r = r0 - lat;
			double th, x, y;
			if (r > 0) {
				th = lon*cos(lat)/r;
				x = r*sin(th);
				y = -r*cos(th);
			}
			else {
				x = 0;
				y = 0;
			}
			
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
					"T-Shirt", "Justin H. Kunimune", "A conformal projection onto a torso",
					Shape.polygon(new double[][] {
							{ 0.000, 1.784},
							{-1.17, 2.38},
							{-2.500, 2.651},
							{-3.83, 2.38},
							{-5.000, 1.784},
							{-5.000, 0.933},
							{-4.071, 0.933},
							{-4.071, -2.500},
							{-0.929, -2.500},
							{-0.929, 0.933},
							{ 0.000, 0.933},
							{ 0.929, 0.933},
							{ 0.929, -2.500},
							{ 4.071, -2.500},
							{ 4.071, 0.933},
							{ 5.000, 0.933},
							{ 5.000, 1.784},
							{ 3.83, 2.38},
							{ 2.500, 2.651},
							{ 1.17, 2.38},
					}),
					false, true, true, false, Type.OTHER, Property.CONFORMAL, 3) {

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
	
	
	public static final Projection CASSINI = new Projection(
			"Cassini", "Cesar-François Cassini de Thury", "A transverse Plate–Carée projection",
			Shape.rectangle(PI, 2*PI), true, true, true, true, Type.CYLINDRICAL, Property.EQUIDISTANT, 2) {
		
		public double[] project(double lat, double lon) {
			double x = asin(cos(lat)*sin(lon));  // I could use obliquifySph() and EQUIRECTANGULAR for this
			double y = atan2(tan(lat), cos(lon));  // but since there are special simple equations I may as well use them
			return new double[] {x, y};  // also this is technically rotated 90° in the plane
		}
		
		public double[] inverse(double x, double y) {
			double lat = asin(sin(y)*cos(x));
			double lon = atan2(tan(x), cos(y));
			return new double[] {lat, lon};
		}
	};
	
	
	public static final Projection LEMONS = new Projection(
			"Gores", "(unknown)", "A heavily interrupted projection that can be pasted onto a globe",
			null, false, true, true, true, Type.OTHER, Property.COMPROMISE, 2,
			new String[] {"Number of gores"}, new double[][] {{4, 72, 12}}) {
		
		private int numLemons;
		private double lemonWidth;
		
		public void initialize(double... params) {
			// read in the number of gores and calculate the longitudinal extent of each
			this.numLemons = (int) Math.round(params[0]);
			this.lemonWidth = 2*PI/numLemons;
			
			// to calculate the shape, start by getting the shape of a generic meridian
			List<Path.Command> polewardSegment = CASSINI.drawLoxodrome(
					0, lemonWidth/2, PI/2, lemonWidth/2, .1);
			// make an equator-to-pole version and a pole-to-equator version
			List<Path.Command> tropicwardSegment = Path.reversed(polewardSegment);
			// remove one endpoint from each so there are no duplicate vertices
			polewardSegment = polewardSegment.subList(0, polewardSegment.size() - 1);
			tropicwardSegment = tropicwardSegment.subList(0, tropicwardSegment.size() - 1);
			// then build up the full shape by transforming the generic segments
			List<Path.Command> envelope = new ArrayList<>(numLemons*4*polewardSegment.size());
			// go east to west in the north hemisphere
			for (int i = numLemons - 1; i >= 0; i --) {
				envelope.addAll(Path.transformed(1, 1, (i - (numLemons - 1)/2.)*lemonWidth, 0,
				                                 polewardSegment));
				envelope.addAll(Path.transformed(-1, 1, (i - (numLemons - 1)/2.)*lemonWidth, 0,
				                                 tropicwardSegment));
			}
			// go west to east in the south hemisphere
			for (int i = 0; i < numLemons; i ++) {
				envelope.addAll(Path.transformed(-1, -1, (i - (numLemons - 1)/2.)*lemonWidth, 0,
				                                 polewardSegment));
				envelope.addAll(Path.transformed(1, -1, (i - (numLemons - 1)/2.)*lemonWidth, 0,
				                                 tropicwardSegment));
			}
			// finally, convert it all to a Shape
			this.shape = Shape.polygon(Path.asArray(envelope));
		}
		
		public double[] project(double lat, double lon) {
			double lemonIndex = floor((lon + PI)/lemonWidth);  // pick a lemon
			if (lon == 2*PI)
				lemonIndex = numLemons - 1;
			final double lemonCenter = (lemonIndex - numLemons/2. + 1/2.)*lemonWidth;
			final double dl = lon - lemonCenter;  // find the relative longitude
			double[] xy = CASSINI.project(lat, dl);  // project on Cassini with that
			xy[0] += lemonCenter;  // shift according to the lemon's x center
			return xy;
		}
		
		public double[] inverse(double x, double y) {
			final int lemonIndex = (int)floor((x + PI)/lemonWidth);  // pick a lemon
			final double lemonCenter = (lemonIndex - numLemons/2. + 1/2.)*lemonWidth;
			final double dx = x - lemonCenter;  // find the relative x
			double[] latLon = CASSINI.inverse(dx, y);  // project from Cassini with that
			if (abs(latLon[1]) > lemonWidth/2)
				latLon[1] += 4*PI;  // mark it as out of bounds if it's out of the bounds of its lemon
			latLon[1] += lemonCenter;  // set the absolute longitude based on which lemon it is
			return latLon;
		}
	};
	
	
	public static final Projection FLAT_EARTH =
			new Projection(
					"Flat Earth", "Samuel B. Rowbotham", "The one true map", Shape.circle(1), true, true, true, true,
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
		public List<Path.Command> drawGraticule(
				double spacing, boolean sparsePole, double precision, double outW, double outH,
				double maxLat, double maxLon, double[] pole) {
			return Azimuthal.POLAR.drawGraticule(spacing, sparsePole, precision, outW, outH, maxLat, maxLon, null);
		}
	};
	
}
