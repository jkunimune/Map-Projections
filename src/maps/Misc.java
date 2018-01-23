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

import org.apache.commons.math3.complex.Complex;

import de.jtem.ellipticFunctions.Jacobi;
import maps.Projection.Property;
import maps.Projection.Type;
import utils.Elliptic;
import utils.Math2;
import utils.SVGMap.Path;

/**
 * All the projections that don't fit into any of the other categories.
 * 
 * @author jkunimune
 */
public class Misc {
	
	public static final Projection PEIRCE_QUINCUNCIAL =
			new Projection(
					"Peirce Quincuncial", "A conformal projection that uses complex elliptic functions.",
					2, 2, 0b1001, Type.OTHER, Property.CONFORMAL, 3) {
		
		private static final double K_RT_HALF = 1.85407; //this is approx K(sqrt(1/2))
		
		public double[] project(double lat, double lon) {
			final double alat = Math.abs(lat);
			final double wMag = Math.tan(Math.PI/4-alat/2);
			final Complex w = new Complex(wMag*Math.sin(lon), -wMag*Math.cos(lon)); //this Complex comes from Apache Commons
			final Complex k = new Complex(Math.sqrt(0.5));
			Complex z = Elliptic.F(w.acos(),k).divide(K_RT_HALF).subtract(Complex.ONE).negate();
			if (z.isInfinite() || z.isNaN())	z = new Complex(0);
			double x = z.getReal(), y = z.getImaginary();
			
			if (lat < 0) {
				if (x >= 0 && y >= 0)
					z = new Complex(1-y, 1-x);
				else if (x >= 0 && y < 0)
					z = new Complex(1+y, -1+x);
				else if (y >= 0)
					z = new Complex(-1+y, 1+x);
				else
					z = new Complex(-1-y, -1-x);
			}
			return new double[] {z.getReal(), z.getImaginary()};
		}
		
		public double[] inverse(double x, double y) {
			de.jtem.mfc.field.Complex u = new de.jtem.mfc.field.Complex(K_RT_HALF*(x+1), K_RT_HALF*y);
			de.jtem.mfc.field.Complex k = new de.jtem.mfc.field.Complex(Math.sqrt(0.5)); //the rest comes from some fancy complex calculus
			de.jtem.mfc.field.Complex ans = Jacobi.cn(u, k);
			double p = 2 * Math.atan(ans.abs());
			double theta = Math.atan2(-ans.getRe(), ans.getIm());
			double lambda = Math.PI/2 - p;
			return new double[] {lambda, theta};
		}
	};
	
	
	public static final Projection GUYOU =
			new Projection(
					"Guyou", "Peirce Quincuncial, rearranged a bit.", 2., 1, 0b1001,
					Type.OTHER, Property.CONFORMAL, 3) {
		
		private static final double K_RT_HALF = 1.854; //this is approx K(sqrt(1/2))
		private final double[] POLE = {0, -Math.PI/2, Math.PI/4};
		
		public double[] project(double lat, double lon) {
			final double[] coords = obliquifySphc(lat,lon, POLE);
			final double alat = Math.abs(coords[0]);
			final double wMag = Math.tan(Math.PI/4-alat/2);
			final Complex w = new Complex(wMag*Math.sin(coords[1]), -wMag*Math.cos(coords[1]));
			final Complex k = new Complex(Math.sqrt(0.5));
			Complex z = Elliptic.F(w.acos(),k).divide(-2*K_RT_HALF)
					.multiply(new Complex(1,1)).add(new Complex(0,0.5));
			if (z.isInfinite() || z.isNaN()) 	z = new Complex(0);
			if (coords[0] < 0) 	z = z.conjugate().negate();
			return new double[] {z.getReal(), z.getImaginary()};
		}
		
		public double[] inverse(double x, double y) {
			de.jtem.mfc.field.Complex u = new de.jtem.mfc.field.Complex(x-y-.5, x+y+.5).times(K_RT_HALF);
			de.jtem.mfc.field.Complex k = new de.jtem.mfc.field.Complex(Math.sqrt(0.5)); //just some fancy complex calculus stuff
			de.jtem.mfc.field.Complex ans = Jacobi.cn(u, k);
			double p = 2 * Math.atan(ans.abs());
			double theta = ans.arg();
			double lambda = Math.PI/2 - p;
			return obliquifyPlnr(new double[] {lambda,theta}, POLE);
		}
	};
	
	
	public static final Projection HAMMER_RETROAZIMUTHAL =
			new Projection(
					"Hammer Retroazimuthal", "The full version of a map where bearing and distance to a reference point is preserved.",
					2*Math.PI, 2*Math.PI, 0b1110, Type.QUASIAZIMUTHAL, Property.RETROAZIMUTHAL, 2,
					new String[] {"Latitude","Longitude"},
					new double[][] {{-89,89,21.4}, {-180,180,39.8}}, false) {
		
		private double phi0, lam0;
		
		public void setParameters(double... params) {
			this.phi0 = Math.toRadians(params[0]);
			this.lam0 = Math.toRadians(params[1]);
		}
		
		public double[] project(double lat, double lon) {
			final double z = Math.acos(Math.sin(phi0)*Math.sin(lat) +
					Math.cos(phi0)*Math.cos(lat)*Math.cos(lon-lam0));
			final double K = z/Math.sin(z);
			final double x = K*Math.cos(phi0)*Math.sin(lon-lam0);
			final double y = -K*(Math.sin(phi0)*Math.cos(lat) -
					Math.cos(phi0)*Math.sin(lat)*Math.cos(lon-lam0));
			if (Math.cos(lon-lam0) < 0)
				return new double[] {-x, -y};
			else
				return new double[] {x, y};
		}
		
		public double[] inverse(double x, double y) {
			double phi1 = Math.PI/2 - Math.hypot(x, y);
			if (phi1 < -Math.PI/2) 	return null;
			double lam1 = Math.atan2(x, -y);
			double phiP = Math.asin(Math.sin(phi0)/Math.hypot(Math.sin(phi1),Math.cos(phi1)*Math.cos(lam1))) - Math.atan2(Math.cos(phi1)*Math.cos(lam1),Math.sin(phi1));
			if (Math.abs(phiP) > Math.PI/2)
				phiP = Math.signum(phiP)*Math.PI - phiP;
			double delL = Math.acos(Math.sin(phi1)/Math.cos(phiP)/Math.cos(phi0) - Math.tan(phiP)*Math.tan(phi0));
			double lamP = lam0 + Math.signum(x)*delL;
			if (Double.isNaN(phiP) || Double.isNaN(lamP)) 	return null;
			if (lamP > Math.PI) 	lamP -= 2*Math.PI;
			if (lamP < -Math.PI) 	lamP += 2*Math.PI;
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
					0, 0, 0b1111, Type.QUASIAZIMUTHAL, Property.EQUIDISTANT, 3,
					new String[] {"Latitude 1","Longitude 1","Latitude 2","Longitude 2"},
					new double[][] {{-90,90,41.9},{-180,180,12.5},{-90,90,34.7},{-180,180,112.4}},
					false) {
		
		private double lat1, lon1, lat2, lon2, D, a, b, c;
		
		public void setParameters(double... params) {
			this.lat1 = Math.toRadians(params[0]); //coordinates of first reference
			this.lon1 = Math.toRadians(params[1]);
			this.lat2 = Math.toRadians(params[2]); //coordinates of second reference
			this.lon2 = Math.toRadians(params[3]);
			this.D = dist(lat1,lon1, lat2,lon2); //distance between references
			this.a = Math.PI - D/2; //semimajor axis
			this.c = D/2; //focal distance
			this.b = Math.sqrt(Math.pow(a, 2) - Math.pow(c, 2)); //semiminor axis
			this.width = 2*a;
			this.height = 2*b;
		}
		
		public double[] project(double lat0, double lon0) {
			final double d1 = dist(lat0,lon0, lat1,lon1);
			final double d2 = dist(lat0,lon0, lat2,lon2);
			final double s =Math.signum(
					Math.tan(lat0)*Math.sin(lon2-lon1) +
					Math.tan(lat1)*Math.sin(lon0-lon2) +
					Math.tan(lat2)*Math.sin(lon1-lon0));
			return new double[] {
					(d1*d1-d2*d2)/(2*D),
					s*Math.sqrt(d1*d1 - Math.pow((d1*d1-d2*d2+D*D)/(2*D), 2)) };
		}
		
		public double[] inverse(double x, double y) {
			final double d1 = Math.hypot(x + c, y);
			final double d2 = Math.hypot(x - c, y);
			if (d1 + d2 > 2*a) 	return null; //TODO find out why it throws a hissy fit when y=0
			final double t1 = -(Math.cos(lat1)*Math.sin(lat2) - Math.sin(lat1)*Math.cos(lat2)*Math.cos(lon1-lon2))/Math.sin(D);
			double t2 = Math.signum(lon1-lon2)*(Math.cos(d1)*Math.cos(D) - Math.cos(d2))/(Math.sin(d1)*Math.sin(D));
			if (Math.abs(t2) > 1) 	t2 = Math.signum(t2);
			final double s0 = Math.signum(lon1-lon2)*Math.signum(y);
			final double casab = t1*t2 +s0* Math.sqrt((t1*t1 - 1)*(t2*t2 - 1));
			final double s1 = Math.signum(Math.sin(Math.acos(t1)-s0*Math.acos(t2)));
			final double PHI = Math.asin(Math.sin(lat1)*Math.cos(d1) - Math.cos(lat1)*Math.sin(d1)*casab);
			final double LAM = lon1 +s1* Math.acos((Math.cos(d1) - Math.sin(lat1)*Math.sin(PHI))/(Math.cos(lat1)*Math.cos(PHI)));
			return new double[] { PHI, Math2.coerceAngle(LAM) };
		}
		
		private double dist(double lat1, double lon1, double lat2, double lon2) {
			return Math.acos(
					Math.sin(lat1)*Math.sin(lat2) +
					Math.cos(lat1)*Math.cos(lat2)*Math.cos(lon1-lon2));
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
	
	
	public static final Projection FLAT_EARTH =
			new Projection(
					"Flat Earth", "The one true map.", 2, 2, 0b1111,
					Type.PLANAR, Property.TRUE, 5, new String[0], new double[0][], false) {
		
		private final double[] CORE_LONGITUDES = {
				Math.toRadians(-60), Math.toRadians(20), Math.toRadians(135) };
		
		public double[] project(double lat, double lon) {
			double a;
			if (lat >= 0 || lat <= -Math.PI/3)
				a = 1;
			else
				a = 1+(4-Math.sqrt(6)+Math.sqrt(2))/4*3/Math.PI*lat;
			double lonR = Double.POSITIVE_INFINITY; //the round earth model was computed using three central longitudes
			for (double lon0: CORE_LONGITUDES)
				if (Math.abs(Math2.coerceAngle(lon-lon0)) < Math.abs(lonR))
					lonR = Math2.coerceAngle(lon-lon0); //pick the one to which you are closest
			
			double r = Math.tan((Math.PI/2-lat)/4); //compute the distance from the center of the Universe
			double th = lonR*a + (lon-lonR); //adjust longitude to remove RET-induced distortion
			return new double[] {r*Math.sin(th), -r*Math.cos(th)};
		}
		
		public double[] inverse(double x, double y) {
			double lat = Math.PI/2 - 4*Math.atan(Math.hypot(x, y)); //simulate the fake concept of latitude
			double a;
			if (lat >= 0 || lat <= -Math.PI/3)
				a = 1;
			else
				a = 1+(4-Math.sqrt(6)+Math.sqrt(2))/4*3/Math.PI*lat;
			
			if (lat < -Math.PI/2) 	lat = -Math.PI/2; //the ice wall extends as far as man knows
			double th = Math.atan2(x, -y);
			
			double th0 = Double.POSITIVE_INFINITY;
			for (double lon: CORE_LONGITUDES)
				if (Double.isNaN(th0) ||
						Math.abs(Math2.coerceAngle(th-lon)) < Math.abs(th-th0))
					th0 = lon; //find the central longitude to which you are closest
			double thR = Math2.coerceAngle(th - th0);
			
			double th1 = Double.NaN; //to fix the empty space
			for (double lon: CORE_LONGITUDES) {
				if (Double.isNaN(th1) ||
						Math2.floorMod((lon-th)*Math.signum(thR), 2*Math.PI) <
						Math2.floorMod((th1-th)*Math.signum(thR), 2*Math.PI))
					th1 = lon; //find the central longitude on your other side
			}
			double thM = Math2.coerceAngle((th0 + th1)/2); //the line between the two
			if (Math.abs(thR)/a > Math.abs(Math2.coerceAngle(thM-th0))) { //if you are farther than the midpoint is, once adjusted
				if (thR>0 == Math2.coerceAngle(thM-th0)>0)
					return new double[] {lat, thM}; //then this point does not exist in the sphere earth model; fill it with guess
				else
					return new double[] {lat, Math2.coerceAngle(thM+Math.PI)};
			}
			
			return new double[] {lat, Math2.coerceAngle(thR/a + th0)};
		}
		
		@Override
		public double[] getDistortionAt(double[] s0) {
			return new double[] {0,0};
		}
		
		@Override
		public Path drawGraticule(double spacing, double precision, double outW, double outH,
			double maxLat, double maxLon, double[] pole) {
			return Azimuthal.POLAR.drawGraticule(spacing, precision, outW, outH, maxLat, maxLon, pole);
		}
	};
	
}
