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

import ellipticFunctions.Jacobi;
import maps.Projection.Property;
import maps.Projection.Type;
import utils.Elliptic;

/**
 * All the projections that don't fit into any of the other categories.
 * 
 * @author jkunimune
 */
public class Misc {
	
	public static final Projection AITOFF =
			new Projection("Aitoff", "A compromise projection shaped like an ellipse",
					2., 0b1011, Type.PSEUDOAZIMUTHAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			final double a = Math.acos(Math.cos(lat)*Math.cos(lon/2));
			return new double[] {
					2*Math.cos(lat)*Math.sin(lon/2)*a/Math.sin(a),
					Math.sin(lat)*a/Math.sin(a)};
		}
		
		public double[] inverse(double x, double y) {
			final double[] intermediate = Azimuthal.POLAR.inverse(x/2, y/2);
			double[] transverse = obliquifyPlnr(intermediate, new double[] {0,0,0});
			if (transverse != null) 	transverse[1] *= 2;
			return transverse;
		}
	};
	
	
	public static final Projection HAMMER =
			new Projection("Hammer", "An equal-area projection shaped like an ellipse",
					2., 0b1111, Type.PSEUDOAZIMUTHAL, Property.EQUAL_AREA) {
		
		public double[] project(double lat, double lon) {
			return new double[] {
					Math.PI*Math.cos(lat)*Math.sin(lon/2)/Math.sqrt(1+Math.cos(lat)*Math.cos(lon/2)),
					Math.PI/2*Math.sin(lat)/Math.sqrt(1+Math.cos(lat)*Math.cos(lon/2)) };
		}
		
		public double[] inverse(double x, double y) {
			final double X = x * Math.sqrt(8);
			final double Y = y * Math.sqrt(2);
			final double z = Math.sqrt(1 - Math.pow(X/4, 2) - Math.pow(Y/2, 2));
			final double shift = (Math.hypot(x, y) > 1) ? 2*Math.PI*Math.signum(x) : 0;
			return new double[] {
					Math.asin(z * Y), 2*Math.atan(0.5*z*X / (2*z*z - 1)) + shift};
		}
	};
	
	
	public static final Projection VAN_DER_GRINTEN =
			new Projection("Van der Grinten", "A circular compromise map that is popular for some reason",
					1., 0b1111, Type.OTHER, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			if (lat == 0) //special case 1: equator
				return new double[] {lon, 0};
			if (lon == 0 || lat >= Math.PI/2 || lat <= -Math.PI/2) //special case 3: prime meridian
				return new double[] {0, Math.PI*Math.tan(Math.asin(2*lat/Math.PI)/2)};
			
			final double t = Math.abs(Math.asin(2*lat/Math.PI));
			final double A = Math.abs(Math.PI/lon - lon/Math.PI)/2;
			final double G = Math.cos(t)/(Math.sin(t)+Math.cos(t)-1);
			final double P = G*(2/Math.sin(t) - 1);
			final double Q = A*A + G;
			return new double[] {
					Math.PI*Math.signum(lon)*(A*(G-P*P)+Math.sqrt(A*A*(G-P*P)*(G-P*P)-(P*P+A*A)*(G*G-P*P)))/(P*P+A*A),
					Math.PI*Math.signum(lat)*(P*Q-A*Math.sqrt((A*A+1)*(P*P+A*A)-Q*Q))/(P*P+A*A)};
		}
		
		public double[] inverse(double x, double y) {
			if (y == 0) // special case 1: equator
				return new double[] {0, x*Math.PI};
			if (x == 0) // special case 3: prime meridian
				return new double[] {Math.PI/2 * Math.sin(2*Math.atan(y)), 0};
			
			double c1 = -Math.abs(y) * (1 + x*x + y*y);
			double c2 = c1 - 2*y*y + x*x;
			double c3 = -2 * c1 + 1 + 2*y*y + Math.pow(x*x + y*y, 2);
			double d = y*y / c3 + 1 / 27.0 * (2*Math.pow(c2 / c3, 3) - 9*c1*c2 / (c3*c3));
			double a1 = 1 / c3*(c1 - c2*c2 / (3*c3));
			double m1 = 2 * Math.sqrt(-a1 / 3);
			double t1 = Math.acos(3*d / (a1 * m1)) / 3;
			return new double[] {
					Math.signum(y) * Math.PI * (-m1 * Math.cos(t1 + Math.PI/3) - c2 / (3*c3)),
					Math.PI*(x*x + y*y - 1 + Math.sqrt(1 + 2*(x*x - y*y) + Math.pow(x*x + y*y, 2)))
							/ (2*x)};
		}
	};
	
	
	public static final Projection PEIRCE_QUINCUNCIAL =
			new Projection("Peirce Quincuncial", "A conformal projection that uses complex elliptic functions",
					1., 0b1001, Type.OTHER, Property.CONFORMAL) {
		
		public double[] project(double lat, double lon) {
			final double alat = Math.abs(lat);
			final double wMag = Math.tan(Math.PI/4-alat/2);
			final Complex w = new Complex(wMag*Math.sin(lon), -wMag*Math.cos(lon));
			final Complex k = new Complex(Math.sqrt(0.5));
			Complex z = Elliptic.F(w.acos(),k).multiply(Math.PI/1.854).subtract(Math.PI).negate();
			if (z.isInfinite() || z.isNaN())	z = new Complex(0);
			double x = z.getReal(), y = z.getImaginary();
			
			if (lat < 0) {
				if (x >= 0 && y >= 0)
					z = new Complex(Math.PI-y, Math.PI-x);
				else if (x >= 0 && y < 0)
					z = new Complex(Math.PI+y, -Math.PI+x);
				else if (y >= 0)
					z = new Complex(-Math.PI+y, Math.PI+x);
				else
					z = new Complex(-Math.PI-y, -Math.PI-x);
			}
			return new double[] {z.getReal(), z.getImaginary()};
		}
		
		public double[] inverse(double x, double y) {
			mfc.field.Complex u = new mfc.field.Complex(1.854*(x+1), 1.854*y); // 1.854 is approx K(sqrt(1/2)
			mfc.field.Complex k = new mfc.field.Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus
			mfc.field.Complex ans = Jacobi.cn(u, k);
			double p = 2 * Math.atan(ans.abs());
			double theta = ans.arg() - Math.PI/2;
			double lambda = Math.PI/2 - p;
			return new double[] {lambda, theta};
		}
	};
	
	
	public static final Projection GUYOU =
			new Projection("Guyou", "Peirce Quincuncial, rearranged a bit", 2., 0b1001,
					Type.OTHER, Property.CONFORMAL) {
		
		private final double[] POLE = {0, -Math.PI/2, Math.PI/4};
		
		public double[] project(double lat, double lon) {
			final double[] coords = obliquifySphc(lat,lon, POLE);
			final double alat = Math.abs(coords[0]);
			final double wMag = Math.tan(Math.PI/4-alat/2);
			final Complex w = new Complex(wMag*Math.sin(coords[1]), -wMag*Math.cos(coords[1]));
			final Complex k = new Complex(Math.sqrt(0.5));
			Complex z = Elliptic.F(w.acos(),k).multiply(new Complex(Math.PI/3.708,Math.PI/3.708)).subtract(new Complex(0,Math.PI/2));
			if (z.isInfinite() || z.isNaN()) 	z = new Complex(0);
			if (coords[0] < 0) 	z = z.conjugate().negate();
			return new double[] {z.getReal(), z.getImaginary()};
		}
		
		public double[] inverse(double x, double y) {
			mfc.field.Complex u = new mfc.field.Complex(1.8558*(x - y/2 - 0.5), 1.8558*(x + y/2 + 0.5)); // don't ask me where 3.7116 comes from
			mfc.field.Complex k = new mfc.field.Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus
			mfc.field.Complex ans = Jacobi.cn(u, k);
			double p = 2 * Math.atan(ans.abs());
			double theta = ans.arg();
			double lambda = Math.PI/2 - p;
			return obliquifyPlnr(new double[] {lambda,theta}, POLE);
		}
	};
	
	
	public static final Projection HAMMER_RETROAZIMUTHAL =
			new Projection("Hammer Retroazimuthal", "The full version of a map where bearing and distance to a reference point is preserved",
					1., 0b1110, Type.QUASIAZIMUTHAL, Property.RETROAZIMUTHAL,
					new String[] {"Latitude","Longitude"},
					new double[][] {{-89,89,21.4},{-180,180,39.8}}, false) {
		
		private double phi0, lam0;
		
		public void setParameters(double... params) {
			this.phi0 = Math.toRadians(params[0]);
			this.lam0 = Math.toRadians(params[1]);
		}
		
		public double[] project(double lat, double lon, double[] pole) {
			return project(lat, lon);
		}
		
		public double[] project(double lat, double lon) {
			final double z = Math.acos(Math.sin(phi0)*Math.sin(lat) + Math.cos(phi0)*Math.cos(lat)*Math.cos(lon-lam0));
			final double K = z/Math.sin(z);
			final double x = K*Math.cos(phi0)*Math.sin(lon-lam0);
			final double y = -K*(Math.sin(phi0)*Math.cos(lat) - Math.cos(phi0)*Math.sin(lat)*Math.cos(lon-lam0));
			if (Math.cos(lon-lam0) < 0)
				return new double[] {-x, -y};
			else
				return new double[] {x, y};
		}
		
		public double[] inverse(double x, double y, double[] pole) {
			return inverse(x, y);
		}
		
		public double[] inverse(double x, double y) {
			final double phi1 = Math.PI/2 - Math.hypot(x, y)*Math.PI;
			if (phi1 < -Math.PI/2) 	return null;
			final double lam1 = Math.atan2(x, -y);
			double phiP = Math.asin(Math.sin(phi0)/Math.hypot(Math.sin(phi1),Math.cos(phi1)*Math.cos(lam1))) - Math.atan2(Math.cos(phi1)*Math.cos(lam1),Math.sin(phi1));
			if (Math.abs(phiP) > Math.PI/2)
				phiP = Math.signum(phiP)*Math.PI - phiP;
			final double delL = Math.acos(Math.sin(phi1)/Math.cos(phiP)/Math.cos(phi0) - Math.tan(phiP)*Math.tan(phi0));
			final double lamP = lam0 + Math.signum(x)*delL;
			if (Double.isNaN(phiP) || Double.isNaN(lamP)) 	return null;
			return new double[] {phiP, lamP};
		}
	};
	
	
	public static final Projection TWO_POINT_EQUIDISTANT =
			new Projection(
					"Two-point Equidistant",
					"A map that preserves distances, but not azimuths, to two arbitrary points",
					1., 0b1111, Type.QUASIAZIMUTHAL, Property.EQUIDISTANT,
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
			this.aspectRatio = a/b;
		}
		
		public double[] project(double lat, double lon, double[] pole) {
			return project(lat, lon);
		}
		
		public double[] project(double lat0, double lon0) {
			final double d1 = dist(lat0,lon0, lat1,lon1);
			final double d2 = dist(lat0,lon0, lat2,lon2);
			final double s =Math.signum(
					Math.tan(lat0)*Math.sin(lon2-lon1) +
					Math.tan(lat1)*Math.sin(lon0-lon2) +
					Math.tan(lat2)*Math.sin(lon1-lon0));
			return new double[] {
					(d1*d1-d2*d2)/(2*D) * Math.PI/a,
					s*Math.sqrt(d1*d1 - Math.pow((d1*d1-d2*d2+D*D)/(2*D), 2)) * Math.PI/a };
		}
		
		public double[] inverse(double x, double y, double[] pole) {
			return inverse(x, y);
		}
		
		public double[] inverse(double x, double y) {
			final double d1 = Math.hypot(a*x + c, b*y);
			final double d2 = Math.hypot(a*x - c, b*y);
			if (d1 + d2 > 2*a) 	return null; //TODO find out why it throws a hissy fit when y=0
			final double t1 = -(Math.cos(lat1)*Math.sin(lat2) - Math.sin(lat1)*Math.cos(lat2)*Math.cos(lon1-lon2))/Math.sin(D);
			double t2 = Math.signum(lon1-lon2)*(Math.cos(d1)*Math.cos(D) - Math.cos(d2))/(Math.sin(d1)*Math.sin(D));
			if (Math.abs(t2) > 1) 	t2 = Math.signum(t2);
			final double s0 = Math.signum(lon1-lon2)*Math.signum(y);
			final double casab = t1*t2 +s0* Math.sqrt((t1*t1 - 1)*(t2*t2 - 1));
			final double s1 = Math.signum(Math.sin(Math.acos(t1)-s0*Math.acos(t2)));
			final double PHI = Math.asin(Math.sin(lat1)*Math.cos(d1) - Math.cos(lat1)*Math.sin(d1)*casab);
			final double LAM = lon1 +s1* Math.acos((Math.cos(d1) - Math.sin(lat1)*Math.sin(PHI))/(Math.cos(lat1)*Math.cos(PHI)));
			return new double[] {PHI,LAM};
		}
		
		private double dist(double lat1, double lon1, double lat2, double lon2) {
			return Math.acos(
					Math.sin(lat1)*Math.sin(lat2) +
					Math.cos(lat1)*Math.cos(lat2)*Math.cos(lon1-lon2));
		}
	};
}
