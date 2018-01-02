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

/**
 * TODO: Write description
 * 
 * @author jkunimune
 */
public class Lenticular {
	
	public static final Projection AITOFF =
			new Projection("Aitoff", "A compromise projection shaped like an ellipse.",
					2*Math.PI, Math.PI, 0b1111, Type.PSEUDOAZIMUTHAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			final double a = Math.acos(Math.cos(lat)*Math.cos(lon/2));
			return new double[] {
					2*Math.cos(lat)*Math.sin(lon/2)*a/Math.sin(a),
					Math.sin(lat)*a/Math.sin(a)};
		}
		
		public double[] inverse(double x, double y) {
			final double[] intermediate = Azimuthal.POLAR.inverse(x/2, y);
			double[] transverse = obliquifyPlnr(intermediate, new double[] {0,0,0});
			if (transverse != null) 	transverse[1] *= 2;
			return transverse;
		}
	};
	
	
	public static final Projection HAMMER =
			new Projection("Hammer", "An equal-area projection shaped like an ellipse.",
					4, 2, 0b1111, Type.PSEUDOAZIMUTHAL, Property.EQUAL_AREA) {
		
		public double[] project(double lat, double lon) {
			return new double[] {
					2*Math.cos(lat)*Math.sin(lon/2)/Math.sqrt(1+Math.cos(lat)*Math.cos(lon/2)),
					Math.sin(lat)/Math.sqrt(1+Math.cos(lat)*Math.cos(lon/2)) };
		}
		
		public double[] inverse(double x, double y) {
			final double z = Math.sqrt(1 - x*x/8 - y*y/2);
			final double shift = (Math.hypot(x/2, y) > 1) ? 2*Math.PI*Math.signum(x) : 0;
			return new double[] {
					Math.asin(z*y*Math.sqrt(2)),
					2*Math.atan(Math.sqrt(.5)*z*x / (2*z*z - 1)) + shift};
		}
	};
	
	
	public static final Projection VAN_DER_GRINTEN =
			new Projection(
					"Van der Grinten", "A circular compromise map that is popular for some reason.",
					2, 2, 0b1111, Type.OTHER, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			if (lat == 0) //special case 1: equator
				return new double[] {lon/Math.PI, 0};
			if (lon == 0 || lat >= Math.PI/2 || lat <= -Math.PI/2) //special case 3: prime meridian
				return new double[] {0, Math.tan(Math.asin(2*lat/Math.PI)/2)};
			
			final double t = Math.abs(Math.asin(2*lat/Math.PI));
			final double A = Math.abs(Math.PI/lon - lon/Math.PI)/2;
			final double G = Math.cos(t)/(Math.sin(t)+Math.cos(t)-1);
			final double P = G*(2/Math.sin(t) - 1);
			final double Q = A*A + G;
			return new double[] {
					Math.signum(lon)*(A*(G-P*P)+Math.sqrt(A*A*(G-P*P)*(G-P*P)-(P*P+A*A)*(G*G-P*P)))/(P*P+A*A),
					Math.signum(lat)*(P*Q-A*Math.sqrt((A*A+1)*(P*P+A*A)-Q*Q))/(P*P+A*A)};
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
	
}
