/**
 * MIT License
 * 
 * Copyright (c) 2018 Justin Kunimune
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

import utils.BoundingBox;
import utils.NumericalAnalysis;

import static java.lang.Double.isNaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.signum;

/**
 * Some parameterized pseudocylindrical/lenticular projections that have been optimized to minimize error
 * 
 * @author Justin Kunimune
 */
public class Gyorffy {
	
	private static final double PI2M2 = pow(PI/2, -2);

	public static final Projection B = new PolynomialProjection(
			"B", "The optimal pseudocylindrical projection.", 2,
			new double[] {0.75762, 2.00000, 4.63375, 0.00264, 1, 0, 0});

	public static final Projection D = new PolynomialProjection(
			"D", "An optimal pointed-polar projection, with an emphasis on polar regions.", 3,
			new double[] {0.71416, 3.79209, 2, 0.00902, 0.87550, 0.01004, 0.00273});

	public static final Projection E = new PolynomialProjection(
			"E", "The optimal pointed-polar projection.", 3,
			new double[] {0.74532, 2, 4.04753, 0.00730, 0.93884, 0.00271, 0.00450});
	
	public static final Projection F = new PolynomialProjection(
			"F", "The optimal pointed-polar projection, with a rounded outline.", 4,
			new double[] {0.77172, 2, 3.26655, 0.00649, 0.88525, 0.00950, 0.00305});
	
	
	private static class PolynomialProjection extends Projection {
		
		private final double[] coefs;
		
		protected PolynomialProjection(
				String letter, String description,
				int rating, double[] coefs) {
			super("Gy\u00F6rffy "+letter, description,
			      new BoundingBox(2*coefs[0]*(PI + coefs[3]*pow(PI, 3)), PI),
					0b1011, (coefs[4]==1 && coefs[5]==0 && coefs[6]==0) ? Type.PSEUDOCYLINDRICAL : Type.OTHER,
					Property.COMPROMISE, rating);
			this.coefs = coefs;
		}

		public double[] project(double lat, double lon) {
			return new double[] {x(lat, lon, coefs), y(lat, lon, coefs)};
		}

		public double[] inverse(double x, double y) {
			double phi0 = y;
			double lam0 = x/pow(1 - pow(y/(PI/2), 2), 1/3.);
			double[] res = NumericalAnalysis.newtonRaphsonApproximation(
					x, y, phi0, lam0, Gyorffy::x, Gyorffy::y,
					Gyorffy::dxdp, Gyorffy::dxdl, Gyorffy::dydp, Gyorffy::dydl,
					1e-4, this.coefs); // this converges surprisingly well atc
			if (res == null || isNaN(res[1]) ||
					(abs(res[1]) < PI && abs(x) > x(y, coefs)))
				return null; // it does have a nasty habit of thinking it's converged outside the map, though
			return res;
		}
	}
	
	private static double x(double y, double[] c) {
		return c[0]*pow(1 - pow(2*abs(y)/PI, c[1]), 1/c[2])*(PI + c[3]*pow(PI, 3));
	}
	
	private static double x(double phi, double lam, double[] c) {
		double lam3 = pow(lam, 3);
		return c[0]*pow(1 - pow(2*abs(y(phi, lam, c))/PI, c[1]), 1/c[2])*(lam + c[3]*lam3);
	}
	
	private static double y(double phi, double lam, double[] c) {
		double phi3 = pow(phi, 3);
		double lam2 = pow(lam, 2);
		double lam4 = pow(lam, 4);
		return c[4]*phi + (1-c[4])*PI2M2*phi3 + (c[5]*lam2 + c[6]*lam4)*(phi - PI2M2*phi3);
	}
	
	private static double dxdp(double phi, double lam, double[] c) {
		double lam3 = pow(lam, 3);
		return -2/PI*c[0]/c[2]*c[1]*pow(2*abs(y(phi, lam, c))/PI, c[1]-1)*
		       pow(1 - pow(2*abs(y(phi, lam, c))/PI, c[1]), 1/c[2]-1)*
		       (lam + c[3]*lam3)*signum(phi)*dydp(phi, lam, c);
	}
	
	private static double dxdl(double phi, double lam, double[] c) {
		double lam2 = pow(lam, 2);
		double lam3 = pow(lam, 3);
		return -2/PI*c[0]/c[2]*c[1]*pow(2*abs(y(phi, lam, c))/PI, c[1]-1)*
		       pow(1 - pow(2*abs(y(phi, lam, c))/PI, c[1]), 1/c[2]-1)*
		       (lam + c[3]*lam3)*signum(phi)*dydl(phi, lam, c) +
		       c[0]*pow(1 - pow(2*abs(y(phi, lam, c))/PI, c[1]), 1/c[2])*(1 + 3*c[3]*lam2);
	}
	
	private static double dydp(double phi, double lam, double[] c) {
		double phi2 = pow(phi, 2);
		double lam2 = pow(lam, 2);
		double lam4 = pow(lam, 4);
		return c[4] + (1-c[4])*3*PI2M2*phi2 + (c[5]*lam2 + c[6]*lam4)*(1 - 3*PI2M2*phi2);
	}
	
	private static double dydl(double phi, double lam, double[] c) {
		double lam3 = pow(lam, 3);
		double phi3 = pow(phi, 3);
		return (2*c[5]*lam + 4*c[6]*lam3)*(phi - PI2M2*phi3);
	}
	
}
