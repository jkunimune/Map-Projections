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

import maps.Projection.Type;
import utils.NumericalAnalysis;

/**
 * TODO: Write description
 * 
 * @author Justin Kunimune
 */
public class Gyorffy {
	
	private static final double PI2M2 = Math.pow(Math.PI/2, -2);
	
	public static final Projection C = new PolynomialProjection(
			"C", "The optimal pseudocylindrical projection.", Type.PSEUDOCYLINDRICAL, 3,
			new double[] {0.76158, 1.67084, 5.17538, 0.00272, 1, 0, 0});
	
	public static final Projection D = new PolynomialProjection(
			"D", "An optimal pointed-polar projection, with an emphasis on polar regions.", Type.OTHER, 3,
			new double[] {0.71416, 3.79209, 2, 0.00902, 0.87550, 0.01004, 0.00273});
	
	public static final Projection E = new PolynomialProjection(
			"E", "The optimal pointed-polar projection.", Type.OTHER, 3,
			new double[] {0.74532, 2, 4.04753, 0.00730, 0.93884, 0.00271, 0.00450});
	
	public static final Projection F = new PolynomialProjection(
			"F", "The optimal pointed-polar projection, with a rounded outline.", Type.OTHER, 4,
			new double[] {0.77172, 2, 3.26655, 0.00649, 0.88525, 0.00950, 0.00305});
	
	
	private static class PolynomialProjection extends Projection {
		
		private double[] coefs;
		
		protected PolynomialProjection(
				String letter, String description, Type type,
				int rating, double[] coefs) {
			super("Gy\u00F6rffy "+letter, description, 2*coefs[0]*(Math.PI + coefs[3]*Math.pow(Math.PI, 3)), Math.PI,
					0b1011, (coefs[4]==1 && coefs[5]==0 && coefs[6]==0) ? Type.PSEUDOCYLINDRICAL : Type.OTHER,
					Property.COMPROMISE, rating);
			this.coefs = coefs;
		}

		public double[] project(double lat, double lon) {
			return new double[] {x(lat, lon, coefs), y(lat, lon, coefs)};
		}

		public double[] inverse(double x, double y) {
			double phi0 = y;
			double lam0 = x/Math.pow(1 - Math.pow(y/(Math.PI/2), 2), 1/3.);
			double[] res = NumericalAnalysis.newtonRaphsonApproximation(
					x, y, phi0, lam0, Gyorffy::x, Gyorffy::y,
					Gyorffy::dxdp, Gyorffy::dxdl, Gyorffy::dydp, Gyorffy::dydl,
					1e-3, this.coefs); // this converges surprisingly well atc
			if (res == null || Double.isNaN(res[1]) ||
					(Math.abs(res[1]) < Math.PI && Math.abs(x) > x(y, coefs)))
				return null; // it does have a nasty habit of thinking it's converged outside the map, though
			return res;
		}
	}
	
	private static final double x(double y, double[] c) {
		return c[0]*Math.pow(1 - Math.pow(2*Math.abs(y)/Math.PI, c[1]), 1/c[2])*(Math.PI + c[3]*Math.pow(Math.PI, 3));
	}
	
	private static final double x(double phi, double lam, double[] c) {
		double lam3 = Math.pow(lam, 3);
		return c[0]*Math.pow(1 - Math.pow(2*Math.abs(y(phi, lam, c))/Math.PI, c[1]), 1/c[2])*(lam + c[3]*lam3);
	}
	
	private static final double y(double phi, double lam, double[] c) {
		double phi3 = Math.pow(phi, 3);
		double lam2 = Math.pow(lam, 2);
		double lam4 = Math.pow(lam, 4);
		return c[4]*phi + (1-c[4])*PI2M2*phi3 + (c[5]*lam2 + c[6]*lam4)*(phi - PI2M2*phi3);
	}
	
	private static final double dxdp(double phi, double lam, double[] c) {
		double lam3 = Math.pow(lam, 3);
		return -2/Math.PI*c[0]/c[2]*c[1]*Math.pow(2*Math.abs(y(phi, lam, c))/Math.PI, c[1]-1)*Math.pow(1 - Math.pow(2*Math.abs(y(phi, lam, c))/Math.PI, c[1]), 1/c[2]-1)*(lam + c[3]*lam3)*Math.signum(phi)*dydp(phi, lam, c);
	}
	
	private static final double dxdl(double phi, double lam, double[] c) {
		double lam2 = Math.pow(lam, 2);
		double lam3 = Math.pow(lam, 3);
		return -2/Math.PI*c[0]/c[2]*c[1]*Math.pow(2*Math.abs(y(phi, lam, c))/Math.PI, c[1]-1)*Math.pow(1 - Math.pow(2*Math.abs(y(phi, lam, c))/Math.PI, c[1]), 1/c[2]-1)*(lam + c[3]*lam3)*Math.signum(phi)*dydl(phi, lam, c)
				+ c[0]*Math.pow(1 - Math.pow(2*Math.abs(y(phi, lam, c))/Math.PI, c[1]), 1/c[2])*(1 + 3*c[3]*lam2);
	}
	
	private static final double dydp(double phi, double lam, double[] c) {
		double phi2 = Math.pow(phi, 2);
		double lam2 = Math.pow(lam, 2);
		double lam4 = Math.pow(lam, 4);
		return c[4] + (1-c[4])*3*PI2M2*phi2 + (c[5]*lam2 + c[6]*lam4)*(1 - 3*PI2M2*phi2);
	}
	
	private static final double dydl(double phi, double lam, double[] c) {
		double lam3 = Math.pow(lam, 3);
		double phi3 = Math.pow(phi, 3);
		return (2*c[5]*lam + 4*c[6]*lam3)*(phi - PI2M2*phi3);
	}
	
}
