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
import utils.NumericalAnalysis;
import utils.Shape;

import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;

/**
 * All the useful Winkel Tripel equations and derivatives.
 * I tried solving for these equations myself, and I think I got them mostly
 * right, but the expressions were just too complicated. I got better results
 * by transcribing the below equations from Ipb&uuml;ker and Bildirici's paper,
 * "A General Algorithm for the Inverse Transformation of Map Projections Using Jacobian Matrices"
 * 
 * Ipb&uuml;ker, Cengizhan; Bildirici, I.&Ouml;ztug (2002). "A General Algorithm for the
 *  	Inverse Transformation of Map Projections Using Jacobian Matrices".
 *  	Proceedings of the Third International Symposium Mathematical &amp;
 *  	Computational Applications. Third International Symposium Mathematical &amp;
 *  	Computational Applications September 4-6, 2002. Konya, Turkey. Selcuk,
 *  	Turkey. pp. 175-182. Archived from the original on 20 October 2014.
 * 
 * @author jkunimune
 */
public final class WinkelTripel {
	
	public static final Projection WINKEL_TRIPEL =
			new Projection("Winkel Tripel", "Oswald Winkel", "National Geographic's compromise projection of choice",
					null, true, true, true, false, Type.OTHER, Property.COMPROMISE, 3,
					new String[] {"Std. Parallel"}, new double[][] {{0, 90, toDegrees(acos(2/PI))}}) {
		
		private double stdParallel;
		
		public void initialize(double... params) {
			this.stdParallel = toRadians(params[0]);
			this.shape = Shape.meridianEnvelope(WINKEL_TRIPEL);
		}
		
		public double[] project(double lat, double lon) {
			return new double[] { f1pX(lat,lon), f2pY(lat,lon) };
		}
		
		public double[] inverse(double x, double y) {
			return NumericalAnalysis.newtonRaphsonApproximation(
					x, y,
					y/2, x*(1 + cos(y*PI/2))/(2 + 2*cos(stdParallel)), //initial guess is Eckert V
					this::f1pX, this::f2pY,
					this::df1dphi, this::df1dlam, this::df2dphi, this::df2dlam, .001);
		}
		
		private double f1pX(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return 2*d/sqrt(c)*cos(phi)*sin(lam/2) + lam*cos(stdParallel);
		}
		
		private double f2pY(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return d/sqrt(c)*sin(phi) + phi;
		}
		
		private double df1dphi(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return sin(lam)*sin(2*phi)/(4*c) - d/pow(c,1.5)*sin(phi)*sin(lam/2);
		}
		
		private double df1dlam(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return pow(cos(phi)*sin(lam/2), 2)/c + d/pow(c,1.5)*cos(phi)*cos(lam/2)*pow(sin(phi),2) + cos(stdParallel);
		}
		
		private double df2dphi(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return pow(sin(phi),2)*cos(lam/2)/c + d/pow(c,1.5)*(1-pow(cos(lam/2),2))*cos(phi) + 1;
		}
		
		private double df2dlam(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return (sin(2*phi)*sin(lam/2)/c - d/pow(c,1.5)*sin(phi)*pow(cos(phi),2)*sin(lam))/4.0;
		}
		
		private double D(double phi, double lam) {
			return acos(cos(phi)*cos(lam/2));
		}
		
		private double C(double phi, double lam) {
			if (phi == 0 && lam == 0) 	return 1; //there's a hole here
			return 1 - pow(cos(phi)*cos(lam/2), 2);
		}
	};
}
