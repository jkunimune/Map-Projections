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
			new Projection("Winkel Tripel", "National Geographic's compromise projection of choice",
					0, 0b1011, Type.OTHER, Property.COMPROMISE,
					new String[] {"Std. Parallel"}, new double[][] {{0,90,50.4598}}) {
		
		public void setParameters(double... params) {
			this.aspectRatio = 1 + Math.cos(Math.toRadians(params[0]));
		}
		
		public double[] project(double lat, double lon) {
			return new double[] { f1pX(lat,lon)/aspectRatio, f2pY(lat,lon)/aspectRatio };
		}
		
		public double[] inverse(double x, double y) {
			return NumericalAnalysis.newtonRaphsonApproximation(
					x*aspectRatio, y, y/2, x*(1 + Math.cos(y*Math.PI/2))/2,
					this::f1pX, this::f2pY,
					this::df1dphi, this::df1dlam, this::df2dphi, this::df2dlam, .002);
		}
		
		private double f1pX(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return 2/Math.PI*d/Math.sqrt(c)*Math.cos(phi)*Math.sin(lam/2)
					+ lam/Math.PI*(aspectRatio-1);
		}
		
		private double f2pY(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return d/Math.sqrt(c)*Math.sin(phi)/Math.PI + phi/Math.PI;
		}
		
		private double df1dphi(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return Math.sin(lam)*Math.sin(2*phi)/(4*c) - d/Math.pow(c,1.5)*Math.sin(phi)*Math.sin(lam/2);
		}
		
		private double df1dlam(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return Math.pow(Math.cos(phi)*Math.sin(lam/2), 2)/c + d/Math.pow(c,1.5)*Math.cos(phi)*Math.cos(lam/2)*Math.pow(Math.sin(phi),2) + aspectRatio-1;
		}
		
		private double df2dphi(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return Math.pow(Math.sin(phi),2)*Math.cos(lam/2)/c + d/Math.pow(c,1.5)*(1-Math.pow(Math.cos(lam/2),2))*Math.cos(phi) + 1;
		}
		
		private double df2dlam(double phi, double lam) {
			final double d = D(phi,lam);
			final double c = C(phi,lam);
			return (Math.sin(2*phi)*Math.sin(lam/2)/c - d/Math.pow(c,1.5)*Math.sin(phi)*Math.pow(Math.cos(phi),2)*Math.sin(lam))/4.0;
		}
		
		private double D(double phi, double lam) {
			return Math.acos(Math.cos(phi)*Math.cos(lam/2));
		}
		
		private double C(double phi, double lam) {
			return 1 - Math.pow(Math.cos(phi)*Math.cos(lam/2), 2);
		}
	};
}
