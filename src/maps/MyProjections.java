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

import maps.Projection.Property;
import maps.Projection.Type;

/**
 * All of the projections I invented, save the tetrahedral ones,
 * because those have so much in common with other tetrahedral projections.
 * 
 * @author jkunimune
 */
public class MyProjections {
	
	public static final Projection MAGNIFIER =
			new Projection("Magnifier",
					"A novelty map projection that blows up the center way out of proportion",
					1., 0b1011, Type.AZIMUTHAL, Property.POINTLESS) {
		
		public double[] project(double lat, double lon) {
			final double p = 1/2.0+lat/Math.PI;
			final double fp = 1 - 0.1*p - 0.9*Math.pow(p,7);
			final double r = Math.PI*fp;
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon) };
		}
		
		public double[] inverse(double x, double y) {
			double R = Math.hypot(x, y);
			if (R <= 1)
				return new double[] {
						Math.PI/2 * (1 - R*.2 - R*R*R*1.8),
						Math.atan2(y, x) + Math.PI/2};
			else
				return null;
		}
	};
	
	
	public static final Projection EXPERIMENT =
			new Projection("Experiment", "What happens when you apply a complex differentiable function to a stereographic projection?",
					1., 0b0000, Type.OTHER, Property.CONFORMAL) {
		
		public double[] project(double lat, double lon) {
			final double wMag = Math.tan(Math.PI/4-lat/2);
			final Complex w = new Complex(wMag*Math.sin(lon), -wMag*Math.cos(lon));
			Complex z = w.asin();
			if (z.isInfinite() || z.isNaN())	z = new Complex(0);
			return new double[] { z.getReal(), z.getImaginary() };
		}
		
		public double[] inverse(double x, double y) {
			Complex z = new Complex((x+1)*Math.PI, y*Math.PI);
			Complex ans = z.sin();
			double p = 2 * Math.atan(ans.abs());
			double theta = ans.getArgument() - Math.PI/2;
			double lambda = Math.PI/2 - p;
			return new double[] {lambda, theta};
		}
	};
	
	
	
	public static final Projection PSEUDOSTEREOGRAPHIC =
			new Projection("Pseudostereographic", "The logical next step after Aitoff and Hammer",
					2, 0b1111, Type.PSEUDOAZIMUTHAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			final double a = Math.PI - Math.acos(Math.cos(lat)*Math.cos(lon/2));
			final double b = Math.acos(Math.sin(lat)/Math.sin(a));
			return new double[] {
					Math.PI*Math.signum(lon)/Math.tan(a/2)*Math.sin(b),
					Math.PI/2/Math.tan(a/2)*Math.cos(b) };
		}
		
		public double[] inverse(double x, double y) {
			double[] transverse = obliquifyPlnr(
					Azimuthal.STEREOGRAPHIC.inverse(x/2, y/2), new double[] {0,0,0});
			if (transverse == null) 	return null;
			else 	return new double[] {transverse[0], 2*transverse[1]};
		}
	};
	
	
	public static final Projection HYPERELLIPOWER =
			new Projection("Hyperellipower", "A parametric projection that I'm still testing",
					2., 0b1111, Type.PSEUDOCYLINDRICAL, Property.COMPROMISE,
					new String[] {"k","n","a"},
					new double[][] {{1,5,4.99},{.5,2.,1.20},{.5,2.,1.13}}) {
		
		private double k, n, a;
		
		public void setParameters(double... params) {
			this.k = params[0];
			this.n = params[1];
			this.a = params[2];
			this.aspectRatio = 2*Math.sqrt(n)/a;
		}
		
		public double[] project(double lat, double lon) {
			final double ynorm = (1-Math.pow(1-Math.abs(lat/(Math.PI/2)), n));
			return new double[] {
					Math.pow(1 - Math.pow(ynorm, k),1/k)*lon,
					ynorm*Math.PI/2/Math.sqrt(n)*a*Math.signum(lat) };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {
					(1 - Math.pow(1-Math.abs(y), 1/n))*Math.PI/2*Math.signum(y),
					x/Math.pow(1 - Math.pow(Math.abs(y),k),1/k)*Math.PI };
		}
	};
	
	
	public static final Projection TWO_POINT_EQUALIZED =
			new Projection("Two-Point Equalized", "A parametric projection that I'm still testing",
					0, 0b1111, Type.OTHER, Property.EQUIDISTANT, new String[] {"Width"},
					new double[][] { {0, 180, 120} }) {
		
		private double theta, a, b;
		
		public void setParameters(double... params) {
			theta = Math.toRadians(params[0])/2;
			this.a = Math.PI - theta; //semimajor axis
			this.b = Math.sqrt(Math.pow(Math.PI-theta, 2) - Math.pow(theta, 2)); //semiminor axis
			if (theta == 0) 				this.aspectRatio = 1;
			else if (theta == Math.PI/2) 	this.aspectRatio = 1.3;
			else 							this.aspectRatio = b/a*Math.sqrt(Math.tan(theta)/theta);
		}
		
		public double[] project(double lat, double lon) {
			if (theta == 0) 	return Azimuthal.POLAR.project(lat, lon);
			final double d1 = Math.acos(
					Math.sin(lat)*Math.cos(theta) - Math.cos(lat)*Math.sin(theta)*Math.cos(lon));
			final double d2 = Math.acos(
					Math.sin(lat)*Math.cos(theta) + Math.cos(lat)*Math.sin(theta)*Math.cos(lon));
			final double k = Math.signum(lon)*Math.sqrt(Math.tan(theta)/theta);
			final double s = Math.PI/(Math.PI-theta/2);
			return new double[] {
					s*k*Math.sqrt(d1*d1 - Math.pow((d1*d1-d2*d2+4*theta*theta)/(4*theta), 2)),
					s*(d2*d2-d1*d1)/(4*theta) };
		}
		
		public double[] inverse(double x, double y) {
			if (theta == 0) 	return Azimuthal.POLAR.inverse(x, y);
			final double d1 = Math.hypot(b*x, a*y - theta);
			final double d2 = Math.hypot(b*x, a*y + theta);
			if (d1 + d2 > 2*a) 	return null;
			final double phi = Math.asin((Math.cos(d1)+Math.cos(d2))/(2*Math.cos(theta)));
			final double lam = Math.signum(x)*
					Math.acos((Math.sin(phi)*Math.cos(theta) - Math.cos(d1))/(Math.cos(phi)*Math.sin(theta)));
			return new double[] { phi, lam };
		}
	};
}
