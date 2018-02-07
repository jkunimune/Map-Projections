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
import utils.Math2;

/**
 * All of the projections I invented, save the tetrahedral ones, because
 * those have so much in common with other tetrahedral projections.
 * 
 * @author jkunimune
 */
public class MyProjections {
	
	public static final Projection MAGNIFIER =
			new Projection("Magnifier",
					"A novelty map projection that blows up the center way out of proportion.",
					2, 2, 0b1011, Type.AZIMUTHAL, Property.POINTLESS, 0) {
		
		public double[] project(double lat, double lon) {
			final double p = 1/2.0+lat/Math.PI;
			final double fp = 1 - 0.1*p - 0.9*Math.pow(p,7);
			return new double[] { fp*Math.sin(lon), -fp*Math.cos(lon) };
		}
		
		public double[] inverse(double x, double y) {
			double R = Math.hypot(x, y);
			if (R <= 1)
				return new double[] {
						Math.PI/2 * (1 - R*.2 - R*R*R*1.8),
						Math.atan2(x, -y)};
			else
				return null;
		}
	};
	
	
	public static final Projection EXPERIMENT = new Projection(
			"Complex Arcsine", "What happens when you apply a complex differentiable function to a stereographic projection?",
			6, 6, 0b0000, Type.OTHER, Property.CONFORMAL, 0) {
		
		public double[] project(double lat, double lon) {
			final double wMag = Math.tan(Math.PI/4-lat/2);
			final Complex w = new Complex(wMag*Math.sin(lon), -wMag*Math.cos(lon));
			Complex z = w.asin();
			if (z.isInfinite() || z.isNaN())	z = new Complex(0);
			return new double[] { z.getReal(), z.getImaginary() };
		}
		
		public double[] inverse(double x, double y) {
			Complex z = new Complex(x, y);
			Complex ans = z.sin();
			double p = 2 * Math.atan(ans.abs());
			double theta = Math2.coerceAngle(ans.getArgument());
			double lambda = Math.PI/2 - p;
			if (x < -Math.PI/2)
				return new double[] {lambda, theta-2*Math.PI};
			else if (x < Math.PI/2)
				return new double[] {lambda, theta};
			else
				return new double[] {lambda, theta+2*Math.PI};
		}
	};
	
	
	
	public static final Projection PSEUDOSTEREOGRAPHIC = new Projection(
			"Pseudostereographic", "The next logical step after Aitoff and Hammer.",
			2, 1, 0b1111, Type.PSEUDOAZIMUTHAL, Property.COMPROMISE, 1) {
		
		public double[] project(double lat, double lon) {
			final double a = Math.PI - Math.acos(Math.cos(lat)*Math.cos(lon/2));
			final double b = Math.acos(Math.sin(lat)/Math.sin(a));
			return new double[] {
					Math.signum(lon)/Math.tan(a/2)*Math.sin(b),
					Math.cos(b)/Math.tan(a/2)/2 };
		}
		
		public double[] inverse(double x, double y) {
			double[] transverse = obliquifyPlnr(
					Azimuthal.STEREOGRAPHIC.inverse(x, 2*y), new double[] {0,0,0});
			if (transverse == null) 	return null;
			else 	return new double[] {transverse[0], 2*transverse[1]};
		}
	};
	
	
	public static final Projection TWO_POINT_EQUALIZED = new Projection("Two-Point Equalised",
			"A projection I invented specifically for viewing small elliptical regions of the Earth.",
			0, 0, 0b1111, Type.OTHER, Property.EQUIDISTANT, 2,
			new String[] {"Width"}, new double[][] { {0, 180, 120} }) {
		
		private double theta;
		
		public void setParameters(double... params) {
			theta = Math.toRadians(params[0])/2;
			this.height = 2*Math.PI - 2*theta; //major axis
			this.width = 2*Math.sqrt(Math.pow(Math.PI-theta, 2) - Math.pow(theta, 2)) *
					Math.sqrt(Math.tan(theta)/theta); //minor axis
			if (theta == 0)
				this.width = this.height = 2*Math.PI;
		}
		
		public double[] project(double lat, double lon) {
			if (theta == 0) 	return Azimuthal.POLAR.project(lat, lon);
			final double d1 = Math.acos(
					Math.sin(lat)*Math.cos(theta) - Math.cos(lat)*Math.sin(theta)*Math.cos(lon));
			final double d2 = Math.acos(
					Math.sin(lat)*Math.cos(theta) + Math.cos(lat)*Math.sin(theta)*Math.cos(lon));
			final double k = Math.signum(lon)*Math.sqrt(Math.tan(theta)/theta);
			return new double[] {
					k*Math.sqrt(d1*d1 - Math.pow((d1*d1-d2*d2+4*theta*theta)/(4*theta), 2)),
					(d2*d2-d1*d1)/(4*theta) };
		}
		
		public double[] inverse(double x, double y) {
			if (theta == 0) 	return Azimuthal.POLAR.inverse(x, y);
			final double d1 = Math.hypot(x/Math.sqrt(Math.tan(theta)/theta), y - theta);
			final double d2 = Math.hypot(x/Math.sqrt(Math.tan(theta)/theta), y + theta);
			if (d1 + d2 > height) 	return null;
			final double phi = Math.asin((Math.cos(d1)+Math.cos(d2))/(2*Math.cos(theta)));
			final double lam = Math.signum(x)*Math.acos(
					(Math.sin(phi)*Math.cos(theta) - Math.cos(d1))/(Math.cos(phi)*Math.sin(theta)));
			return new double[] { phi, lam };
		}
	};
}
