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
import utils.Shape;

import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.cos;
import static java.lang.Math.hypot;
import static java.lang.Math.pow;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.toRadians;

/**
 * All of the projections I invented, save the tetrahedral ones, because
 * those have so much in common with other tetrahedral projections.
 * 
 * @author jkunimune
 */
public class MyProjections {

	public static final Projection TWO_POINT_EQUALIZED = new Projection("Two-Point Equalised",
			"A projection I invented specifically for viewing small elliptical regions of the Earth.",
			null, 0b1111, Type.OTHER, Property.EQUIDISTANT, 2,
			new String[] {"Width"}, new double[][] { {0, 180, 120} }) {
		
		private double theta;
		
		public void initialize(double... params) {
			theta = toRadians(params[0])/2;
			if (theta != 0)
				this.shape = Shape.ellipse(
						sqrt(pow(PI - theta, 2) - pow(theta, 2))*
						sqrt(tan(theta)/theta), //minor axis
						PI - theta //major axis
				);
			else
				this.shape = Shape.circle(PI);
		}
		
		public double[] project(double lat, double lon) {
			if (theta == 0) 	return Azimuthal.POLAR.project(lat, lon);
			final double d1 = acos(
					sin(lat)*cos(theta) - cos(lat)*sin(theta)*cos(lon));
			final double d2 = acos(
					sin(lat)*cos(theta) + cos(lat)*sin(theta)*cos(lon));
			final double k = signum(lon)*sqrt(tan(theta)/theta);
			return new double[] {
					k*sqrt(d1*d1 - pow((d1*d1-d2*d2+4*theta*theta)/(4*theta), 2)),
					(d2*d2-d1*d1)/(4*theta) };
		}
		
		public double[] inverse(double x, double y) {
			if (theta == 0) 	return Azimuthal.POLAR.inverse(x, y);
			final double d1 = hypot(x/sqrt(tan(theta)/theta), y - theta);
			final double d2 = hypot(x/sqrt(tan(theta)/theta), y + theta);
			if (d1 + d2 > 2*shape.yMax) 	return null;
			final double phi = asin((cos(d1)+cos(d2))/(2*cos(theta)));
			final double lam = signum(x)*acos(
					(sin(phi)*cos(theta) - cos(d1))/(cos(phi)*sin(theta)));
			return new double[] { phi, lam };
		}
	};
}
