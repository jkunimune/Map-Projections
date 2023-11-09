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
import utils.BoundingBox;

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
			theta = Math.toRadians(params[0])/2;
			if (theta != 0)
				this.bounds = new BoundingBox(
						2*Math.PI - 2*theta, //major axis
				        2*Math.sqrt(Math.pow(Math.PI-theta, 2) - Math.pow(theta, 2)) *
				        Math.sqrt(Math.tan(theta)/theta) //minor axis
				);
			else
				this.bounds = new BoundingBox(2*Math.PI, 2*Math.PI);
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
			if (d1 + d2 > 2*bounds.yMax) 	return null;
			final double phi = Math.asin((Math.cos(d1)+Math.cos(d2))/(2*Math.cos(theta)));
			final double lam = Math.signum(x)*Math.acos(
					(Math.sin(phi)*Math.cos(theta) - Math.cos(d1))/(Math.cos(phi)*Math.sin(theta)));
			return new double[] { phi, lam };
		}
	};
}
