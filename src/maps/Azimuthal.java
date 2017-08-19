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

public class Azimuthal {
	
	public static final Projection STEREOGRAPHIC =
			new Projection("Stereographic", 1., 0b0111, Type.AZIMUTHAL, Property.CONFORMAL,
					"mathematically important") {
		
		public double[] project(double lat, double lon) {
			final double r = Math.PI/2/(Math.tan(lat/2 + Math.PI/4));
			return new double[] {r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.PI/2 - 2*Math.atan(2*Math.hypot(x, y)),
					Math.atan2(y, x) + Math.PI/2 };
		}
	};
	
	
	public static final Projection POLAR =
			new Projection("Polar", 1., 0b1111, Type.AZIMUTHAL, Property.EQUIDISTANT) {
		
		public double[] project(double lat, double lon) {
			final double r = Math.PI/2 - lat;
			return new double[] {r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			double phi = Math.PI/2 - Math.PI * Math.hypot(x, y);
			if (phi > -Math.PI/2)
				return new double[] {phi, Math.atan2(y, x) + Math.PI/2};
			else
				return null;
		}
	};
	
	
	public static final Projection EQUAL_AREA =
			new Projection(
					"Azimuthal Equal-Area", 1., 0b1111, Type.AZIMUTHAL, Property.EQUAL_AREA) {
		
		public double[] project(double lat, double lon) {
			final double r = Math.PI*Math.cos((Math.PI/2+lat)/2);
			return new double[] {r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			double R = Math.hypot(x, y);
			if (R <= 1)
				return new double[] {Math.asin(1-2*R*R), Math.atan2(y,x)+Math.PI/2};
			else
				return null;
		}
	};
	
	
	public static final Projection ORTHOGRAPHIC =
			new Projection("Orthographic", "A projection that mimics the Earth viewed from a great distance",
					1., 0b1111, Type.AZIMUTHAL, Property.PERSPECTIVE) {
		
		public double[] project(double lat, double lon) {
			if (lat < 0)	lat = 0;
			final double r = Math.PI*Math.cos(lat);
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon) };
		}
		
		public double[] inverse(double x, double y) {
			double R = Math.hypot(x, y);
			if (R <= 1)
				return new double[] { Math.acos(R), Math.atan2(y, x) + Math.PI/2 };
			else
				return null;
		}
	};
	
	
	public static final Projection GNOMONIC =
			new Projection(
					"Gnomonic", "A projection that draws all great circles as straight lines",
					1., 0b0111, Type.AZIMUTHAL, Property.GNOMONIC) {
		
		public double[] project(double lat, double lon) {
			if (lat <= 0)	lat = 1e-5;
			final double r = Math.tan(Math.PI/2 - lat);
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.PI/2 - Math.atan(2*Math.hypot(x, y)),
					Math.atan2(y, x) + Math.PI/2 };
		}
	};
	
	
	public static final Projection PERSPECTIVE =
			new Projection(
					"Perspective", "A projection that mimics the actual appearance of the Earth",
					1., 0b1111, Type.AZIMUTHAL, Property.PERSPECTIVE,
					new String[] {"Percentage"}, new double[][] {{1,99}}) {
		
		private double amount;
		
		public void setParameters(double... params) {
			this.amount = params[0]/100.;
		}

		public double[] project(double lat, double lon) {
			// TODO: Projection wishlist
			return null;
		}

		public double[] inverse(double x, double y) {
			// TODO: Projection wishlist
			return null;
		}
	};
}
