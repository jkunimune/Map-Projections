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
			new Projection(
					"Stereographic", 4, 4, 0b0111, Type.AZIMUTHAL, Property.CONFORMAL, 2,
					"mathematically important") {
		
		public double[] project(double lat, double lon) {
			if (lat < -1.5) 	lat = -1.5;
			final double r = 1/(Math.tan(lat/2 + Math.PI/4));
			return new double[] {r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.PI/2 - 2*Math.atan(Math.hypot(x, y)), Math.atan2(x, -y) };
		}
	};
	
	
	public static final Projection POLAR =
			new Projection(
					"Azimuthal Equidistant", 2*Math.PI, 2*Math.PI, 0b1111,
					Type.AZIMUTHAL, Property.EQUIDISTANT, 2) {
		
		public double[] project(double lat, double lon) {
			final double r = Math.PI/2 - lat;
			return new double[] {r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			double phi = Math.PI/2 - Math.hypot(x, y);
			if (phi > -Math.PI/2)
				return new double[] {phi, Math.atan2(x, -y)};
			else
				return null;
		}
	};
	
	
	public static final Projection EQUAL_AREA =
			new Projection(
					"Azimuthal Equal-Area", 2, 2, 0b1111, Type.AZIMUTHAL, Property.EQUAL_AREA, 1) {
		
		public double[] project(double lat, double lon) {
			final double r = Math.cos((Math.PI/2+lat)/2);
			return new double[] {r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			double r = Math.hypot(x, y);
			if (r <= 1)
				return new double[] { Math.asin(1-2*r*r), Math.atan2(x, -y) };
			else
				return null;
		}
	};
	
	
	public static final Projection GNOMONIC =
			new Projection(
					"Gnomonic", "A projection that draws all great circles as straight lines.",
					4, 4, 0b0111, Type.AZIMUTHAL, Property.GNOMONIC, 2) {
		
		public double[] project(double lat, double lon) {
			if (lat < 0.01) 	lat = 0.01;
			final double r = Math.tan(Math.PI/2 - lat);
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.PI/2 - Math.atan(Math.hypot(x, y)), Math.atan2(x, -y) };
		}
	};
	
	
	public static final Projection ORTHOGRAPHIC =
			new Projection(
					"Orthographic", "A projection that mimics the Earth viewed from a great distance.",
					2, 2, 0b0111, Type.AZIMUTHAL, Property.PERSPECTIVE, 3) {
		
		public double[] project(double lat, double lon) {
			if (lat < 0)	lat = 0;
			return new double[] { Math.cos(lat)*Math.sin(lon), -Math.cos(lat)*Math.cos(lon) };
		}
		
		public double[] inverse(double x, double y) {
			double R = Math.hypot(x, y);
			if (R <= 1)
				return new double[] { Math.acos(R), Math.atan2(x, -y) };
			else
				return null;
		}
	};
	
	
	public static final Projection PERSPECTIVE =
			new Projection(
					"Perspective", "A projection that mimics the actual appearance of the Earth.",
					0, 0, 0b0111, Type.AZIMUTHAL, Property.PERSPECTIVE, 4,
					new String[] {"Percentage"}, new double[][] {{1,99,33.3}}) {
		
		private double d; //viewing distance in sphere radii
		
		public void setParameters(double... params) {
			this.d = 1/(1 - 2*params[0]/100);
			this.width = this.height = 2/Math.sqrt(d*d-1);
		}
		
		public double[] project(double lat, double lon) {
			if (Double.isInfinite(d)) 	return ORTHOGRAPHIC.project(lat, lon);
			if (lat < Math.asin(1/d)) 	lat = Math.asin(1/d);
			final double r = Math.abs(Math.cos(lat)/(d - Math.sin(lat)));
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon) };
		}
		
		public double[] inverse(double x, double y) {
			if (Double.isInfinite(d)) 	return ORTHOGRAPHIC.inverse(x, y);
			final double h = Math.hypot(x, y);
			if (h > this.width/2) 	return null;
			final double theta = Math.atan(h);
			final double phi = Math.acos(d*Math.sin(theta)) + theta;
			if (phi < Math.PI/2)
				return new double[] { phi, Math.atan2(x, -y) };
			else
				return new double[] { Math.PI - phi, Math.atan2(x, -y) };
		}
	};
}
