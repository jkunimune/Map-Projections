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

public class Azimuthal {
	
	public static final Projection STEREOGRAPHIC = new Projection(
			"Stereographic", new BoundingBox(4, 4), 0b0111, Type.AZIMUTHAL, Property.CONFORMAL, 2,
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
	
	
	public static final Projection POLAR = new Projection(
			"Azimuthal Equidistant", new BoundingBox(2*Math.PI, 2*Math.PI), 0b1111, Type.AZIMUTHAL,
			Property.EQUIDISTANT, 2) {
		
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
	
	
	public static final Projection EQUAL_AREA = new Projection(
			"Azimuthal Equal-Area", new BoundingBox(2, 2), 0b1111, Type.AZIMUTHAL, Property.EQUAL_AREA, 1) {
		
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
	
	
	public static final Projection GNOMONIC = new Projection(
			"Gnomonic", "A projection that draws all great circles as straight lines.",
			new BoundingBox(4, 4), 0b0111, Type.AZIMUTHAL, Property.GNOMONIC, 2) {
		
		public double[] project(double lat, double lon) {
			if (lat < 0.2) 	lat = 0.2;
			final double r = Math.tan(Math.PI/2 - lat);
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.PI/2 - Math.atan(Math.hypot(x, y)), Math.atan2(x, -y) };
		}
	};
	
	
	public static final Projection ORTHOGRAPHIC = new Projection(
			"Orthographic", "A projection that mimics the Earth viewed from a great distance.",
			new BoundingBox(2, 2), 0b0111, Type.AZIMUTHAL, Property.PERSPECTIVE, 3) {
		
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
	
	
	public static final Projection PERSPECTIVE = new Projection(
			"Perspective", "A projection that mimics the actual appearance of the Earth.",
			null, 0b0111, Type.AZIMUTHAL, Property.PERSPECTIVE, 4,
			new String[] {"Percentage"}, new double[][] {{1,99,33.3}}) {
		
		private double d; //viewing distance in sphere radii
		
		public void initialize(double... params) {
			this.d = 1/(1 - 2*params[0]/100);
			double radius = (Double.isFinite(d)) ? 1/Math.sqrt(d*d-1) : ORTHOGRAPHIC.bounds.xMax;
			this.bounds = new BoundingBox(-radius, radius, -radius, radius);
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
			if (h > bounds.xMax) 	return null;
			final double theta = Math.atan(h);
			final double phi = Math.acos(d*Math.sin(theta)) + theta;
			if (phi < Math.PI/2)
				return new double[] { phi, Math.atan2(x, -y) };
			else
				return new double[] { Math.PI - phi, Math.atan2(x, -y) };
		}
	};
	
	
	public static final Projection MAGNIFIER = new Projection(
			"Magnifying glass", "A projection that dilates its center to great scales.",
			new BoundingBox(2, 2), 0b1111, Type.AZIMUTHAL, Property.POINTLESS, 2,
			new String[] {"Actual size", "Apparent size"},
			new double[][] {{1, 60, 20}, {0.5, 1.0, 0.5}}) {
		
		private double p0, r0; // scale factor of magnified portion
		
		public void initialize(double... params) {
			this.p0 = Math.toRadians(params[0]);
			this.r0 = params[1];
		}
		
		public double[] project(double lat, double lon) {
			double p = Math.PI/2 - lat;
			double r;
			if (p < p0)
				r = r0*Math.sin(p/2)/Math.sin(p0/2);
			else
				r = Math.sqrt(1 - (1 - r0*r0)*Math.pow(Math.cos(p/2)/Math.cos(p0/2), 2));
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			double r = Math.hypot(x, y);
			double th = Math.atan2(x, -y);
			double p;
			if (r <= r0)
				p = 2*Math.asin(Math.sin(p0/2)*r/r0);
			else if (r <= 1)
				p = 2*Math.acos(Math.cos(p0/2)*Math.sqrt((1 - r*r)/(1 - r0*r0)));
			else
				return null;
			return new double[] { Math.PI/2 - p, th};
		}
	};
}
