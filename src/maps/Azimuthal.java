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

import static java.lang.Double.isInfinite;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.hypot;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.toRadians;

public class Azimuthal {
	
	public static final Projection STEREOGRAPHIC = new Projection(
			"Stereographic", "A mathematically important conformal azimuthal projection",
			Shape.rectangle(4, 4), true, false, true, true, Type.AZIMUTHAL, Property.CONFORMAL, 2) {
		
		public double[] project(double lat, double lon) {
			if (lat < -1.5) 	lat = -1.5;
			final double r = 1/(tan(lat/2 + PI/4));
			return new double[] {r*sin(lon), -r*cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { PI/2 - 2*atan(hypot(x, y)), atan2(x, -y) };
		}
	};
	
	
	public static final Projection POLAR = new Projection(
			"Azimuthal Equidistant", "An equidistant azimuthal projection",
			Shape.circle(PI), true, true, true, true, Type.AZIMUTHAL,
			Property.EQUIDISTANT, 2) {
		
		public double[] project(double lat, double lon) {
			final double r = PI/2 - lat;
			return new double[] {r*sin(lon), -r*cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			double phi = PI/2 - hypot(x, y);
			if (phi > -PI/2)
				return new double[] {phi, atan2(x, -y)};
			else
				return null;
		}
	};
	
	
	public static final Projection EQUAL_AREA = new Projection(
			"Azimuthal Equal-Area", "An equal-area azimuthal projection",
			Shape.circle(1), true, true, true, true, Type.AZIMUTHAL, Property.EQUAL_AREA, 1) {

		public double[] project(double lat, double lon) {
			final double r = cos((PI/2+lat)/2);
			return new double[] {r*sin(lon), -r*cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			double r = hypot(x, y);
			if (r <= 1)
				return new double[] { asin(1-2*r*r), atan2(x, -y) };
			else
				return null;
		}
	};
	
	
	public static final Projection GNOMONIC = new Projection(
			"Gnomonic", "A projection that draws all great circles as straight lines",
			Shape.rectangle(4, 4), true, false, true, true, Type.AZIMUTHAL, Property.GNOMONIC, 2) {

		public double[] project(double lat, double lon) {
			if (lat < 0.2) 	lat = 0.2;
			final double r = tan(PI/2 - lat);
			return new double[] { r*sin(lon), -r*cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { PI/2 - atan(hypot(x, y)), atan2(x, -y) };
		}
	};
	
	
	public static final Projection ORTHOGRAPHIC = new Projection(
			"Orthographic", "A projection that mimics the Earth viewed from a great distance",
			Shape.circle(1), true, false, true, true, Type.AZIMUTHAL, Property.PERSPECTIVE, 3) {

		public double[] project(double lat, double lon) {
			if (lat < 0)	lat = 0;
			return new double[] { cos(lat)*sin(lon), -cos(lat)*cos(lon) };
		}
		
		public double[] inverse(double x, double y) {
			double R = hypot(x, y);
			if (R <= 1)
				return new double[] { acos(R), atan2(x, -y) };
			else
				return null;
		}
	};
	
	
	public static final Projection PERSPECTIVE = new Projection(
			"Perspective", "A projection that mimics the actual appearance of the Earth",
			null, true, false, true, true, Type.AZIMUTHAL, Property.PERSPECTIVE, 4,
			new String[] {"Percentage"}, new double[][] {{1,99,33.3}}) {
		
		private double d; //viewing distance in sphere radii
		
		public void initialize(double... params) {
			this.d = 1/(1 - 2*params[0]/100);
			this.shape = Shape.circle((Double.isFinite(d)) ? 1/Math.sqrt(d*d - 1) : 1);
		}
		
		public double[] project(double lat, double lon) {
			if (isInfinite(d)) 	return ORTHOGRAPHIC.project(lat, lon);
			if (lat < asin(1/d)) 	lat = asin(1/d);
			final double r = abs(cos(lat)/(d - sin(lat)));
			return new double[] { r*sin(lon), -r*cos(lon) };
		}
		
		public double[] inverse(double x, double y) {
			if (Double.isInfinite(d)) 	return ORTHOGRAPHIC.inverse(x, y);
			final double h = Math.hypot(x, y);
			if (h > this.shape.xMax) 	return null;
			final double theta = Math.atan(h);
			final double phi = Math.acos(d*Math.sin(theta)) + theta;
			if (phi < Math.PI/2)
				return new double[] { phi, Math.atan2(x, -y) };
			else
				return new double[] { PI - phi, atan2(x, -y) };
		}
	};
	
	
	public static final Projection MAGNIFIER = new Projection(
			"Magnifying glass", "A projection that dilates its center to great scales",
			Shape.circle(1), true, true, true, true, Type.AZIMUTHAL, Property.POINTLESS, 2,
			new String[] {"Actual size", "Apparent size"},
			new double[][] {{1, 60, 20}, {0.5, 1.0, 0.5}}) {
		
		private double p0, r0; // scale factor of magnified portion
		
		public void initialize(double... params) {
			this.p0 = toRadians(params[0]);
			this.r0 = params[1];
		}
		
		public double[] project(double lat, double lon) {
			double p = PI/2 - lat;
			double r;
			if (p < p0)
				r = r0*sin(p/2)/sin(p0/2);
			else
				r = sqrt(1 - (1 - r0*r0)*pow(cos(p/2)/cos(p0/2), 2));
			return new double[] { r*sin(lon), -r*cos(lon)};
		}
		
		public double[] inverse(double x, double y) {
			double r = hypot(x, y);
			double th = atan2(x, -y);
			double p;
			if (r <= r0)
				p = 2*asin(sin(p0/2)*r/r0);
			else if (r <= 1)
				p = 2*acos(cos(p0/2)*sqrt((1 - r*r)/(1 - r0*r0)));
			else
				return null;
			return new double[] { PI/2 - p, th};
		}
	};
}
