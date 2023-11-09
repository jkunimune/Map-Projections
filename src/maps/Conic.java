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
import utils.BoundingBox;
import utils.Math2;

public class Conic {
	
	public static final Projection LAMBERT =
			new ConicProjection("Conformal Conic", 0b0111, Property.CONFORMAL, 2) {
		
		private double n; //the scaling factor for angles

		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Mercator; indicate with n=0
				this.n = 0;
			else if (lat1 == lat2) //equation becomes indeterminate; use limit
				this.n = Math.sin(lat1);
			else //normal conic
				this.n = Math.log(Math.cos(lat1)/Math.cos(lat2))/Math.log(Math.tan(Math.PI/4+lat2/2)/Math.tan(Math.PI/4+lat1/2));

			// calculate the radii of the standard parallels
			final double r1 = Math.pow(Math.tan(Math.PI/4+Math2.max(lat1,lat2,0)/2), -n);
			final double r2 = Math.pow(Math.tan(Math.PI/4+Math2.min(lat1,lat2,0)/2), -n);
			// and expand those to get a minor and major radius
			final double r = Math.max(2*r1 - r2, 0);
			final double R = 2*r2 - r1;
			// if the angle is greater than 180°
			if (n > 0.5) {
				// the upper bound is set by the outer radius
				this.bounds = new BoundingBox(-R, R, -R, -R*Math.cos(Math.PI*n));
			}
			// if the angle is less than 180°
			else if (n > 0) {
				// the upper bound is set by the inner radius
				double yCenter = -(R+r)/2*Math.cos(Math.PI*n);
				double xMax = R*Math.sin(Math.PI*n);
				double yMax = -r*Math.cos(Math.PI*n);
				this.bounds = new BoundingBox(
						-xMax, xMax,
						Math.min(-R, yCenter - xMax),
						Math.min(Math.max(yMax, yCenter + xMax), 0));
			}
			// if the angle is 0°
			else {
				// this changes to a Mercator projection
				double width = 2*Math.PI;
				double height = Math.max(6*Math.log(Math.tan(Math.PI/4+Math.abs(lat1)/2)), width);
				this.bounds = new BoundingBox(width, height);
			}

			// reverse the bounds if this is inverted
			if (reversed)
				this.bounds = new BoundingBox(
						-bounds.xMax, -bounds.xMin,
						-bounds.yMax, -bounds.yMin);
		}
		
		public double[] project(double lat, double lon) {
			if (n == 0) 	return Cylindrical.MERCATOR.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			if (lat < -1.5) 	lat = -1.5; //remove polar infinite values
			final double s = reversed ? -1 : 1;
			final double r = Math.pow(Math.tan(Math.PI/4+lat/2), -n);
			return new double[] { s*r*Math.sin(n*lon), -s*r*Math.cos(n*lon) };
		}
		
		public double[] inverse(double x, double y) {
			if (n == 0) 	return Cylindrical.MERCATOR.inverse(x, y);
			else if (reversed) {
				x = -x;
				y = -y;
			}
			final double r = Math.hypot(x, -y);
			final double phi = 2*Math.atan(Math.pow(r, -1/n)) - Math.PI/2;
			final double lam = Math.atan2(x, -y)/n;
			if (Math.abs(lam) > Math.PI) 	return null;
			else if (reversed) 				return new double[] {-phi, -lam};
			else 							return new double[] {phi, lam};
		}
	};
	
	
	public static final Projection EQUIDISTANT =
			new ConicProjection("Equidistant Conic", 0b1111, Property.EQUIDISTANT, 2) {
		
		private double m; //the scaling factor for radii
		private double n; //the scaling factor for angles

		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Equirectangular; indicate with m=0
				this.m = 0;
			else if (lat1 == lat2) //equation becomes indeterminate; use limit
				this.m = 1/(1/Math.tan(lat1)/Math.PI + lat1/Math.PI + .5);
			else //normal conic
				this.m = (1/Math.cos(lat2)-1/Math.cos(lat1))/((-lat1/Math.PI-.5)/Math.cos(lat1)-(-lat2/Math.PI-.5)/Math.cos(lat2));
			
			this.n = m*Math.cos(lat1)/((-lat1/Math.PI-.5)*m+1)/Math.PI;

			// if the angle is greater than 180°
			if (n > 0.5) {
				// the upper bound is set by the outer radius
				this.bounds = new BoundingBox(-1, 1, -1, -Math.cos(Math.PI*n));
			}
			// if the angle is less than 180°
			else if (n > 0) {
				// the upper bound is set by the inner radius
				this.bounds = new BoundingBox(
						-Math.sin(Math.PI*n), Math.sin(Math.PI*n),
						-1, -(1-m)*Math.cos(Math.PI*n));
			}
			// if the angle is 0°
			else {
				// this changes to an equirectangular projection
				Cylindrical.EQUIRECTANGULAR.initialize(Math.toDegrees(lat1));
				this.bounds = Cylindrical.EQUIRECTANGULAR.bounds;
			}

			// reverse the bounds if this is inverted
			if (reversed)
				this.bounds = new BoundingBox(
						-bounds.xMax, -bounds.xMin,
						-bounds.yMax, -bounds.yMin);
		}
		
		public double[] project(double lat, double lon) {
			if (m == 0) 	return Cylindrical.EQUIRECTANGULAR.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			final double s = reversed ? -1 : 1;
			final double r = 1 - m*lat/Math.PI - m/2;
			return new double[] { s*r*Math.sin(n*lon), s*(-r*Math.cos(n*lon)) };
		}
		
		public double[] inverse(double x, double y) {
			if (m == 0) 	return Cylindrical.EQUIRECTANGULAR.inverse(x, y);
			if (reversed) {
				x = -x;
				y = -y;
			}
			final double r = Math.hypot(x, y);
			final double phi = (1 - m/2 - r)*Math.PI/m;
			final double lam = Math.atan2(x, -y)/n;
			if (Math.abs(lam) > Math.PI || Math.abs(phi) > Math.PI/2)
				return null;
			else if (reversed) 	return new double[] {-phi, -lam};
			else 				return new double[] {phi, lam};
		}
	};
	
	
	public static final Projection ALBERS =
			new ConicProjection("Albers", 0b1111, Property.EQUAL_AREA, 2) {
		
		private double n; //the scaling factor for angles
		private double C; //a scaling factor for radii

		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Equirectangular; indicate with n=0
				this.n = 0;
			else //normal conic
				this.n = (Math.sin(lat1) + Math.sin(lat2))/2;
			
			this.C = Math.pow(Math.cos(lat1), 2) + 2*n*Math.sin(lat1);
			
			final double r = Math.sqrt(C - 2*n);
			final double R = Math.sqrt(C + 2*n);
			// if the angle is greater than 180°
			if (n > 0.5) {
				// the upper bound is set by the outer radius
				this.bounds = new BoundingBox(
						-R, R,
						-R, -R*Math.cos(Math.PI*n));
			}
			// if the angle is less than 180°
			else if (n > 0) {
				// the upper bound is set by the inner radius
				this.bounds = new BoundingBox(
						-R*Math.sin(Math.PI*n), R*Math.sin(Math.PI*n),
						-R, -r*Math.cos(Math.PI*n));
			}
			// if the angle is zero
			else {
				// this becomes a cylindrical equal area projection
				Cylindrical.EQUAL_AREA.initialize(Math.toDegrees(lat1));
				this.bounds = Cylindrical.EQUAL_AREA.bounds;
			}

			// reverse the bounds if this is inverted
			if (reversed)
				this.bounds = new BoundingBox(
						-bounds.xMax, -bounds.xMin,
						-bounds.yMax, -bounds.yMin);
		}
		
		public double[] project(double lat, double lon) {
			if (n == 0) 	return Cylindrical.EQUAL_AREA.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			final double r = Math.sqrt(C - 2*n*Math.sin(lat));
			final double x = r*Math.sin(n*lon);
			final double y = -r*Math.cos(n*lon);
			if (reversed) 	return new double[] {-x,-y};
			else 			return new double[] { x, y};
		}
		
		public double[] inverse(double x, double y) {
			if (n == 0) 	return Cylindrical.EQUAL_AREA.inverse(x, y);
			if (reversed) {
				x = -x;
				y = -y;
			}
			final double r = Math.hypot(x, y);
			final double phi = Math.asin((C - Math.pow(r,2))/(2*n));
			final double lam = Math.atan2(x, -y)/n;
			if (Math.abs(lam) > Math.PI || Double.isNaN(phi))
				return null;
			else if (reversed) 	return new double[] {-phi, -lam};
			else 				return new double[] {phi, lam};
		}
	};
	
	
	/**
	 * A base for all conic projections
	 * @author jkunimune
	 */
	private static abstract class ConicProjection extends Projection {
		
		protected double lat1, lat2;
		protected boolean reversed;
		
		ConicProjection(String name, int fisc, Property property, int rating) {
			super(name, "The "+property+" conic projection.", null, fisc, Type.CONIC, property,
					rating, new String[] {"Std. Parallel 1", "Std. Parallel 2"},
					new double[][] {{-89,89,15},{-89,89,45}});
		}
		
		public final void initialize(double... params) {
			this.lat1 = Math.toRadians(params[0]);
			this.lat2 = Math.toRadians(params[1]);
			this.reversed = lat1 + lat2 < 0;
			if (reversed) {
				lat1 = -lat1;
				lat2 = -lat2;
			}
			setSpecificParameters(); //this is where subclasses interpret lat1 and lat2
		}
		
		protected abstract void setSpecificParameters(); //a way to require subclasses to set lat1 and lat2
	}
}
