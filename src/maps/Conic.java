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
import utils.Shape;

import static java.lang.Double.isNaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.asin;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.hypot;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;
import static utils.Math2.max;
import static utils.Math2.min;

public class Conic {
	
	public static final Projection LAMBERT =
			new ConicProjection("Conformal Conic", false, true, true, Property.CONFORMAL, 2) {
		
		private double n; //the scaling factor for angles

		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Mercator; indicate with n=0
				this.n = 0;
			else if (lat1 == lat2) //equation becomes indeterminate; use limit
				this.n = sin(lat1);
			else //normal conic
				this.n = log(cos(lat1)/cos(lat2))/log(tan(PI/4+lat2/2)/tan(PI/4+lat1/2));

			// calculate the radii of the standard parallels
			final double r1 = pow(tan(PI/4 + max(lat1, lat2, 0)/2), -n);
			final double r2 = pow(tan(PI/4 + min(lat1, lat2, 0)/2), -n);
			// and expand those to get a minor and major radius
			final double r = max(2*r1 - r2, 0);
			final double R = 2*r2 - r1;
			// if the angle is greater than 180°
			if (n > 0.5) {
				// the upper bound is set by the outer radius
				this.shape = Shape.polygon(new double[][] {
						{0, 0},
						{-R*sin(PI*n), -R*cos(PI*n)},
						{-R, -R*cos(PI*n)},
						{-R, -R},
						{R, -R},
						{R, -R*cos(PI*n)},
						{R*sin(PI*n), -R*cos(PI*n)},
				});
			}
			// if the angle is less than 180°
			else if (n > 0) {
				// the upper bound is set by the inner or central radius
				double yCenter = -(R+r)/2*cos(PI*n);
				double width = 2*R*sin(PI*n);
				double yMin = min(-R, yCenter - width/2);
				double yMax = min(max(-r*cos(PI*n), yCenter + width/2), 0);
				this.shape = Shape.polygon(new double[][] {
						{yMax*tan(PI*n), yMax},
						{-width/2, -R*cos(PI*n)},
						{-width/2, yMin},
						{width/2, yMin},
						{width/2, -R*cos(PI*n)},
						{-yMax*tan(PI*n), yMax},
				});
			}
			// if the angle is 0°
			else {
				// this changes to a Mercator projection
				double width = 2*PI;
				double height = max(6*log(tan(PI/4+abs(lat1)/2)), width);
				this.shape = Shape.rectangle(width, height);
			}

			// reverse the shape if this is inverted
			if (reversed)
				this.shape = Shape.scaled(shape, 1, -1);
		}
		
		public double[] project(double lat, double lon) {
			if (n == 0) 	return Cylindrical.MERCATOR.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			if (lat < -1.5) 	lat = -1.5; //remove polar infinite values
			final double s = reversed ? -1 : 1;
			final double r = pow(tan(PI/4+lat/2), -n);
			return new double[] { s*r*sin(n*lon), -s*r*cos(n*lon) };
		}
		
		public double[] inverse(double x, double y) {
			if (n == 0)
				return Cylindrical.MERCATOR.inverse(x, y);
			else if (reversed) {
				x = -x;
				y = -y;
			}
			final double r = hypot(x, -y);
			final double phi = 2*atan(pow(r, -1/n)) - PI/2;
			final double lam = atan2(x, -y)/n;
			if (abs(lam) > PI)
				return null;
			else if (reversed)
				return new double[] {-phi, -lam};
			else
				return new double[] {phi, lam};
		}
	};
	
	
	public static final Projection EQUIDISTANT =
			new ConicProjection("Equidistant Conic", true, true, true, Property.EQUIDISTANT, 2) {
		
		private double m; //the scaling factor for radii
		private double n; //the scaling factor for angles

		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Equirectangular; indicate with m=0
				this.m = 0;
			else if (lat1 == lat2) //equation becomes indeterminate; use limit
				this.m = 1/(1/tan(lat1)/PI + lat1/PI + .5);
			else //normal conic
				this.m = (1/cos(lat2)-1/cos(lat1))/((-lat1/PI-.5)/cos(lat1)-(-lat2/PI-.5)/cos(lat2));
			
			this.n = m*cos(lat1)/((-lat1/PI-.5)*m+1)/PI;

			// if this is well-formed, use the annular sector shape
			if (n > 0) {
				this.shape = Shape.annularSector(1 - m, 1, 2*PI*n, reversed);
			}
			// if the angle is 0°, this changes to an equirectangular projection
			else {
				Cylindrical.EQUIRECTANGULAR.initialize(toDegrees(lat1));
				this.shape = Cylindrical.EQUIRECTANGULAR.shape;
			}
		}
		
		public double[] project(double lat, double lon) {
			if (m == 0) 	return Cylindrical.EQUIRECTANGULAR.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			final double s = reversed ? -1 : 1;
			final double r = 1 - m*lat/PI - m/2;
			return new double[] { s*r*sin(n*lon), s*(-r*cos(n*lon)) };
		}
		
		public double[] inverse(double x, double y) {
			if (m == 0) 	return Cylindrical.EQUIRECTANGULAR.inverse(x, y);
			if (reversed) {
				x = -x;
				y = -y;
			}
			final double r = hypot(x, y);
			final double phi = (1 - m/2 - r)*PI/m;
			final double lam = atan2(x, -y)/n;
			if (abs(lam) > PI || abs(phi) > PI/2)
				return null;
			else if (reversed) 	return new double[] {-phi, -lam};
			else 				return new double[] {phi, lam};
		}
	};
	
	
	public static final Projection ALBERS =
			new ConicProjection("Albers", true, true, true, Property.EQUAL_AREA, 2) {
		
		private double n; //the scaling factor for angles
		private double C; //a scaling factor for radii

		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Equirectangular; indicate with n=0
				this.n = 0;
			else //normal conic
				this.n = (sin(lat1) + sin(lat2))/2;
			
			this.C = pow(cos(lat1), 2) + 2*n*sin(lat1);

			// if this is well-formed, use the annular sector shape
			if (n > 0) {
				this.shape = Shape.annularSector(
						sqrt(C - 2*n), sqrt(C + 2*n), 2*PI*n, reversed);
			}
			// if the angle is zero, this becomes a cylindrical equal area projection
			else {
				Cylindrical.EQUAL_AREA.initialize(toDegrees(lat1));
				this.shape = Cylindrical.EQUAL_AREA.shape;
			}
		}
		
		public double[] project(double lat, double lon) {
			if (n == 0) 	return Cylindrical.EQUAL_AREA.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			final double r = sqrt(C - 2*n*sin(lat));
			final double x = r*sin(n*lon);
			final double y = -r*cos(n*lon);
			if (reversed) 	return new double[] {-x,-y};
			else 			return new double[] { x, y};
		}
		
		public double[] inverse(double x, double y) {
			if (n == 0) 	return Cylindrical.EQUAL_AREA.inverse(x, y);
			if (reversed) {
				x = -x;
				y = -y;
			}
			final double r = hypot(x, y);
			final double phi = asin((C - pow(r,2))/(2*n));
			final double lam = atan2(x, -y)/n;
			if (abs(lam) > PI || isNaN(phi))
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
		
		ConicProjection(String name, boolean comprehensive, boolean solvable, boolean invertible,
		                Property property, int rating) {
			super(name, "The " + property + " conic projection", null, true, comprehensive, solvable,
			      invertible, Type.CONIC, property, rating,
			      new String[] {"Std. Parallel 1", "Std. Parallel 2"},
			      new double[][] {{-89,89,15},{-89,89,45}});
		}
		
		public final void initialize(double... params) {
			this.lat1 = toRadians(params[0]);
			this.lat2 = toRadians(params[1]);
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
