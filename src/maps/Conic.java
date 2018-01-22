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
import utils.Math2;

public class Conic {
	
	public static final Projection LAMBERT =
			new ConicProjection("Conformal Conic", 0b0111, Property.CONFORMAL) {
		
		private double n; //the scaling factor for angles
		private double r0; //the primary radius based on the standard axes
		
		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Mercator; indicate with n=0
				this.n = 0;
			else if (lat1 == lat2) //equation becomes indeterminate; use limit
				this.n = Math.sin(lat1);
			else //normal conic
				this.n = Math.log(Math.cos(lat1)/Math.cos(lat2))/Math.log(Math.tan(Math.PI/4+lat2/2)/Math.tan(Math.PI/4+lat1/2));
			
			final double r1 = Math.pow(Math.tan(Math.PI/4+Math2.max(lat1,lat2,0)/2), -n);
			final double r2 = Math.pow(Math.tan(Math.PI/4+Math2.min(lat1,lat2,0)/2), -n);
			final double r = Math.max(2*r1 - r2, 0);
			final double R = 2*r2 - r1;
			if (n > 0.5) {
				this.width = 2*R;
				this.height = R - R*Math.cos(Math.PI*n);
				this.r0 = (R + R*Math.cos(Math.PI*n))/2;
			}
			else if (n > 0) {
				this.width = 2*R*Math.sin(Math.PI*n);
				this.height = Math.max(R - r*Math.cos(Math.PI*n), width);
				this.r0 = (R+r)/2*Math.cos(Math.PI*n);
				if (height > 2*r0) {
					final double err = height/2-r0;
					this.height -= err;
					this.r0 += err/2;
				}
			}
			else {
				this.width = Cylindrical.MERCATOR.getWidth(); //TODO: fix this; make it consistent
				this.height = Cylindrical.MERCATOR.getHeight();
			}
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
			return new double[] { s*r*Math.sin(n*lon), s*(r0 - r*Math.cos(n*lon)) };
		}
		
		public double[] inverse(double x, double y) {
			if (n == 0) 	return Cylindrical.MERCATOR.inverse(x, y);
			else if (reversed) {
				x = -x;
				y = -y;
			}
			final double r = Math.hypot(x, r0-y);
			final double phi = 2*Math.atan(Math.pow(r, -1/n)) - Math.PI/2;
			final double lam = Math.atan2(x, r0-y)/n;
			if (Math.abs(lam) > Math.PI) 	return null;
			else if (reversed) 				return new double[] {-phi, -lam};
			else 							return new double[] {phi, lam};
		}
	};
	
	
	public static final Projection EQUIDISTANT =
			new ConicProjection("Equidistant Conic", 0b1111, Property.EQUIDISTANT) {
		
		private double m; //the scaling factor for radii
		private double n; //the scaling factor for angles
		private double y0; //the centered position
		
		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Equirectangular; indicate with m=0
				this.m = 0;
			else if (lat1 == lat2) //equation becomes indeterminate; use limit
				this.m = 1/(1/Math.tan(lat1)/Math.PI + lat1/Math.PI + .5);
			else //normal conic
				this.m = (1/Math.cos(lat2)-1/Math.cos(lat1))/((-lat1/Math.PI-.5)/Math.cos(lat1)-(-lat2/Math.PI-.5)/Math.cos(lat2));
			
			this.n = m*Math.cos(lat1)/((-lat1/Math.PI-.5)*m+1)/Math.PI;
			
			if (n > 0.5) {
				this.width = 2;
				this.height = 1 - Math.cos(Math.PI*n);
				this.y0 = (1 + Math.cos(Math.PI*n))/2;
			}
			else if (n > 0) {
				this.width = 2*Math.sin(Math.PI*n);
				this.height = 1 - (1-m)*Math.cos(Math.PI*n);
				this.y0 = (1 + (1-m)*Math.cos(Math.PI*n))/2;
			}
			else {
				Cylindrical.EQUIRECTANGULAR.setParameters(Math.toDegrees(lat1));
				this.width = Cylindrical.EQUIRECTANGULAR.getWidth();
				this.height = Cylindrical.EQUIRECTANGULAR.getHeight();
			}
		}
		
		public double[] project(double lat, double lon) {
			if (m == 0) 	return Cylindrical.EQUIRECTANGULAR.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			final double s = reversed ? -1 : 1;
			final double r = 1 - m*lat/Math.PI - m/2;
			return new double[] { s*r*Math.sin(n*lon), s*(y0 - r*Math.cos(n*lon)) };
		}
		
		public double[] inverse(double x, double y) {
			if (m == 0) 	return Cylindrical.EQUIRECTANGULAR.inverse(x, y);
			if (reversed) {
				x = -x;
				y = -y;
			}
			final double r = Math.hypot(x, y-y0);
			final double phi = (1 - m/2 - r)*Math.PI/m;
			final double lam = Math.atan2(x, y0-y)/n;
			if (Math.abs(lam) > Math.PI || Math.abs(phi) > Math.PI/2)
				return null;
			else if (reversed) 	return new double[] {-phi, -lam};
			else 				return new double[] {phi, lam};
		}
	};
	
	
	public static final Projection ALBERS =
			new ConicProjection("Albers", 0b1111, Property.EQUAL_AREA) {
		
		private double n; //the scaling factor for angles
		private double C; //a scaling factor for radii
		private double y0; //the centering y-shift
		
		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Equirectangular; indicate with n=0
				this.n = 0;
			else //normal conic
				this.n = (Math.sin(lat1) + Math.sin(lat2))/2;
			
			this.C = Math.pow(Math.cos(lat1), 2) + 2*n*Math.sin(lat1);
			
			final double r = Math.sqrt(C - 2*n);
			final double R = Math.sqrt(C + 2*n);
			if (n > 0.5) {
				this.width = 2*R;
				this.height = R - R*Math.cos(Math.PI*n);
				this.y0 = (R + R*Math.cos(Math.PI*n))/2;
			}
			else if (n > 0) {
				this.width = 2*R*Math.sin(Math.PI*n);
				this.height = R - r*Math.cos(Math.PI*n);
				this.y0 = (R + r*Math.cos(Math.PI*n))/2;
			}
			else {
				Cylindrical.EQUAL_AREA.setParameters(Math.toDegrees(lat1));
				this.width = Cylindrical.EQUAL_AREA.getWidth();
				this.height = Cylindrical.EQUAL_AREA.getHeight();
			}
		}
		
		public double[] project(double lat, double lon) {
			if (n == 0) 	return Cylindrical.EQUAL_AREA.project(lat, lon);
			if (reversed) {
				lat = -lat;
				lon = -lon;
			}
			final double s = reversed ? -1 : 1;
			final double r = Math.sqrt(C - 2*n*Math.sin(lat));
			return new double[] { s*r*Math.sin(n*lon), s*(y0 - r*Math.cos(n*lon)) };
		}
		
		public double[] inverse(double x, double y) {
			if (n == 0) 	return Cylindrical.EQUAL_AREA.inverse(x, y);
			if (reversed) {
				x = -x;
				y = -y;
			}
			final double r = Math.hypot(x, y-y0);
			final double phi = Math.asin((C - Math.pow(r,2))/(2*n));
			final double lam = Math.atan2(x, y0-y)/n;
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
		
		ConicProjection(String name, int fisc, Property property) {
			super(name, "The "+property+" conic projection.", 0,0, fisc, Type.CONIC, property,
					new String[] {"Std. Parallel 1", "Std. Parallel 2"},
					new double[][] {{-89,89,15},{-89,89,45}});
		}
		
		public final void setParameters(double... params) {
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
