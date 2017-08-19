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

public class Conic {
	
	public static final Projection LAMBERT =
			new ConicProjection("Conformal Conic", 0b0111, Property.CONFORMAL) {
		
		private double n, F, r0, d;
		
		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Mercator; indicate with n=0
				this.n = 0;
			else if (lat1 == lat2) //equation becomes indeterminate; use limit
				this.n = Math.sin(lat1);
			else //normal conic
				this.n = Math.log(Math.cos(lat1)/Math.cos(lat2))/Math.log(Math.tan(Math.PI/4+lat2/2)/Math.tan(Math.PI/4+lat1/2));
			
			this.F = Math.cos(lat1)*Math.pow(Math.tan(Math.PI/4+lat1/2), n)/n/Math.PI;
			this.r0 = F*Math.pow(Math.tan(Math.PI/4+(lat1+lat2)/4), -n);
			
			if (n >= .5 && -Math.tan(Math.PI*n)*(1-r0) < 1)
				this.d = 0;
			else if (n >= .5)
				this.d = 1 - r0 + 1/Math.tan(Math.PI*n);
			else if (Math.tan(Math.PI*n)*(1+r0) < 1)
				this.d = (1 - 1/Math.tan(Math.PI*n))/r0;
			else if (r0 >= 1)
				this.d = 0;
			else
				this.d = 1 - r0; //TODO I can do better
			this.aspectRatio = 2/(2 - Math.max(d, 0));
		}
		
		public double[] project(double lat, double lon) {
			if (n == 0) 	return Cylindrical.MERCATOR.project(lat, lon);
			if (reversed) {
				lat = -lat; lon = -lon;
			}
			final double r = F*Math.pow(Math.tan(Math.PI/4+lat/2), -n);
			final double s = reversed ? -1 : 1;
			final double x = s*Math.PI*r*Math.sin(n*lon);
			final double y = s*Math.PI*(r0 - r*Math.cos(n*lon));
			if (d > 0) 			return new double[] {x, y+Math.PI*d/2};
			else if (d < 0) 	return new double[] {-d*x, -d*y};
			else 				return new double[] {x, y};
		}
		
		public double[] inverse(double x, double y) {
			if (n == 0) 	return Cylindrical.MERCATOR.inverse(x, y);
			else if (reversed) {
				x = -x; y = -y;
			}
			if (d > 0)
				y = y*(1-d/2) - d/2;
			else if (d < 0) {
				x /= -d; y /= -d;
			}
			final double r = Math.hypot(x, r0-y);
			final double phi = 2*Math.atan(Math.pow(F/r, 1/n)) - Math.PI/2;
			final double lam = Math.atan2(x, r0-y)/n;
			if (Math.abs(lam) > Math.PI) 	return null;
			else if (reversed) 				return new double[] {-phi, -lam};
			else 							return new double[] {phi, lam};
		}
	};
	
	
	public static final Projection EQUIDISTANT =
			new ConicProjection("Equidistant Conic", 0b1111, Property.EQUIDISTANT) {
		
		private double m, n, s;
		
		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Equirectangular; indicate with m=0
				this.m = 0;
			else if (lat1 == lat2) //equation becomes indeterminate; use limit
				this.m = 1/(1/Math.tan(lat1)/Math.PI + lat1/Math.PI + .5);
			else //normal conic
				this.m = (1/Math.cos(lat2)-1/Math.cos(lat1))/((-lat1/Math.PI-.5)/Math.cos(lat1)-(-lat2/Math.PI-.5)/Math.cos(lat2));
			
			this.n = m*Math.cos(lat1)/((-lat1/Math.PI-.5)*m+1)/Math.PI;
			
			if (m == 0) {
				this.aspectRatio = 2*Math.cos(lat1);
			}
			else if (n >= .5) {
				this.aspectRatio = 2/(1 - Math.cos(Math.PI*n));
				this.s = 1;
			}
			else {
				this.aspectRatio = 2*Math.sin(Math.PI*n)/(1 - (1-m)*Math.cos(Math.PI*n));
				if (aspectRatio >= 1) 	this.s = 1/Math.sin(Math.PI*n);
				else 					this.s = 2/(1 - (1-m)*Math.cos(Math.PI*n));
			}
		}
		
		public double[] project(double lat, double lon) {
			if (m == 0) 	return Cylindrical.EQUIRECTANGULAR.project(lat, lon);
			if (reversed) {
				lat = -lat; lon = -lon;
			}
			final double r = Math.PI - m*lat - Math.PI*m/2;
			final double x = s*r*Math.sin(n*lon);
			final double y = Math.PI*(s-1/Math.max(aspectRatio,1)) - s*r*Math.cos(n*lon);
			if (reversed) 	return new double[] { -x, -y };
			else 			return new double[] { x, y };
		}
		
		public double[] inverse(double x, double y) {
			if (m == 0) 	return Cylindrical.EQUIRECTANGULAR.inverse(x, y);
			if (reversed) {
				x = -x; y = -y;
			}
			x = x/s;
			y = (y+1)/s/Math.max(aspectRatio, 1) - 1;
			final double r = Math.hypot(x, y);
			final double phi = (1 - m/2 - r)*Math.PI/m;
			final double lam = Math.atan2(x, -y)/n;
			if (Math.abs(lam) > Math.PI || Math.abs(phi) > Math.PI/2)
				return null;
			else if (reversed) 				return new double[] {-phi, -lam};
			else 							return new double[] {phi, lam};
		}
	};
	
	
	public static final Projection ALBERS =
			new ConicProjection("Albers", 0b1111, Property.EQUAL_AREA) {
		
		private double n, C, s;
		
		public void setSpecificParameters() {
			if (lat1 == -lat2) //degenerates into Equirectangular; indicate with n=0
				this.n = 0;
			else //normal conic
				this.n = (Math.sin(lat1) + Math.sin(lat2))/2;
			
			this.C = Math.pow(Math.cos(lat1), 2) + 2*n*Math.sin(lat1);
			if (n == 0) {
				this.aspectRatio = Math.PI*Math.pow(Math.cos(lat1), 2);
			}
			else if (n >= .5) {
				this.aspectRatio = 2/(1 - Math.cos(Math.PI*n));
				this.s = n/Math.sqrt(C+2*n);
			}
			else {
				this.aspectRatio = 2*Math.sin(Math.PI*n)/(1 - Math.sqrt(C-2*n)/Math.sqrt(C+2*n)*Math.cos(Math.PI*n));
				if (aspectRatio >= 1)
					this.s = n/Math.sqrt(C+2*n)/Math.sin(Math.PI*n);
				else
					this.s = 2/(Math.sqrt(C+2*n)/n - Math.sqrt(C-2*n)/n*Math.cos(Math.PI*n));
			}
		}
		
		public double[] project(double lat, double lon) {
			if (n == 0) 	return Cylindrical.EQUAL_AREA.project(lat, lon);
			if (reversed) {
				lat = -lat; lon = -lon;
			}
			final double r = Math.sqrt(C - 2*n*Math.sin(lat))/n;
			final double x = Math.PI*s*r*Math.sin(n*lon);
			final double y = Math.PI*(s*Math.sqrt(C+2*n)/n-1/Math.max(aspectRatio,1)) - Math.PI*s*r*Math.cos(n*lon);
			if (reversed) 	return new double[] { -x, -y };
			else 			return new double[] { x, y };
		}
		
		public double[] inverse(double x, double y) {
			if (n == 0) 	return Cylindrical.EQUAL_AREA.inverse(x, y);
			if (reversed) {
				x = -x; y = -y;
			}
			x = x/s;
			y = (y+1)/s/Math.max(aspectRatio, 1) - Math.sqrt(C+2*n)/n;
			final double r = Math.hypot(x, y);
			final double phi = Math.asin((C - Math.pow(n*r,2))/(2*n));
			final double lam = Math.atan2(x, -y)/n;
			if (Math.abs(lam) > Math.PI || Double.isNaN(phi))
				return null;
			else if (reversed) 				return new double[] {-phi, -lam};
			else 							return new double[] {phi, lam};
		}
	};
	
	
	
	private static abstract class ConicProjection extends Projection {
		
		protected double lat1, lat2;
		protected boolean reversed;
		
		ConicProjection(String name, int fisc, Property property) {
			super(name, "The "+property+" conic projection.", 0, fisc, Type.CONIC, property,
					new String[] {"Std. Parallel 1", "Std. Parallel 2"},
					new double[][] {{-89,89,15},{-89,89,45}});
		}
		
		public final void setParameters(double... params) {
			this.lat1 = Math.toRadians(params[0]);
			this.lat2 = Math.toRadians(params[1]);
			this.reversed = lat1 + lat2 < 0;
			if (reversed) {
				lat1 = -lat1; lat2 = -lat2;
			}
			setSpecificParameters();
		}
		
		protected abstract void setSpecificParameters(); //a way to require subclasses to set lat1 and lat2
	}
}
