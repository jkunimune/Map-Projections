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

/**
 * Map projections where x is a linear function of longitude.
 * 
 * @author jkunimune
 */
public class Cylindrical {
	
	public static final Projection MERCATOR =
			new Projection(
					"Mercator", 1., 0b0111, Type.CYLINDRICAL, Property.CONFORMAL, "very popular") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon/Math.PI, Math.log(Math.tan(Math.PI/4+lat/2))/Math.PI};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {Math.atan(Math.sinh(y*Math.PI)), x*Math.PI};
		}
	};
	
	
	public static final Projection PLATE_CARREE=
			new Projection("Plate Carr\u00E9e", 2., 0b1111,
					Type.CYLINDRICAL, Property.EQUIDISTANT, null, "focused on the equator"){
	
		public double[] project(double lat, double lon) {
			return new double[] {lon/Math.PI, lat/Math.PI};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {y*Math.PI/2, x*Math.PI};
		}
	};
	
	
	public static final Projection EQUIRECTANGULAR =
			new Projection("Equirectangular",
					"A linear mapping from longitude and latitude to x and y",
					2, 0b1111, Type.CYLINDRICAL, Property.EQUIDISTANT,
					new String[]{"Std. parallel"}, new double[][]{{0, 85, 0}}) {
		
		public void setParameters(double... params) {
			this.aspectRatio = 2*Math.cos(Math.toRadians(params[0]));
		}
		
		public double[] project(double lat, double lon) {
			if (aspectRatio >= 1)
				return new double[] {lon/Math.PI, 2*lat/Math.PI/aspectRatio};
			else
				return new double[] {lon/Math.PI*aspectRatio, 2*lat/Math.PI};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {y*Math.PI/2, x*Math.PI};
		}
	};
	
	
	public static final Projection GALL_PETERS =
			new Projection("Gall-Peters", Math.PI/2, 0b1111,
					Type.CYLINDRICAL, Property.EQUAL_AREA, "somewhat controversial") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon/Math.PI, Math.sin(lat)/aspectRatio};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
	};
	
	
	public static final Projection HOBO_DYER =
			new Projection("Hobo-Dyer", 1.977, 0b1111, Type.CYLINDRICAL, Property.EQUAL_AREA, null,
					"with least distortion at 37.5\u00B0") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon/Math.PI, Math.sin(lat)/aspectRatio};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
	};
	
	
	public static final Projection BEHRMANN =
			new Projection("Behrmann", 2.356, 0b1111, Type.CYLINDRICAL, Property.EQUAL_AREA, null,
					"with least distortion at 30\u00B0") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon/Math.PI, Math.sin(lat)/aspectRatio};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
	};
	
	
	public static final Projection LAMBERT =
			new Projection("Lambert cylindrical", Math.PI, 0b1111, Type.CYLINDRICAL,
					Property.EQUAL_AREA, null, "with least distortion along the equator") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon/Math.PI, Math.sin(lat)/aspectRatio};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
	};
	
	
	public static final Projection EQUAL_AREA =
			new Projection("Equal-area cylindrical", "A generalized equal-area cylindrical projection",
					0, 0b1111, Type.CYLINDRICAL, Property.EQUAL_AREA,
					new String[]{"Std. parallel"}, new double[][]{{0, 85, 30}}) {
		
		public void setParameters(double... params) {
			final double stdParallel = Math.toRadians(params[0]);
			this.aspectRatio = Math.PI*Math.pow(Math.cos(stdParallel), 2);
		}
		
		public double[] project(double lat, double lon) {
			if (aspectRatio >= 1)
				return new double[] {lon/Math.PI, Math.sin(lat)/aspectRatio};
			else
				return new double[] {lon*aspectRatio/Math.PI, Math.sin(lat)};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
	};
	
	
	public static final Projection GALL =
			new Projection(
					"Gall Stereographic", 4/3., 0b1111, Type.CYLINDRICAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon/Math.PI, Math.tan(lat/2)*(1+Math.sqrt(2))/Math.PI};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { 2*Math.atan(y), x*Math.PI };
		}
	};
	
	
	public static final Projection MILLER =
			new Projection("Miller", 4*Math.PI/5/Math.log(Math.tan(9*Math.PI/20)), 0b1111,
					Type.CYLINDRICAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon/Math.PI, 1.25*Math.log(Math.tan(Math.PI/4+.8*lat/2))/Math.PI};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {1.25*Math.atan(Math.sinh(y*Math.log(Math.tan(9*Math.PI/20)))), x*Math.PI};
		}
	};
}
