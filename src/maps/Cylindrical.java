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
	
	public static final Projection MERCATOR = new Projection(
			"Mercator", 2*Math.PI, 2*Math.PI, 0b0111, Type.CYLINDRICAL, Property.CONFORMAL, 1,
			"very popular") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, Math.log(Math.tan(Math.PI/4+lat/2))};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {Math.atan(Math.sinh(y)), x};
		}
	};
	
	
	public static final Projection PLATE_CARREE = new Projection(
			"Plate Carr\u00E9e", 2*Math.PI, Math.PI, 0b1111, Type.CYLINDRICAL,
			Property.EQUIDISTANT, 2, null, "focused on the equator"){
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, lat};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {y, x};
		}
	};
	
	
	public static final Projection EQUIRECTANGULAR = new Projection(
			"Equirectangular", "A linear mapping from longitude and latitude to x and y.",
			2*Math.PI, 0., 0b1111, Type.CYLINDRICAL, Property.EQUIDISTANT, 2,
			new String[]{"Std. parallel"}, new double[][]{{0, 89, 0}}) {
		
		private double stdParallel;
		
		public void setParameters(double... params) {
			this.stdParallel = Math.toRadians(params[0]);
			this.height = Math.PI/Math.cos(stdParallel);
		}
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, lat/Math.cos(stdParallel)};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {y*Math.cos(stdParallel), x};
		}
	};
	
	
	public static final Projection GALL_ORTHOGRAPHIC = new Projection(
			"Gall-Peters", 2*Math.PI, 4, 0b1111, Type.CYLINDRICAL, Property.EQUAL_AREA, 0,
			"somewhat controversial", "with least distortion at 45\u00B0") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, Math.sin(lat)*height/2};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y*2/height), x};
		}
	};
	
	
	public static final Projection HOBO_DYER = new Projection(
			"Hobo-Dyer", 2*Math.PI, 3.178, 0b1111, Type.CYLINDRICAL, Property.EQUAL_AREA, 2,
			null, "with least distortion at 37.5\u00B0") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, Math.sin(lat)*height/2};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y*2/height), x };
		}
	};
	
	
	public static final Projection BEHRMANN = new Projection(
			"Behrmann", 2*Math.PI, 3, 0b1111, Type.CYLINDRICAL, Property.EQUAL_AREA, 3,
			null, "with least distortion at 30\u00B0") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, Math.sin(lat)*height/2};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y*2/height), x };
		}
	};
	
	
	public static final Projection LAMBERT = new Projection(
			"Lambert cylindrical", 2*Math.PI, 2, 0b1111, Type.CYLINDRICAL, Property.EQUAL_AREA, 2,
			null, "with least distortion along the equator") {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, Math.sin(lat)*height/2};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y*2/height), x };
		}
	};
	
	
	public static final Projection EQUAL_AREA = new Projection(
			"Cylindrical Equal-area", "A generalized equal-area cylindrical projection.",
			2*Math.PI, 0, 0b1111, Type.CYLINDRICAL, Property.EQUAL_AREA, 2,
			new String[]{"Std. parallel"}, new double[][]{{0, 89, 30}}) {
		
		public void setParameters(double... params) {
			this.height = 2/Math.pow(Math.cos(Math.toRadians(params[0])), 2);
		}
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, Math.sin(lat)*height/2};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { Math.asin(y*2/height), x };
		}
	};
	
	
	public static final Projection GALL_STEREOGRAPHIC = new Projection(
			"Gall Stereographic", 2*Math.PI, 1.5*Math.PI, 0b1111, Type.CYLINDRICAL,
			Property.COMPROMISE, 2) {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, Math.tan(lat/2)*(1+Math.sqrt(2))};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { 2*Math.atan(y/(1+Math.sqrt(2))), x };
		}
	};
	
	
	public static final Projection MILLER = new Projection(
			"Miller", 2*Math.PI, 2.5*Math.log(Math.tan(9*Math.PI/20)), 0b1111, Type.CYLINDRICAL,
			Property.COMPROMISE, 2) {
		
		public double[] project(double lat, double lon) {
			return new double[] {lon, Math.log(Math.tan(Math.PI/4+.8*lat/2))/.8};
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {Math.atan(Math.sinh(y*.8))/.8, x};
		}
	};
}
