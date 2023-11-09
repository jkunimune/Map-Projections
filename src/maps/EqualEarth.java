/**
 * MIT License
 * 
 * Copyright (c) 2018 Justin Kunimune
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
import utils.NumericalAnalysis;

/**
 * A fairly simple projection, but I have it in its own file, because otherwise the way the height and width are
 * defined is inconvenient.
 * 
 * @author Justin Kunimune
 */
public class EqualEarth {

	private static final double[] A = {1.340264, -.081106, .000893, .003796};
	private static final double B = Math.sqrt(3)/2;
	private static final double Y_SCALE = poly9(Math.PI/3)/(Math.PI/3);
	
	
	public static final Projection EQUAL_EARTH = new Projection(
			"Equal Earth",
			new BoundingBox(2/B/poly8(0)*Math.PI, 2*poly9(Math.PI/3)),
			0b1111, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3, null,
			"specifically designed to woo Gall-Peters supporters away from that horrid projection") {
		
		public double[] project(double lat, double lon) {
			double th = Math.asin(B*Math.sin(lat));
			return new double[] { Math.cos(th)/B/poly8(th)*lon, poly9(th) };
		}
		
		public double[] inverse(double x, double y) {
			double th = NumericalAnalysis.newtonRaphsonApproximation(
					y, y/Y_SCALE, EqualEarth::poly9, EqualEarth::poly8, 1e-6);
			return new double[] { Math.asin(Math.sin(th)/B), x*B/Math.cos(th)*poly8(th) };
		}
	};
	
	
	
	private static double poly9(double x) {
		return   A[3]*Math.pow(x,9) +   A[2]*Math.pow(x,7) +   A[1]*Math.pow(x,3) + A[0]*x;
	}
	
	private static double poly8(double x) {
		return 9*A[3]*Math.pow(x,8) + 7*A[2]*Math.pow(x,6) + 3*A[1]*Math.pow(x,2) + A[0];
	}
	
}
