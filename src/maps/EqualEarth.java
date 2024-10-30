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
import utils.NumericalAnalysis;
import utils.Shape;

import static java.lang.Math.PI;
import static java.lang.Math.asin;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

/**
 * A fairly simple projection, but I have it in its own file, because otherwise the way the height and width are
 * defined is inconvenient.
 * 
 * @author Justin Kunimune
 */
public class EqualEarth {

	private static final double[] A = {1.340264, -.081106, .000893, .003796};
	private static final double B = sqrt(3)/2;
	private static final double Y_SCALE = poly9(PI/3)/(PI/3);
	
	
	public static final Projection EQUAL_EARTH = new Projection(
			"Equal Earth", "B. Savric, T. Patterson, and B. Jenny", "An equal-area pseudocylindrical projection specifically designed to woo Gall-Peters supporters away from that horrid thing",
			null, true, true, true, true, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3) {
		
		public double[] project(double lat, double lon) {
			double th = asin(B*sin(lat));
			return new double[] { cos(th)/B/poly8(th)*lon, poly9(th) };
		}
		
		public double[] inverse(double x, double y) {
			double th = NumericalAnalysis.newtonRaphsonApproximation(
					y, y/Y_SCALE, EqualEarth::poly9, EqualEarth::poly8, 1e-6);
			return new double[] { asin(sin(th)/B), x*B/cos(th)*poly8(th) };
		}
		
		public void initialize(double... params) throws IllegalArgumentException {
			this.shape = Shape.meridianEnvelope(this);
		}
	};
	
	
	
	private static double poly9(double x) {
		return   A[3]*pow(x,9) +   A[2]*pow(x,7) +   A[1]*pow(x,3) + A[0]*x;
	}
	
	private static double poly8(double x) {
		return 9*A[3]*pow(x,8) + 7*A[2]*pow(x,6) + 3*A[1]*pow(x,2) + A[0];
	}
	
}
