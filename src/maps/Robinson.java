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

import java.util.Arrays;

import maps.Projection.Property;
import maps.Projection.Type;
import utils.NumericalAnalysis;

/**
 * A class with static methods and a matrix constant that pertain to the
 * Robinson projection.
 * 
 * @author jkunimune
 */
public class Robinson {
	
	private static final double Y_MAX = 0.5072;
	
	private static final double[][] TABLE = {
			{ -90,-85,-80,-75,-70,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-05, 00,
				05, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90 }, //Latitude
			{  0.5322, 0.5722, 0.6213, 0.6732, 0.7186, 0.7597, 0.7986, 0.8350, 0.8679, 0.8962, 0.9216, 0.9427, 0.9600, 0.9730, 0.9822, 0.9900, 0.9954, 0.9986, 1.0000,
				0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322 }, //PLEN
			{ -1.0000,-0.9761,-0.9394,-0.8936,-0.8435,-0.7903,-0.7346,-0.6769,-0.6176,-0.5571,-0.4958,-0.4340,-0.3720,-0.3100,-0.2480,-0.1860,-0.1240,-0.0620, 0.0000,
				0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000 } //PDFE
	};
	
	private static final int ORDER = 2; //half the order of the polynomials used
	
	
	public static final Projection ROBINSON =
			new Projection("Robinson", "A visually pleasing piecewise compromise map",
					2, 2*Y_MAX, 0b1111, Type.PSEUDOCYLINDRICAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			return new double[] {
					lon/Math.PI*smartInterpolate(Math.toDegrees(lat), TABLE[0], TABLE[1], ORDER),
					Y_MAX*smartInterpolate(Math.toDegrees(lat), TABLE[0], TABLE[2], ORDER) };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {
					Math.toRadians(smartInterpolate(y/Y_MAX, TABLE[2], TABLE[0], ORDER)),
					Math.PI*x/smartInterpolate(y/Y_MAX, TABLE[2], TABLE[1], ORDER) };
		}
	};
	
	
	public static final double smartInterpolate(double x, double[] X, double[] f, int k) {
		int i = Arrays.binarySearch(X, x);
		if (i < 0)	i = -i - 1; //if you couldn't find it, don't worry about it
		return NumericalAnalysis.aitkenInterpolate(x, X, f,
				Math.max(i-k,0), Math.min(i+k,X.length)); //just call aitken with the correct bounds
	}

}
