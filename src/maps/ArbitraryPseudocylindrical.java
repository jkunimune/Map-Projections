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

import utils.Shape;

import static java.lang.Math.PI;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;
import static utils.NumericalAnalysis.polynomialInterpolate;

/**
 * A class specifically for pseudocylindrical projections that use arbitrary tables of numbers.
 * 
 * @author jkunimune
 */
public class ArbitraryPseudocylindrical {
	
	private static final int ORDER = 3; // the order of the polynomials used
	
	
	public static final Projection ROBINSON = new ArbitraryProjection(
			"Robinson", "Arthur H. Robinson", 0.5072, new double[][] {
				{ -90,-85,-80,-75,-70,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10, -5,  0,
					 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90 }, //Latitude
				{  0.5322, 0.5722, 0.6213, 0.6732, 0.7186, 0.7597, 0.7986, 0.8350, 0.8679, 0.8962, 0.9216, 0.9427, 0.9600, 0.9730, 0.9822, 0.9900, 0.9954, 0.9986, 1.0000,
					0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322 }, //PLEN
				{ -1.0000,-0.9761,-0.9394,-0.8936,-0.8435,-0.7903,-0.7346,-0.6769,-0.6176,-0.5571,-0.4958,-0.4340,-0.3720,-0.3100,-0.2480,-0.1860,-0.1240,-0.0620, 0.0000,
					0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000 } //PDFE
			});

	
	public static final Projection NATURAL_EARTH = new ArbitraryProjection(
			"Natural Earth", "Tom Patterson", 0.520, new double[][] {
				{ -90,-85,-80,-75,-70,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10, -5,  0,
					 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90 }, //Latitude
				{  0.5630, 0.6270, 0.6754, 0.7160, 0.7525, 0.7874, 0.8196, 0.8492, 0.8763, 0.9006, 0.9222, 0.9409, 0.9570, 0.9703, 0.9811, 0.9894, 0.9953, 0.9988, 1.0000,
					0.9988, 0.9953, 0.9894, 0.9811, 0.9703, 0.9570, 0.9409, 0.9222, 0.9006, 0.8763, 0.8492, 0.8196, 0.7874, 0.7525, 0.7160, 0.6754, 0.6270, 0.5630 }, //PLEN
				{ -1.0000,-0.9761,-0.9394,-0.8936,-0.8435,-0.7903,-0.7346,-0.6769,-0.6176,-0.5571,-0.4958,-0.4340,-0.3720,-0.3100,-0.2480,-0.1860,-0.1240,-0.06200, 0.0000,
					0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000 } //PDFE
			});
	
	
	
	private static class ArbitraryProjection extends Projection {
		
		private final double[][] table;
		private final double yScale;
		
		public ArbitraryProjection(String title, String inventor, double yScale, double[][] table) {
			super(title, inventor, "A compromise pseudocylindrical projection designed by " + inventor,
			      null, true, true, true, true, Type.PSEUDOCYLINDRICAL, Property.COMPROMISE, 3);
			this.table = table;
			this.yScale = yScale;
			this.shape = Shape.meridianEnvelope(this);
		}

		public double[] project(double lat, double lon) {
			return new double[] {
					lon/PI*polynomialInterpolate(toDegrees(lat), table[0], table[1], ORDER),
					yScale*polynomialInterpolate(toDegrees(lat), table[0], table[2], ORDER) };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {
					toRadians(polynomialInterpolate(y/yScale, table[2], table[0], ORDER)),
					PI*x/polynomialInterpolate(y/yScale, table[2], table[1], ORDER) };
		}
		
	}
	
}
