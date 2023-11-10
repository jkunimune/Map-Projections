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
import utils.BoundingBox;
import utils.NumericalAnalysis;

import java.util.Arrays;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.asin;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static utils.Math2.linInterp;

/**
 * A class of values and functions used to approximate the Tobler projection
 * 
 * @author jkunimune
 */
public class Tobler {
	
	public static final Projection TOBLER =
			new Projection(
					"Tobler hyperelliptical", "An equal-area projection shaped like a hyperellipse.",
					new BoundingBox(2*PI, PI), 0b1001, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 4,
					new String[]{"alpha","K"},
					new double[][] {{0,1,0.0}, {1,5,2.5}}) {
		
		private static final int N = 20000;
		private double alpha, kappa, epsilon; //epsilon is related to gamma, but defined somewhat differently
		private double[] Z; //Z[i] = sin(phi) when y = i/(Z.length-1)
		
		public void initialize(double... params) {
			this.alpha = params[0];
			this.kappa = params[1];
			this.epsilon = NumericalAnalysis.simpsonIntegrate(
					0, 1, this::hyperEllipse, 1./N);
			this.Z = NumericalAnalysis.simpsonODESolve(
					1, N, this::dZdY, 1./N);
		}
		
		public double[] project(double lat, double lon) {
			final double z0 = abs(sin(lat));
			final int i = Arrays.binarySearch(Z, z0);
			final double y;
			if (i >= 0)
				y = i/(Z.length-1.);
			else if (-i-1 >= Z.length)
				y = Z[Z.length-1];
			else
				y = linInterp(z0, Z[-i-2], Z[-i-1], -i-2, -i-1)/
						(Z.length-1.);
			return new double[] {
					lon * abs(alpha + (1-alpha)*hyperEllipse(y)),
					y * signum(lat)*PI/2 };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {
					asin(Z[(int)round(abs(2*y/PI)*(Z.length-1))])*signum(y),
					x / abs(alpha + (1-alpha)*hyperEllipse(2*y/PI)) };
		}
		
		public double dZdY(double y) {
			return abs((alpha + (1-alpha)*hyperEllipse(y))/
					(alpha + (1-alpha)*epsilon));
		}
		
		public double hyperEllipse(double y) {
			return pow(1 - pow(abs(y),kappa), 1/kappa);
		}
	};
}
