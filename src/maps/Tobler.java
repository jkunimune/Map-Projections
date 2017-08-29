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
import utils.Math2;
import utils.NumericalAnalysis;

/**
 * A class of values and functions used to approximate the Tobler projection
 * 
 * @author jkunimune
 */
public class Tobler {
	
	public static final Projection TOBLER =
			new Projection(
					"Tobler", "An equal-area projection shaped like a hyperellipse (in case you're wondering about gamma, it's calculated automatically)",
					2., 0b1001, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA,
					new String[]{"Std. Parallel","alpha","K"},
					new double[][] {{0,89,37}, {0,1,0}, {1,8,3.6}}) { //optimal parameters are 30.6,.50,3.63, but these defaults are more recognizably Tobler
		
		private static final int N = 10000;
		private double alpha, kappa, epsilon; //epsilon is related to gamma, but defined somewhat differently
		private double[] Z; //Z[i] = sin(phi) when y = i/(Z.length-1)
		
		public void setParameters(double... params) {
			this.aspectRatio = Math.pow(Math.cos(Math.toRadians(params[0])),2)*Math.PI;
			this.alpha = params[1];
			this.kappa = params[2];
			this.epsilon = NumericalAnalysis.simpsonIntegrate(
					0, 1, this::hyperEllipse, 1./N);
			this.Z = NumericalAnalysis.simpsonODESolve(
					1, N, this::dZdY, 1./N);
		}
		
		public double[] project(double lat, double lon) {
			final double z0 = Math.abs(Math.sin(lat));
			final int i = Arrays.binarySearch(Z, z0);
			final double y;
			if (i >= 0)
				y = i/(Z.length-1.);
			else if (-i-1 >= Z.length)
				y = Z[Z.length-1];
			else
				y = Math2.linInterp(z0, Z[-i-2], Z[-i-1], -i-2, -i-1)/
						(Z.length-1.);
			return new double[] {
					lon * Math.abs(alpha + (1-alpha)*hyperEllipse(y))/Math.PI,
					y * Math.signum(lat)/aspectRatio };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] {
					Math.asin(Z[(int)Math.round(Math.abs(y)*(Z.length-1))])*Math.signum(y),
					Math.PI * x / Math.abs(alpha + (1-alpha)*hyperEllipse(y)) };
		}
		
		public double dZdY(double y) {
			return Math.abs((alpha + (1-alpha)*hyperEllipse(y))/
					(alpha + (1-alpha)*epsilon));
		}
		
		public double hyperEllipse(double y) {
			return Math.pow(1 - Math.pow(Math.abs(y),kappa), 1/kappa);
		}
	};
}
