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
package util;

import java.util.function.BinaryOperator;

/**
 * A class of some useful Math functions that seem like they could be in Math
 * 
 * @author Justin Kunimune
 */
public class Math2 {

	public static final double mean(double[][] values) {
		double s = 0, n = 0;
		for (double[] row: values) {
			for (double x: row) {
				if (Double.isFinite(x)) { //ignore NaN values in the average
					s += x;
					n += 1;
				}
			}
		}
		return s/n;
	}
	
	
	public static final double stdDev(double[][] values) {
		double s = 0, ss = 0, n = 0;
		for (double[] row: values) {
			for (double x: row) {
				if (Double.isFinite(x)) {
					s += x;
					ss += x*x;
					n += 1;
				}
			}
		}
		return Math.sqrt(ss/n - s*s/n*n);
	}
	
	
	public static final double combine(double n, int k) {
		double output = 1;
		for (int i = k; i > 0; i --) {
			output *= (n+i-k)/i;
		}
		return output;
	}
	
	
	public static final double mod(double x, double y) {
		return x - Math.floor(x/y)*y;
	}
	
	
	public static final double[] newtonRaphsonApproximation(double x0, double y0, BinaryOperator<Double> f1pX, BinaryOperator<Double> f2pY, BinaryOperator<Double> df1dp,
			BinaryOperator<Double> df1dl, BinaryOperator<Double> df2dp, BinaryOperator<Double> df2dl, double tolerance) {
		double x = x0;
		double y = y0;
		double phi = y;
		double lam = x; // I used equirectangular for my initial guess
		double error = Math.PI;
		
		for (int i = 0; i < 6 && error > tolerance; i++) {
			final double f1 = f1pX.apply(phi, lam) - x;
			final double f2 = f2pY.apply(phi, lam) - y;
			final double df1dP = df2dp.apply(phi, lam);
			final double df1dL = df1dl.apply(phi, lam);
			final double df2dP = df2dp.apply(phi, lam);
			final double df2dL = df2dl.apply(phi, lam);
			
			phi -= (f1*df2dL - f2*df1dL) / (df1dP*df2dL - df2dP*df1dL);
			lam -= (f2*df1dP - f1*df2dP) / (df1dP*df2dL - df2dP*df1dL);
			
			error = Math.hypot(f1, f2);
		}
		
		if (error > tolerance) // if it aborted due to timeout
			return null;
		else // if it aborted due to convergence
			return new double[] {phi, lam};
	}

}
