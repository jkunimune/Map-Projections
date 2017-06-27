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

import java.util.Arrays;
import java.util.function.BinaryOperator;

/**
 * A whole class just for numeric approximation methods
 * 
 * @author jkunimune
 */
public class NumericalAnalysis {

	/**
	 * Applies Newton's method in two dimensions to solve for phi and lam given
	 * desired x and y values, x(phi,lam), y(phi,lam), and some derivatives
	 * @param x0 The x value that the functions need to match
	 * @param y0 The y value that the functions need to match
	 * @param f1pX x in terms of phi and lam
	 * @param f2pY y in terms of phi and lam
	 * @param df1dp The partial derivative of x with respect to phi
	 * @param df1dl The partial derivative of x with respect to lam
	 * @param df2dp The partial derivative of y with respect to phi
	 * @param df2dl The partial derivative of y with respect to lam
	 * @param tolerance The maximum error that this can return
	 * @return the values of phi and lam that put f1pX and f2pY near x0 and y0
	 */
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
	
	
	/**
	 * Applies aitken interpolation to an array of tabulated values
	 * @param x The input value
	 * @param X The sorted array of inputs on which to interpolate
	 * @param f The sorted array of outputs on which to interpolate
	 * @return f(x), approximately
	 */
	public static final double aitkenInterpolate(double x, double[] X, double[] f) {
		return aitkenInterpolate(x, X, f, 0, X.length);
	}
	
	/**
	 * Applies aitken interpolation to a subset of an array of tabulated values
	 * @param x The input value
	 * @param X The sorted array of inputs on which to interpolate
	 * @param f The sorted array of outputs on which to interpolate
	 * @param from The index of the arrays at which to start
	 * @param to The index of the arrays at which to stop
	 * @return f(x), approximately
	 */
	public static final double aitkenInterpolate(double x,
			double[] X, double[] f, int from, int to) { //map from X to f using elements start to end
		final int N = to - from;
		final double[][] fx = new double[N][]; // the table of successive approximations
		
		fx[0] = Arrays.copyOfRange(f, from, to); //fill in the zeroth row
		
		for (int i = 1; i < N; i ++) { //i+1 is the number of points interpolated on
			fx[i] = new double[N];
			for (int j = i; j < N; j ++) { //the points will be 0, ..., i-1, j
				fx[i][j] = 1/(X[from+j] - X[from+i-1])*Math2.determ(
						fx[i-1][i-1], fx[i-1][j],
						X[from+i-1] - x, X[from+j] - x); //I won't bother to explain this; go look up Aitken interpolation
			}
		}
		
		return fx[N-1][N-1];
	}
	
	
	public static final void main(String[] args) {
		System.out.println("Testing Aitken Interpolation:");
		System.out.println("Input points are (-1,1), (-.5,-1), (0,0), (.5,1), and (1,-1)");
		double[] X = {-1, -.5, 0, .5, 1};
		double[] Y = {1, -1, 0, 1, -1};
		System.out.println("Interpolating points from -1 to 1:");
		for (double x = -1; x <= 1; x += 1/1024.)
			System.out.println(x+", "+aitkenInterpolate(x, X, Y)+";");
	}

}
