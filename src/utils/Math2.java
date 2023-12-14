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
package utils;

import de.jtem.mfc.field.Complex;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Double.isFinite;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.toRadians;

/**
 * A class of some useful Math functions that seem like they could be in Math
 *
 * @author Justin Kunimune
 */
public class Math2 {

	/**
	 * Compute the arithmetic mean.
	 * @param values the numbers to be averaged
	 * @return       the average of values
	 */
	public static double mean(double[][] values) {
		double s = 0, n = 0;
		for (double[] row: values) {
			for (double x: row) {
				if (isFinite(x)) { //ignore NaN values in the average
					s += x;
					n += 1;
				}
			}
		}
		return s/n;
	}


	/**
	 * Compute the standard deviation.
	 * @param values the numbers that deviate
	 * @return       the standard deviation of values
	 */
	public static double stdDev(double[][] values) {
		double s = 0, ss = 0, n = 0;
		for (double[] row: values) {
			for (double x: row) {
				if (isFinite(x)) {
					s += x;
					ss += x*x;
					n += 1;
				}
			}
		}
		return sqrt(ss/n - s*s/(n*n));
	}


	/**
	 * Compute the root-mean-square.
	 * @param values the values to evaluate
	 * @return       the square-root of the mean of the square of values
	 */
	public static double rms(double[][] values) {
		double ss = 0, n = 0;
		for (double[] row: values) {
			for (double x: row) {
				if (isFinite(x)) {
					ss += x*x;
					n += 1;
				}
			}
		}
		return sqrt(ss/n);
	}


	/**
	 * Compute the combination
	 * @param n the number of elements from which to combine
	 * @param k the number of slots into which to combine
	 * @return  the number of ways the given elements can fill the given slots
	 */
	public static double combine(double n, int k) {
		double output = 1;
		for (int i = k; i > 0; i --) {
			output *= (n+i-k)/i;
		}
		return output;
	}


	/**
	 * Compute the modulus in a more predictable way from Java's primitive one.
	 * @param x the dividend
	 * @param y the divisor
	 * @return  the modulus x%y, but with proper handling of negative numbers
	 */
	public static double floorMod(double x, double y) {
		return x - floor(x / y) * y;
	}


	/**
	 * Like signum, but with sigone(0) = 1.
	 * @param x the value whose sign to check
	 * @return  1 if x >= 0, -1 otherwise
	 */
	public static double sigone(double x) {
		if (x >= 0)
			return 1;
		else
			return -1;
	}


	public static double coerceAngle(double ang) {
		return floorMod(ang+PI, 2*PI) - PI;
	}


	public static double determ(double a, double b, double c, double d) {
		return a*d - b*c;
	}


	public static double linInterp(double x, double a0, double a1, double b0, double b1) {
		return (x-a0)*(b1-b0)/(a1-a0) + b0;
	}
	
	
	public static double hypot(double[] a, double[] b) {
		return Math.hypot(a[0] - b[0], a[1] - b[1]);
	}


	/**
	 * like Math.max() but takes an arbitrary number of arguments
	 */
	public static double max(double... ds) {
		double m = NEGATIVE_INFINITY;
		for (double d: ds)
			if (d > m)
				m = d;
		return m;
	}


	/**
	 * like Math.min() but takes an arbitrary number of arguments
	 */
	public static double min(double... ds) {
		double m = POSITIVE_INFINITY;
		for (double d: ds)
			if (d < m)
				m = d;
		return m;
	}


	/**
	 * like Math.round(), but lets you specify a numer of decimal places to be preserved
	 */
	public static double round(double x, int numPlaces) {
		return Math.round(x*pow(10, numPlaces))/pow(10, numPlaces);
	}
	
	
	public static double sind(double angdeg) {
		return sin(toRadians(angdeg));
	}
	
	
	public static double cosd(double angdeg) {
		return cos(toRadians(angdeg));
	}
	
	
	public static double tand(double angdeg) {
		return tan(toRadians(angdeg));
	}
	
	
	public static double cotd(double angdeg) {
		return 1/tan(toRadians(angdeg));
	}
	
	
	public static double secd(double angdeg) {
		return 1/cos(toRadians(angdeg));
	}


	public static Complex atan(Complex z) {
		return new Complex(0, 1).plus(z).divide(new Complex(0, 1).minus(z)).log().times(new Complex(0, 0.5));
	}

}
