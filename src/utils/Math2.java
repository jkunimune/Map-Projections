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

/**
 * A class of some useful Math functions that seem like they could be in Math
 * 
 * @author Justin Kunimune
 */
public class Math2 {

	public static final double PHI = (1+Math.sqrt(5))/2;
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
		return Math.sqrt(ss/n - s*s/(n*n));
	}
	
	
	public static final double rms(double[][] values) {
		double ss = 0, n = 0;
		for (double[] row: values) {
			for (double x: row) {
				if (Double.isFinite(x)) {
					ss += x*x;
					n += 1;
				}
			}
		}
		return Math.sqrt(ss/n);
	}
	
	
	public static final double combine(double n, int k) {
		double output = 1;
		for (int i = k; i > 0; i --) {
			output *= (n+i-k)/i;
		}
		return output;
	}
	
	
	public static final double floorMod(double x, double y) {
		return x - Math.floor(x / y) * y;
	}
	
	
	public static final double determ(double a, double b, double c, double d) {
		return a*d - b*c;
	}
	
	
	public static final double linInterp(double x, double a0, double a1, double b0, double b1) {
		return (x-a0)*(b1-b0)/(a1-a0) + b0;
	}
	
	public static double[] linInterp(double[] xs, double[][] A, double[][] B) {
		double[] out = new double[xs.length];
		for (int i = 0; i < xs.length; i ++)
			out[i] = linInterp(xs[i], A[0][i], A[1][i], B[0][i], B[1][i]);
		return out;
	}
	
	
	public static double
			lineSegmentDistance(double X, double Y, double x0, double y0, double x1, double y1) {
		double A = X - x0;
		double B = Y - y0;
		double C = x1 - x0;
		double D = y1 - y0;
		
		double dot = A * C + B * D;
		double len_sq = C * C + D * D;
		double param = -1;
		if (len_sq != 0) // in case of 0 length line
			param = dot / len_sq;
		
		double xc, yc;
		if (param < 0) {
			xc = x0;
			yc = y0;
		}
		else if (param > 1) {
			xc = x1;
			yc = y1;
		}
		else {
			xc = x0 + param * C;
			yc = y0 + param * D;
		}
		return Math.hypot(X - xc, Y - yc);
	}


	public static final double[] invArr(double[] ds) {
		if (ds == null) 	return null; //this is kind of a weird place to put this catch, but whatever
		for (int i = 0; i < ds.length; i ++)
			ds[i] *= -1;
		return ds;
	}
	
	
	public static double max(double... ds) {
		double m = Double.NEGATIVE_INFINITY;
		for (double d: ds)
			if (d > m)
				m = d;
		return m;
	}
	
	
	public static double min(double... ds) {
		double m = Double.POSITIVE_INFINITY;
		for (double d: ds)
			if (d < m)
				m = d;
		return m;
	}
	
	
	public static double round(double x, int numPlaces) {
		return Math.round(x*Math.pow(10, numPlaces))/Math.pow(10, numPlaces);
	}
	
	
	public static double sind(double angdeg) {
		return Math.sin(Math.toRadians(angdeg));
	}
	
	
	public static double cosd(double angdeg) {
		return Math.cos(Math.toRadians(angdeg));
	}


	public static double tand(double angdeg) {
		return Math.tan(Math.toRadians(angdeg));
	}


	public static double cotd(double angdeg) {
		return 1/Math.tan(Math.toRadians(angdeg));
	}


	public static double removeRoundErrorZeros(double d) { //this is a really specific function to deal with round-off error for doubles
//		if (Math.abs(d) <= Math.pow(2, -53))
//			return 0;
//		else
			return d;
	}

}
