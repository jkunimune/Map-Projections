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

import org.apache.commons.math3.complex.Complex;

import static java.lang.Math.abs;
import static utils.Math2.combine;

/**
 * Actually just the one incomplete integral. I honestly don't remember where I got this sequence from
 * 
 * @author jkunimune
 */
public class Elliptic {

	public static Complex F(Complex phi, final Complex k) { //series solution to incomplete elliptic integral of the first kind
		final double TOLERANCE = 1e-3;
		
		Complex sum = Complex.ZERO;
		Complex i_n = phi;
		Complex delt;
		
		int n = 0;
		do {
			if (n > 0)
				i_n = i_n.multiply((2.0 * n - 1) / (2.0 * n))
						.subtract(phi.cos().multiply(phi.sin().pow(2.0 * n - 1)).divide(2.0 * n));
			delt = i_n.multiply(abs(combine(-.5, n))).multiply(k.pow(2.0 * n));
			sum = sum.add(delt);
			n ++;
		} while (delt.abs() > TOLERANCE);
		
		return sum;
	}

}
