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

/**
 * A class of values and functions used to approximate the Tobler projection
 * 
 * @author jkunimune
 */
public class Tobler {

	public static final double lam(double x, double y, double[] params) {
		final double a = params[0];
		return Math.PI * x / Math.abs(a + (1-a)*hyperEllipse(y, params));
	}
	
	
	public static final double X(double y, double lam, double[] params) {
		final double a = params[0];
		return lam * Math.abs(a + (1-a)*hyperEllipse(y, params));
	}
	
	
	public static final double dZdY(double y, double[] params) {
		final double a = params[0], e = params[2];
		return Math.abs((a + (1-a)*hyperEllipse(y, params))/
				(a + (1-a)*e));
	}
	
	
	public static final double hyperEllipse(double y, double[] params) {
		final double k = params[1];
		return Math.pow(1 - Math.pow(Math.abs(y),k), 1/k);
	}

}
