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
 * 
 * 
 * @author jkunimune
 */
public class Aitoff {
	
	public static final double f1pY(double phi, double lam) {
		return Math.sin(phi)/Math.sqrt(P(phi,lam));
	}
	
	public static final double f2pX(double phi, double lam) {
		return 2*Math.cos(phi)*Math.sin(lam/2)/Math.sqrt(P(phi,lam));
	}
	
	public static final double df1dphi(double phi, double lam) {
		final double p = P(phi,lam);
		return 1/2.*((1+p)*Math.cos(phi) + Math.cos(lam/2))/Math.pow(p,1.5);
	}
	
	public static final double df1dlam(double phi, double lam) {
		final double p = P(phi,lam);
		return 1/4.*Math.sin(phi)*Math.cos(phi)*Math.sin(lam/2)/Math.pow(p,1.5);
	}
	
	public static final double df2dphi(double phi, double lam) {
		final double p = P(phi,lam);
		return -(1+p)*Math.sin(phi)*Math.sin(lam/2)/Math.pow(p,1.5);
	}
	
	public static final double df2dlam(double phi, double lam) {
		final double p = P(phi,lam);
		return 1/2.*((1+p)*Math.cos(phi)*Math.cos(lam/2) + Math.pow(Math.cos(phi),2))/Math.pow(p,1.5);
	}
	
	private static final double P(double phi, double lam) {
		return 1 + Math.cos(phi)*Math.cos(lam/2);
	}

}
