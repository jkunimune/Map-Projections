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

/**
 * A class with a few handy Dixon elliptic functions as they pertain to the Lee conformal projection.
 * All the algorithms here came directly from L.P. Lee's paper,
 * "Some Conformal Projections Based on Elliptic Functions"
 * 
 * Lee, L. P. "Some Conformal Projections Based on Elliptic Functions."
 *  	Geographical Review, vol. 55, no. 4, 1965, pp. 563-580. JSTOR, 
 *  	www.jstor.org/stable/212415.
 * 
 * @author jkunimune
 */
public class Dixon {

	private static final double TOLERANCE = 1e-3;
	
	
	public static Complex leeFunc(Complex w) { //the 28th order McLaurin polynomial for 2sm(w/2)cm(w/2)
		return w.plus(w.pow( 4).times(.625000e-1))
				.plus( w.pow( 7).times(.223214e-2))
				.plus(w.pow(10).times(.069754e-3))
				.plus( w.pow(13).times(.020121e-4))
				.plus(w.pow(16).times(.005539e-5));
//				.plus( w.pow(19).times(.001477e-6))
//				.plus(w.pow(22).times(.000385e-7))
//				.plus( w.pow(25).times(.000099e-8))
//				.plus(w.pow(28).times(.000025e-9));
	}
	
	
	public static Complex invFunc(Complex z) { //the iterative algorithm specifically suggested by Lee
		Complex wi;
		Complex wf = z;
		
		do {
			wi = wf;
			wf = z.plus(wi.minus(leeFunc(wi)));
		} while (wf.minus(wi).abs() > TOLERANCE);
		
		return wf;
	}

}
