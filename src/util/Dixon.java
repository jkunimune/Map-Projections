/**
 * 
 */
package util;

import mfc.field.Complex;

/**
 * A class with a few handy Dixon elliptic functions as they pertain to the Lee conformal projection.
 * All the algorithms here came directly from L.P. Lee's paper,
 * "Some Conformal Projections Based on Elliptic Functions"
 * 
 * Lee, L. P. “Some Conformal Projections Based on Elliptic Functions.”
 *  	Geographical Review, vol. 55, no. 4, 1965, pp. 563–580. JSTOR, 
 *  	www.jstor.org/stable/212415.
 * 
 * @author jkunimune
 */
public class Dixon {

	private static final double TOLERANCE = 1e-2;//TODO
	
	
	public static Complex leeFunc(Complex w) { //the 28th order McLaurin polynomial for 2sm(w/2)cm(w/2)
		return w.minus(w.pow( 4).times(.625000e-1))
				.plus( w.pow( 7).times(.223214e-2))
				.minus(w.pow(10).times(.069754e-3))
				.plus( w.pow(13).times(.020121e-4))
				.minus(w.pow(16).times(.005539e-5))
				.plus( w.pow(19).times(.001477e-6))
				.minus(w.pow(22).times(.000385e-7))
				.plus( w.pow(25).times(.000099e-8))
				.minus(w.pow(28).times(.000025e-9));
	}
	
	
	public static Complex invFunc(Complex z) { //the iterative algorithm specifically suggested by Lee
		Complex wi;
		Complex wf = z;
		
		//do {
		for (int i = 0; i < 5; i ++) {
			wi = wf;
			wf = z.plus(wi.minus(leeFunc(wi)));
		}
		//} while (wf.minus(wi).abs() > TOLERANCE);
		
		return wf;
	}

}
