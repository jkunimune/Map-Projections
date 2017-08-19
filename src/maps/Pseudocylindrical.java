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

import maps.Projection.Property;
import maps.Projection.Type;

/**
 * Projections where y is a function of latitude
 * 
 * @author jkunimune
 */
public class Pseudocylindrical {
	
	public static final Projection SINUSOIDAL =
			new Projection("Sinusoidal", "An equal-area map shaped like a sine-wave",
					2., 0b1111, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA) {
		
		public double[] project(double lat, double lon) {
			return new double[] { Math.cos(lat)*lon, lat };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { y*Math.PI/2, x*Math.PI/Math.cos(y*Math.PI/2) };
		}
	};
	
	
	public static final Projection MOLLWEIDE =
			new Projection("Mollweide", "An equal-area projection shaped like an ellipse",
					2., 0b1101, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA) {
		
		public double[] project(double lat, double lon) {
			double tht = lat;
			for (int i = 0; i < 10; i ++)
				tht -= (2*tht+Math.sin(2*tht)-Math.PI*Math.sin(lat))/
						(2+2*Math.cos(2*tht));
			return new double[] { lon*Math.cos(tht), Math.PI/2*Math.sin(tht) };
		}
		
		public double[] inverse(double x, double y) {
			double tht = Math.asin(y);
			return new double[] {
					Math.asin((2*tht + Math.sin(2*tht)) / Math.PI),
					Math.PI * x / Math.cos(tht)};
		}
	};
	
	
	public static final Projection LEMONS =
			new Projection(
					"Lemons", "BURN LIFE'S HOUSE DOWN!!!", 2., 0b1110,
					Type.PSEUDOCYLINDRICAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			return null; //TODO: projection wishlist
		}
		
		public double[] inverse(double x, double y) {
			x = x+2;
			final double lemWdt = 1/6.0;
			if (Math.abs(x % lemWdt - lemWdt / 2.0) <= Math.cos(y*Math.PI/2) * lemWdt/2.0) // if it is in
				return new double[] { y*Math.PI/2,	// a sine curve
						Math.PI * (x%lemWdt - lemWdt/2.0) / (Math.cos(y*Math.PI/2))
								+ (int)(x/lemWdt) * Math.PI/6 };
			else
				return null;
		}
	};
}
