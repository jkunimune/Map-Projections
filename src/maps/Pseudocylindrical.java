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
			new Projection("Sinusoidal", "An equal-area map shaped like a sine-wave.",
					2*Math.PI, Math.PI, 0b1111, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA) {
		
		public double[] project(double lat, double lon) {
			return new double[] { Math.cos(lat)*lon, lat };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { y, x/Math.cos(y) };
		}
	};
	
	
	public static final Projection MOLLWEIDE =
			new Projection("Mollweide", "An equal-area projection shaped like an ellipse.",
					2*Math.PI, Math.PI, 0b1101, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA) {
		
		public double[] project(double lat, double lon) {
			double tht = lat;
			for (int i = 0; i < 10; i ++)
				tht -= (2*tht+Math.sin(2*tht)-Math.PI*Math.sin(lat))/
						(2+2*Math.cos(2*tht));
			return new double[] { lon*Math.cos(tht), Math.sin(tht)*Math.PI/2 };
		}
		
		public double[] inverse(double x, double y) {
			double tht = Math.asin(2*y/Math.PI);
			return new double[] {
					Math.asin((2*tht + Math.sin(2*tht)) / Math.PI),
					x / Math.cos(tht)};
		}
	};
	
	
	public static final Projection LEMONS =
			new Projection(
					"Lemons", "BURN LIFE'S HOUSE DOWN!!!", 2*Math.PI, Math.PI, 0b1110,
					Type.OTHER, Property.COMPROMISE) {
		
		private static final int NUM_LEMONS = 12; //number of lemons
		private static final double LEM_WIDTH = 2*Math.PI/NUM_LEMONS; //longitude span of 1 lemon
		
		public double[] project(double lat, double lon) {
			final int lemNum = (int)Math.floor(lon/LEM_WIDTH);
			final double dl = (lon+2*Math.PI) % LEM_WIDTH - LEM_WIDTH/2;
			return new double[] {
					Math.asin(Math.cos(lat)*Math.sin(dl)) + (lemNum+.5)*LEM_WIDTH,
					Math.asin(Math.sin(lat)/Math.sqrt(1-Math.pow(Math.cos(lat)*Math.sin(dl), 2)))};
		}
		
		public double[] inverse(double x, double y) {
			final int lemNum = (int)Math.floor(x/LEM_WIDTH);
			final double dx = (x+2*Math.PI) % LEM_WIDTH - LEM_WIDTH/2;
			final double dl = Math.asin(
					Math.sin(dx)/Math.sqrt(1-Math.pow(Math.cos(dx)*Math.sin(y), 2)));
			if (Math.abs(dl) > LEM_WIDTH/2)
				return null;
			else
				return new double[] {
						Math.asin(Math.cos(dx)*Math.sin(y)), dl + (lemNum+.5)*LEM_WIDTH };
		}
	};
}
