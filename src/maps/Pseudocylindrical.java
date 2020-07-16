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
import utils.NumericalAnalysis;

/**
 * Projections where y is a function of latitude
 * 
 * @author jkunimune
 */
public class Pseudocylindrical {
	
	public static final Projection SINUSOIDAL = new Projection(
			"Sinusoidal", "An equal-area map shaped like a sine-wave.",
			2*Math.PI, Math.PI, 0b1111, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 1) {
		
		public double[] project(double lat, double lon) {
			return new double[] { Math.cos(lat)*lon, lat };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { y, x/Math.cos(y) };
		}
	};
	
	
	public static final Projection MOLLWEIDE = new Projection(
			"Mollweide", "An equal-area projection shaped like an ellipse.",
			4, 2, 0b1101, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3) {
		
		public double[] project(double lat, double lon) {
			double tht = NumericalAnalysis.newtonRaphsonApproximation(
					Math.PI*Math.sin(lat), lat,
					(t) -> (2*t + Math.sin(2*t)),
					(t) -> (2 + 2*Math.cos(2*t)), 1e-6);
			if (Double.isNaN(tht))
				tht = Math.PI/2*Math.signum(lat);
			return new double[] { lon/Math.PI*2*Math.cos(tht), Math.sin(tht) };
		}
		
		public double[] inverse(double x, double y) {
			double tht = Math.asin(y);
			return new double[] {
					Math.asin((2*tht + Math.sin(2*tht))/Math.PI),
					x/Math.cos(tht)*Math.PI/2 };
		}
	};
	
	
	public static final Projection HOMOLOSINE = new Projection(
			"Homolosine (uninterrupted)", "A combination of the sinusoidal and Mollweide projections.",
			2*Math.PI, 2.72282, 0b1101, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3) {
		
		private final double phiH = 0.71098;
		private final double scale = Math.sqrt(2);
		private final double yH = MOLLWEIDE.project(phiH, 0)[1]*scale;
		
		public double[] project(double lat, double lon) {
			if (Math.abs(lat) <= phiH) {
				return SINUSOIDAL.project(lat, lon);
			}
			else {
				double[] xy = MOLLWEIDE.project(lat, lon);
				if (lat > 0)
					return new double[] {xy[0]*scale, xy[1]*scale + phiH - yH};
				else
					return new double[] {xy[0]*scale, xy[1]*scale - phiH + yH};
			}
		}
		
		public double[] inverse(double x, double y) {
			if (Math.abs(y) <= phiH)
				return SINUSOIDAL.inverse(x, y);
			else if (y > 0)
				return MOLLWEIDE.inverse(x/scale, (y - phiH + yH)/scale);
			else
				return MOLLWEIDE.inverse(x/scale, (y + phiH - yH)/scale);
		}
	};
	
	
	public static final Projection HOMOLOSINE_INTERRUPTED = new Projection(
			"Good Homolosine", "An interrupted combination of the sinusoidal and Mollweide projections.",
			2*Math.PI, 2.72282, 0b1100, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3) {
		
		private final double[][] edges = {
				{Math.toRadians(-40), Math.toRadians(180)},
				{Math.toRadians(-100), Math.toRadians(-20), Math.toRadians(80), Math.toRadians(180)}};
		private final double[][] centers = {
				{Math.toRadians(-100), Math.toRadians(30)},
				{Math.toRadians(-160), Math.toRadians(-60), Math.toRadians(20), Math.toRadians(140)}};
		
		public double[] project(double lat, double lon) {
			int i = (lat > 0) ? 0 : 1;
			for (int j = 0; j < edges[i].length; j ++) {
				if (lon <= edges[i][j]) {
					double[] xy = HOMOLOSINE.project(lat, lon - centers[i][j]);
					return new double[] {xy[0] + centers[i][j], xy[1]};
				}
			}
			return null;
		}
		
		public double[] inverse(double x, double y) {
			int i = (y > 0) ? 0 : 1;
			for (int j = 0; j < edges[i].length; j ++) {
				if (x <= edges[i][j]) {
					double[] phiLam = HOMOLOSINE.inverse(x - centers[i][j], y);
					if ((j < edges[i].length-1 && phiLam[1] + centers[i][j] > edges[i][j]) ||
							(j > 0 && phiLam[1] + centers[i][j] < edges[i][j-1]))
						return null;
					else
						return new double[] {phiLam[0], phiLam[1] + centers[i][j]};
				}
			}
			return null;
		}
	};
	
	
	public static final Projection ECKERT_IV = new Projection(
			"Eckert IV", "An equal-area projection released in a set of six (I'm only giving you the one because the others are not good).",
			4, 2, 0b1101, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3) {
		
		public double[] project(double lat, double lon) {
			double tht = NumericalAnalysis.newtonRaphsonApproximation(
					(2+Math.PI/2)*Math.sin(lat), lat,
					(t) -> (t + Math.sin(2*t)/2 + 2*Math.sin(t)),
					(t) -> (1 + Math.cos(2*t) + 2*Math.cos(t)), 1e-4);
			return new double[] { lon/Math.PI*(1+Math.cos(tht)), Math.sin(tht)};
		}
		
		public double[] inverse(double x, double y) {
			double tht = Math.asin(y);
			return new double[] {
					Math.asin((tht + Math.sin(2*tht)/2 + 2*Math.sin(tht))/(2+Math.PI/2)),
					x/(1 + Math.cos(tht))*Math.PI };
		}
		
	};
	
	
	public static final Projection WAGNER_II = new Projection(
			"Wagner II", "A compromise projection with sinusoidal meridians.",
			2*2.9054, 2*1.4527, 0b1111, Type.OTHER, Property.COMPROMISE, 2) {
		
		private final double c0 = 0.92483, c1 = 1.38725,
				c2 = 0.88022, c3 = 0.8855;
		
		public double[] project(double lat, double lon) {
			double psi = Math.asin(c2*Math.sin(c3*lat));
			return new double[] {c0*lon*Math.cos(psi), c1*psi};
		}
		
		public double[] inverse(double x, double y) {
			double psi = y/c1;
			return new double[] {Math.asin(Math.sin(psi)/c2)/c3, x/c0/Math.cos(psi)};
		}
	};
	
	
	public static final Projection WAGNER_V = new Projection(
			"Wagner V", "A compromise projection with elliptical meridians.",
			2*2.8581, 2*1.4291, 0b1111, Type.OTHER, Property.COMPROMISE, 3) {
		
		private final double c0 = 0.909771, c1 = 1.650142,
				c2 = 3.008957, c3 = 0.8855;
		
		public double[] project(double lat, double lon) {
			double psi = NumericalAnalysis.newtonRaphsonApproximation(
					c2*Math.sin(c3*lat), Math.sin(lat)*Math.PI/3, (ps)->(2*ps + Math.sin(2*ps)),
					(ps)->(2 + 2*Math.cos(2*ps)), 1e-5);
			return new double[] {c0*lon*Math.cos(psi), c1*Math.sin(psi)};
		}
		
		public double[] inverse(double x, double y) {
			double psi = Math.asin(y/c1);
			return new double[] {Math.asin((2*psi + Math.sin(2*psi))/c2)/c3, x/c0/Math.cos(psi)};
		}
	};
	
	
	public static final Projection KAVRAYSKIY_VII = new Projection(
			"Kavrayskiy VII", Math.PI*Math.sqrt(3), Math.PI, 0b1111, Type.PSEUDOCYLINDRICAL,
			Property.COMPROMISE, 2, null, "mostly popular in the former Soviet Union") {
		
		public double[] project(double lat, double lon) {
			return new double[] { 1.5*lon*Math.sqrt(1/3.-Math.pow(lat/Math.PI, 2)), lat };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { y, x/1.5/Math.sqrt(1/3.-Math.pow(y/Math.PI, 2)) };
		}
		
	};
	
	
	public static final Projection LEMONS = new Projection(
			"Lemons", "BURN LIFE'S HOUSE DOWN!!!", 2*Math.PI, Math.PI, 0b1110,
			Type.CYLINDRICAL, Property.COMPROMISE, 2) {
		
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
