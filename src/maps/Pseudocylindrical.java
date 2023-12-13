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

import image.Path;
import maps.Projection.Property;
import maps.Projection.Type;
import utils.NumericalAnalysis;
import utils.Shape;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Double.isNaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.asin;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.toRadians;

/**
 * Projections where y is a function of latitude
 * 
 * @author jkunimune
 */
public class Pseudocylindrical {
	
	public static final Projection SINUSOIDAL = new Projection(
			"Sinusoidal", "An equal-area map shaped like a sine-wave",
			null, true, true, true, true, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 1) {
		public void initialize(double... params) {
			this.shape = Shape.meridianEnvelope(this);
		}

		public double[] project(double lat, double lon) {
			return new double[] { cos(lat)*lon, lat };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { y, x/cos(y) };
		}
	};


	public static final Projection MOLLWEIDE = new Projection(
			"Mollweide", "An equal-area projection shaped like an ellipse",
			Shape.ellipse(2, 1), true, true, false, true,
			Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3) {

		public double[] project(double lat, double lon) {
			double tht = NumericalAnalysis.newtonRaphsonApproximation(
					PI*sin(lat), lat,
					(t) -> (2*t + sin(2*t)),
					(t) -> (2 + 2*cos(2*t)), 1e-6);
			if (isNaN(tht))
				tht = PI/2*signum(lat);
			return new double[] { lon/PI*2*cos(tht), sin(tht) };
		}
		
		public double[] inverse(double x, double y) {
			double tht = asin(y);
			return new double[] {
					asin((2*tht + sin(2*tht))/PI),
					x/cos(tht)*PI/2 };
		}
	};
	
	
	public static final Projection HOMOLOSINE = new Projection(
			"Homolosine (uninterrupted)", "A combination of the sinusoidal and Mollweide projections",
			null, true, true, false, true, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3) {
		
		private final double phiH = 0.71098;
		private final double scale = sqrt(2);
		private final double yH = MOLLWEIDE.project(phiH, 0)[1]*scale;

		public void initialize(double... params) {
			this.shape = Shape.meridianEnvelope(this);
		}

		public double[] project(double lat, double lon) {
			if (abs(lat) <= phiH) {
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
			if (abs(y) <= phiH)
				return SINUSOIDAL.inverse(x, y);
			else if (y > 0)
				return MOLLWEIDE.inverse(x/scale, (y - phiH + yH)/scale);
			else
				return MOLLWEIDE.inverse(x/scale, (y + phiH - yH)/scale);
		}
	};
	
	
	public static final Projection HOMOLOSINE_INTERRUPTED = new Projection(
			"Goode Homolosine", "An interrupted combination of the sinusoidal and Mollweide projections",
			null, false, true, false, true, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3) {
		
		private final double[][] edges = {
				{toRadians(-40), toRadians(180)},
				{toRadians(-100), toRadians(-20), toRadians(80), toRadians(180)}};
		private final double[][] centers = {
				{toRadians(-100), toRadians(30)},
				{toRadians(-160), toRadians(-60), toRadians(20), toRadians(140)}};

		public void initialize(double... params) {
			// to calculate the shape, start by getting the shape of a generic meridian
			List<Path.Command> polewardSegment = HOMOLOSINE.drawLoxodrome(0, 1, PI/2, 1, .1);
			// make an equator-to-pole version and a pole-to-equator version
			List<Path.Command> tropicwardSegment = Path.reversed(polewardSegment);
			// remove one endpoint from each so there are no duplicate vertices
			polewardSegment = polewardSegment.subList(0, polewardSegment.size() - 1);
			tropicwardSegment = tropicwardSegment.subList(0, tropicwardSegment.size() - 1);
			// then build up the full shape by transforming the generic segments
			List<Path.Command> envelope = new ArrayList<>((edges[0].length + edges[1].length + 2)*polewardSegment.size());
			// go east to west in the north hemisphere
			for (int i = edges[0].length - 1; i >= 0; i --) {
				double westEdge = (i > 0) ? edges[0][i - 1] : -PI;
				double eastEdge = edges[0][i];
				envelope.addAll(Path.transformed(eastEdge - centers[0][i], 1, centers[0][i], 0, polewardSegment));
				envelope.addAll(Path.transformed(westEdge - centers[0][i], 1, centers[0][i], 0, tropicwardSegment));
			}
			// go west to east in the south hemisphere
			for (int i = 0; i < edges[1].length; i ++) {
				double westEdge = (i > 0) ? edges[1][i - 1] : -PI;
				double eastEdge = edges[1][i];
				envelope.addAll(Path.transformed(westEdge - centers[1][i], -1, centers[1][i], 0, polewardSegment));
				envelope.addAll(Path.transformed(eastEdge - centers[1][i], -1, centers[1][i], 0, tropicwardSegment));
			}
			// finally, convert it all to a Shape
			this.shape = Shape.polygon(Path.asArray(envelope));
		}

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
			"Eckert IV", "An equal-area projection released in a set of six (I'm only giving you the one because the others are not good)",
			new Shape(-2, 2, -1, 1, List.of(
					new Path.Command('M', -1, 1),
					new Path.Command('A', 1, 1, 0, 0, 1, -1, -1),
					new Path.Command('L', 1, -1),
					new Path.Command('A', 1, 1, 0, 0, 1, 1, 1),
					new Path.Command('Z'))),
			true, true, false, true, Type.PSEUDOCYLINDRICAL, Property.EQUAL_AREA, 3) {
		
		public double[] project(double lat, double lon) {
			double tht = NumericalAnalysis.newtonRaphsonApproximation(
					(2+PI/2)*sin(lat), lat,
					(t) -> (t + sin(2*t)/2 + 2*sin(t)),
					(t) -> (1 + cos(2*t) + 2*cos(t)), 1e-4);
			return new double[] { lon/PI*(1+cos(tht)), sin(tht)};
		}
		
		public double[] inverse(double x, double y) {
			double tht = asin(y);
			return new double[] {
					asin((tht + sin(2*tht)/2 + 2*sin(tht))/(2+PI/2)),
					x/(1 + cos(tht))*PI };
		}
		
	};
	
	
	public static final Projection WAGNER_II = new Projection(
			"Wagner II", "A compromise projection with sinusoidal meridians",
			null, true, true, true, true, Type.OTHER, Property.COMPROMISE, 2) {
		public void initialize(double... params) {
			this.shape = Shape.meridianEnvelope(this);
		}

		private final double c0 = 0.92483, c1 = 1.38725,
				c2 = 0.88022, c3 = 0.8855;
		
		public double[] project(double lat, double lon) {
			double psi = asin(c2*sin(c3*lat));
			return new double[] {c0*lon*cos(psi), c1*psi};
		}
		
		public double[] inverse(double x, double y) {
			double psi = y/c1;
			return new double[] {asin(sin(psi)/c2)/c3, x/c0/cos(psi)};
		}
	};
	
	
	public static final Projection WAGNER_V = new Projection(
			"Wagner V", "A compromise projection with elliptical meridians",
			null, true, true, true, true, Type.OTHER, Property.COMPROMISE, 3) {
		public void initialize(double... params) {
			this.shape = Shape.meridianEnvelope(this);
		}

		private final double c0 = 0.909771, c1 = 1.650142,
				c2 = 3.008957, c3 = 0.8855;
		
		public double[] project(double lat, double lon) {
			double psi = NumericalAnalysis.newtonRaphsonApproximation(
					c2*sin(c3*lat), sin(lat)*PI/3, (ps)->(2*ps + sin(2*ps)),
					(ps)->(2 + 2*cos(2*ps)), 1e-5);
			return new double[] {c0*lon*cos(psi), c1*sin(psi)};
		}
		
		public double[] inverse(double x, double y) {
			double psi = asin(y/c1);
			return new double[] {asin((2*psi + sin(2*psi))/c2)/c3, x/c0/cos(psi)};
		}
	};
	
	
	public static final Projection KAVRAYSKIY_VII = new Projection(
			"Kavrayskiy VII", "A compromise pseudocylindrical projection mostly popular in the former Soviet Union",
			null, true, true, true, true, Type.PSEUDOCYLINDRICAL, Property.COMPROMISE, 2) {
		public void initialize(double... params) {
			this.shape = Shape.meridianEnvelope(this);
		}

		public double[] project(double lat, double lon) {
			return new double[] { 1.5*lon*sqrt(1/3.-pow(lat/PI, 2)), lat };
		}
		
		public double[] inverse(double x, double y) {
			return new double[] { y, x/1.5/sqrt(1/3.-pow(y/PI, 2)) };
		}
		
	};
	
}
