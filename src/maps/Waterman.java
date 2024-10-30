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
import utils.Shape;

import static java.lang.Double.NaN;
import static java.lang.Math.PI;
import static java.lang.Math.hypot;
import static java.lang.Math.sqrt;
import static utils.Math2.linInterp;

/**
 * A class of methods pertaining to the Waterman projection.
 *
 * <a href="http://www.watermanpolyhedron.com/methodp.html">watermanpolyhedron.com/methodp</a>
 *
 * @author jkunimune
 */
public class Waterman {
	
	private static final double[] yELD =
		{-(sqrt(3)-1)/8, -sqrt(3)/8, -3*sqrt(3)/8, NaN};
	private static final double[] dXdL = {0, 1/(2*PI), 3/(2*PI), (2+sqrt(2))/(2*PI)};
	
	private static final double lonDiag = PI/(4+sqrt(8));
	private static final double sin15 = (sqrt(3)-1)/sqrt(8);  // sin(15°)
	private static final double cos15 = (sqrt(3)+1)/sqrt(8);  // cos(15°)
	
	
	public static final Projection FACE = new Projection(
			"Waterman (face)", "Steve Waterman", "A single face of Waterman's octahedral projection",
			Shape.polygon(new double[][] {{0., 0.}, {0., -sqrt(3)/2.}, {1/2., -sqrt(3)/2.}}),
			true, true, true, false, Type.OTHER, Property.COMPROMISE, 2) {
		
		public double[] project(double lat, double lon) {
			double[] xELD = jointPositions(lon);
			double desirLength = linInterp(lat, PI/2, 0, 0, totalLength(xELD, yELD));
			
			for (int i = 1; i < xELD.length; i ++) { //now that we've established the meridian, lets place our point on it
				double l = hypot(xELD[i]-xELD[i-1], yELD[i]-yELD[i-1]); //inch up the meridian
				if (l >= desirLength || i == xELD.length-1) //if it fits on this segment
					return new double[] {
							linInterp(desirLength, 0, l, xELD[i-1], xELD[i]),
							linInterp(desirLength, 0, l, yELD[i-1], yELD[i])}; //interpolate and return
				else //if we have further to go
					desirLength -= l; //mark off this length and continue
			}
			return null;
		}
		
		
		public double[] inverse(double x, double y) {
			if (x > -y-(sqrt(3)-1)/8 || x > -y/sqrt(3) || x > (2-sqrt(3))*(3/4.-y) ||
			    x > (7/4.+sqrt(3))+(2+sqrt(3))*y  || y < -sqrt(3)/2) //this describes the footprint of the octant
				return null;

			double longitude;
			int i = 0;
			while (i < yELD.length-1 && yELD[i] > y)
				i ++; //figure out in which segment it is
			if (i <= 2) { //well-behaved segments?
				longitude = x/linInterp(y, yELD[i-1], yELD[i], dXdL[i-1], dXdL[i]);
			}
			else { //the well-behaved part of the last segment?
				longitude = x/linInterp(y, yELD[i-1], -sqrt(3)/2, dXdL[i-1], dXdL[i]);
				if (longitude > lonDiag) { //the diagonal part of the last segment?
					double a = dXdL[2]*dXdL[3]*sin15; //surprisingly, the equation becomes quadratic here
					double b = (dXdL[3]*cos15-dXdL[2])*(y-yELD[2]) - dXdL[3]*sin15*x - dXdL[2]*(dXdL[3]*lonDiag*sin15-yELD[1]);
					double c = (dXdL[3]*lonDiag*sin15-yELD[1])*x + (1/4.-dXdL[3]*lonDiag*cos15)*(y-yELD[2]);
					longitude = (-b - sqrt(b*b - 4*a*c))/(2*a);
				}
			}
			
			double[] xELD = jointPositions(longitude);
			double totalLength = totalLength(xELD, yELD);
			double pointLength = hypot(x-xELD[i-1], y-yELD[i-1]);
			for (int j = 1; j < i; j ++) //add in the lengths of all previous segments
				pointLength += hypot(xELD[j]-xELD[j-1], yELD[j]-yELD[j-1]);
			double latitude = linInterp(pointLength, 0, totalLength, PI/2, 0);

			return new double[] {latitude, longitude};
		}
		
		
		private double[] jointPositions(double lon) { //POSTCONDITION: this also sets the value of yELD[3] based on the longitude. An odd place for a side-effect? Bad coding practice to manipulate static variables? meh. I documented it, didn't I?
			double[] xELD = new double[yELD.length]; //the positions of the intersections of this meridian with the Equal Line Delineations
			for (int i = 0; i < xELD.length; i ++)
				xELD[i] = dXdL[i]*lon;
			if (lon <= lonDiag) { //if it hits the equatorial ELD on the horizontal (hexagon) part
				yELD[3] = -sqrt(3)/2.;
			}
			else { //if it hits it on the diagonal (square) part
				xELD[3] = (lon-lonDiag)*dXdL[3]*cos15 + 1/4.;
				yELD[3] = (lon-lonDiag)*dXdL[3]*sin15 - sqrt(3)/2.;
			}
			return xELD;
		}
		
		
		private double totalLength(double[] xs, double[] ys) {
			double totalLength = 0; //compute meridian length
			for (int i = 1; i < xs.length; i ++)
				totalLength += hypot(xs[i]-xs[i-1], ys[i]-ys[i-1]);
			return totalLength;
		}
	};
	
}
