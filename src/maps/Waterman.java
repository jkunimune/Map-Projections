package maps;

import static java.lang.Double.NaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.hypot;
import static java.lang.Math.sqrt;
import static utils.Math2.linInterp;

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

/**
 * A class of methods pertaining to the Waterman projection.
 * http://www.watermanpolyhedron.com/methodp.html
 * 
 * @author jkunimune
 */
public class Waterman {
	
	private static final double[] xELD =
		{(sqrt(3)-1)/2, .5*sqrt(3), 1.5*sqrt(3), NaN};
	private static final double[] dYdL = {0, 2/PI, 6/PI, (4+sqrt(8))/PI};
	
	private static final double lonDiag = PI/(4+sqrt(8));
	private static final double sin15 = (sqrt(3)-1)/sqrt(8);
	private static final double cos15 = (sqrt(3)+1)/sqrt(8);
	
	
	public static double[] faceProject(double lat, double lon) {
		double[] yELD = jointHeights(lon);
		double desirLength = linInterp(lat, PI/2, 0, 0, totalLength(xELD, yELD));
		
		for (int i = 1; i < xELD.length; i ++) { //now that we've established the meridian, lets place our point on it
			double l = hypot(xELD[i]-xELD[i-1], yELD[i]-yELD[i-1]); //inch up the meridian
			if (l >= desirLength || i == xELD.length-1) //if it fits on this segment
				return new double[] {
						linInterp(desirLength, 0, l, xELD[i-1], xELD[i]),
						linInterp(desirLength, 0, l, yELD[i-1], yELD[i]) }; //interpolate and return
			else //if we have further to go
				desirLength -= l; //mark off this length and continue
		}
		return null;
	}

	public static double[] faceInverse(double x, double y) {
		if (y > x-(sqrt(3)-1)/2 || y > x/sqrt(3) || y > (2-sqrt(3))*(x+3) ||
				y > (7+4*sqrt(3))-(2+sqrt(3))*x  || x > 2*sqrt(3)) //this describes the footprint of the octant
			return null;
		
		double longitude;
		int i = 0;
		while (i < xELD.length-1 && xELD[i] < x)
			i ++; //figure out in which segment it is
		if (i <= 2) { //well-behaved segments?
			longitude = y/linInterp(x, xELD[i-1], xELD[i], dYdL[i-1], dYdL[i]);
		}
		else { //the well-behaved part of the last segment?
			longitude = y/linInterp(x, xELD[i-1], 2*sqrt(3), dYdL[i-1], dYdL[i]);
			if (longitude > lonDiag) { //the diagonal part of the last segment?
				double a = dYdL[2]*dYdL[3]*sin15; //surprisingly, the equation becomes quadratic here
				double b = (dYdL[3]*cos15-dYdL[2])*(1.5*sqrt(3)-x) - dYdL[3]*sin15*y - dYdL[2]*(sqrt(3)/2+dYdL[3]*lonDiag*sin15);
				double c = (sqrt(3)/2+dYdL[3]*lonDiag*sin15)*y + (1-dYdL[3]*lonDiag*cos15)*(1.5*sqrt(3)-x);
				longitude = (-b - sqrt(b*b - 4*a*c))/(2*a);
			}
		}
		
		double[] yELD = jointHeights(longitude);
		double totalLength = totalLength(xELD, yELD);
		double pointLength = hypot(x-xELD[i-1], y-yELD[i-1]);
		for (int j = 1; j < i; j ++) //add in the lengths of all previous segments
			pointLength += hypot(xELD[j]-xELD[j-1], yELD[j]-yELD[j-1]);
		double latitude = linInterp(pointLength, 0, totalLength, PI/2, 0);
		
		return new double[] {latitude, longitude};
	}
	
	
	private static double[] jointHeights(double lon) { //POSTCONDITION: this also sets the value of xELD[3] based on the longitude. An odd place for a side-effect? Bad coding practice to manipulate static variables? meh. I documented it, didn't I?
		double[] yELD = new double[xELD.length]; //the height of the intersections of this meridian with the Equal Line Delineations
		for (int i = 0; i < xELD.length; i ++)
			yELD[i] = dYdL[i]*lon;
		if (abs(lon) <= lonDiag) { //if it hits the equatorial ELD on the vertical (hexagon) part
			xELD[3] = 2*sqrt(3);
		}
		else { //if it hits it on the diagonal (square) part
			xELD[3] = -(abs(lon)-lonDiag)*dYdL[3]*sin15 + 2*sqrt(3);
			yELD[3] =  (abs(lon)-lonDiag)*dYdL[3]*cos15 + 1;
		}
		return yELD;
	}
	
	
	private static double totalLength(double[] xs, double[] ys) {
		double totalLength = 0; //compute meridian length
		for (int i = 1; i < xELD.length; i ++)
			totalLength += hypot(xs[i]-xs[i-1], ys[i]-ys[i-1]);
		return totalLength;
	}
	
}
