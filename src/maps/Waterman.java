package maps;
import java.util.Arrays;

import maps.Projection.Property;
import maps.Projection.Type;
import utils.Math2;

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
 * A class devoted to the Waterman butterfly, because it doesn't really fit into Polyhedral.java or
 * CahillKeyes.java, and Misc.java has too many projections already.
 * http://www.watermanpolyhedron.com/methodp.html
 * 
 * @author jkunimune
 */
public class Waterman {
	
	private static final double[] xELD =
		{(Math.sqrt(3)-1)/2, .5*Math.sqrt(3), 1.5*Math.sqrt(3), Double.NaN};
	private static final double[] dYdL = {0, 2/Math.PI, 6/Math.PI, (4+Math.sqrt(8))/Math.PI};
	
	private static final double lossX = (Math.sqrt(3)-1)/2; //horizontal dimension interruption loss
	private static final double lossY = lossX*Math.sqrt(3)/2; //vertical dimension interruption loss
	private static final double lonDiag = Math.PI/(4+Math.sqrt(8));
	private static final double sin15 = (Math.sqrt(3)-1)/Math.sqrt(8);
	private static final double cos15 = (Math.sqrt(3)+1)/Math.sqrt(8);
	
	
	public static final Projection BUTTERFLY =
			new Projection(
					"Waterman", "An aesthetically pleasing octohedral map arrangement.",
					8*Math.sqrt(3)-2*lossX, 8-2*lossY, 0b1010,
					Type.TETRADECAHEDRAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			double centralLon = Math.floor(lon/(Math.PI/2))*Math.PI/2 + Math.PI/4; //get the octant longitude
			if (lon == Math.PI) 	centralLon = 3*Math.PI/4; //a hack to fix a minor issue with the IDL
			
			double[] mjCoords = faceProject(Math.abs(lat), Math.abs(lon - centralLon));
			double mjX = mjCoords[0]; //"mj" stands for "Mary Jo Graca". I know she had nothing to do with the Waterman projection, but it means the same thing here as it does in Cahill-Keyes
			double mjY = mjCoords[1];
			
			if (lat < 0) //if it's in the southern hemisphere
				mjX = 4*Math.sqrt(3) - mjX; //flip it around
			if (lon < centralLon) //if the relative longitude is negative
				mjY = -mjY; //flip it the other way
			return new double[] {
					mjX*Math.sin(centralLon*2/3) + mjY*Math.cos(centralLon*2/3),
					-mjX*Math.cos(centralLon*2/3) + mjY*Math.sin(centralLon*2/3) + 2 };
		}
		
		public double[] inverse(double x, double y) {
			y -= 2;
			double tht = Math.atan2(x, -y);
			if (Math.abs(tht) > 5*Math.PI/6) 	return null; //only show a little bit of extra
			
			double quadrAngle = (Math.floor(tht/(Math.PI/3))+2)*Math.PI/3; //the angle of the centre of the quadrant, measured widdershins from -y
			double centralLon = quadrAngle*1.5 - 3*Math.PI/4; //the central meridian of this quadrant
			double mjX = -x*Math.cos(quadrAngle) - y*Math.sin(quadrAngle);
			double mjY =  x*Math.sin(quadrAngle) - y*Math.cos(quadrAngle);
			
			double[] relCoords = faceInverse(Math.min(mjX, 4*Math.sqrt(3)-mjX), Math.abs(mjY));
			if (relCoords == null)
				return null;
			
			if (mjY < 0) 	relCoords[1] *= -1; //the left half of the octant gets shifted west
			if (mjX > 2*Math.sqrt(3)) 	relCoords[0] *= -1; //the outer rim of the map is the southern hemisphere
			return new double[] { relCoords[0], relCoords[1] + centralLon };
		}
	};
	
	
	
	private static double[] faceProject(double lat, double lon) {
		double[] yELD = jointHeights(lon);
		double desirLength = Math2.linInterp(lat, Math.PI/2, 0, 0, totalLength(xELD, yELD));
		
		for (int i = 1; i < xELD.length; i ++) { //now that we've established the meridian, lets place our point on it
			double l = Math.hypot(xELD[i]-xELD[i-1], yELD[i]-yELD[i-1]); //inch up the meridian
			if (l >= desirLength || i == xELD.length-1) //if it fits on this segment
				return new double[] {
						Math2.linInterp(desirLength, 0, l, xELD[i-1], xELD[i]),
						Math2.linInterp(desirLength, 0, l, yELD[i-1], yELD[i]) }; //interpolate and return
			else //if we have further to go
				desirLength -= l; //mark off this length and continue
		}
		return null;
	}

	private static double[] faceInverse(double x, double y) {
		if (y > x-(Math.sqrt(3)-1)/2 || y > x/Math.sqrt(3) || y > (2-Math.sqrt(3))*(x+3) ||
				y > (7+4*Math.sqrt(3))-(2+Math.sqrt(3))*x  || x > 2*Math.sqrt(3)) //this describes the footprint of the octant
			return null;
		
		double longitude;
		int i = 0;
		while (i < xELD.length-1 && xELD[i] < x)
			i ++; //figure out in which segment it is
		if (i <= 2) { //well-behaved segments?
			longitude = y/Math2.linInterp(x, xELD[i-1], xELD[i], dYdL[i-1], dYdL[i]);
		}
		else { //the well-behaved part of the last segment?
			longitude = y/Math2.linInterp(x, xELD[i-1], 2*Math.sqrt(3), dYdL[i-1], dYdL[i]);
			if (longitude > lonDiag) { //the diagonal part of the last segment?
				double a = dYdL[2]*dYdL[3]*sin15; //surprisingly, the equation is quadratic here
				double b = (dYdL[3]*cos15-dYdL[2])*(1.5*Math.sqrt(3)-x) - dYdL[3]*sin15*y - Math.sqrt(3)/2*dYdL[2];
				double c = Math.sqrt(3)/2*y - x + 1.5*Math.sqrt(3);
				longitude = (-b - Math.sqrt(b*b - 4*a*c))/(2*a);
			}
		}
		
		double[] yELD = jointHeights(longitude);
		double totalLength = totalLength(xELD, yELD);
		double pointLength = Math.hypot(x-xELD[i-1], y-yELD[i-1]);
		for (int j = 1; j < i; j ++) //add in the lengths of all previous segments
			pointLength += Math.hypot(xELD[j]-xELD[j-1], yELD[j]-yELD[j-1]);
		double latitude = Math2.linInterp(pointLength, 0, totalLength, Math.PI/2, 0);
		
		return new double[] {latitude, longitude};
	}
	
	
	private static double[] jointHeights(double lon) { //POSTCONDITION: this also sets the value of xELD[3] based on the longitude. An odd place for a side-effect? Bad coding practice to manipulate static variables? meh. I documented it, didn't I?
		double[] yELD = new double[xELD.length]; //the height of the intersections of this meridian with the Equal Line Delineations
		for (int i = 0; i < xELD.length; i ++)
			yELD[i] = dYdL[i]*lon;
		if (Math.abs(lon) <= lonDiag) { //if it hits the equatorial ELD on the vertical (hexagon) part
			xELD[3] = 2*Math.sqrt(3);
		}
		else { //if it hits it on the diagonal (square) part
			xELD[3] = -(Math.abs(lon)-lonDiag)*dYdL[3]*sin15 + 2*Math.sqrt(3);
			yELD[3] =  (Math.abs(lon)-lonDiag)*dYdL[3]*cos15 + 1;
		}
		return yELD;
	}
	
	
	private static double totalLength(double[] xs, double[] ys) {
		double totalLength = 0; //compute meridian length
		for (int i = 1; i < xELD.length; i ++)
			totalLength += Math.hypot(xs[i]-xs[i-1], ys[i]-ys[i-1]);
		return totalLength;
	}
	
}
