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

import utils.Math2;
import utils.NumericalAnalysis;

/**
 * A truncated octohedral map. The projection is Cahill-Keyes, because it was invented after Cahill,
 * so it is presumably better (not that Tobler World in a Square wasn't invented after Lambert EAC).
 * http://www.genekeyes.com/CKOG-OOo/7-CKOG-illus-&-coastline.html
 * 
 * @author jkunimune
 */
public class CahillKeyes {
	
	public static final double lMA = 940; //the distance from triangle vertex M to octant vertex A (the red length)
	public static final double lMG = 10000; //the altitude of the triangle
	private static final double lNG = lMG/Math.sqrt(3); //the height of the triangle
	private static final double lENy = lMA*Math.sqrt(3)/2; //the height difference between the triangle and the octant
	private static final double lMB = lMA*2/(Math.sqrt(3)-1); //the distance from triangle vertex M to octant vertex B (the blue length)
	private static final double lGF = lMG/Math.sqrt(3) - lMB; //the distance from triangle vertex G to octant vertex F
	private static final double lFE = lMA*Math.sqrt(2)/(Math.sqrt(3)-1); //the distance from octant vertex A to octant vertex B (the green length)
	private static final double lCV = 5760.8557; //the radius of arc CDV
	private static final double rC = 1560; //15 deg * 104 mm/deg = radius of zone c
	private static final double rE = 1760; //rC + 2 deg * 100 mm/deg = radius of zone e
	private static final double tF = 45*lGF/(lGF+lFE); //the longitude of point F
	private static final double dMEq = lGF/tF; //the number of mm per degree of longitude on the equator
	private static final double dMEqx = dMEq*(Math.sqrt(3)-1)/Math.sqrt(8); //dMEq*sind(15), for use on line segment EF
	private static final double dMEqy = dMEq*(Math.sqrt(3)+1)/Math.sqrt(8); //dMEq*cosd(15), for use on line segment EF
	private static final double bDE = (2*lNG-lMB)*(2-Math.sqrt(3)); //the y-intercept of the extension of line DE
	private static final double xU = lMA + Math.sqrt(3)/2*rE;
	private static final double xC = 2786.8887; //the x coordinate of the centre of arc CDV
	private static final double yU = rE/2;
	private static final double yC = 1609.0110; //the y coordinate of the centre of arc CDV
	
	private static final double TOLERANCE = 5; //this is a reasonable tolerance when you recall that we're operating on the order of 10,000 units


	/**
	 * convert adjusted lat and lon in degrees to Mary Jo's coordinates
	 */
	public static final double[] faceProjectD(double latD, double lonD) {
		final double[][] mer = meridian(lonD);
		if (latD >= 75) { //zone c (frigid zone)
			return new double[] { lMA + 104*(90-latD)*Math2.cosd(lonD),
					104*(90-latD)*Math2.sind(lonD) };
		}
		else if (latD >= 73 && lonD <= 30) { //zone e ()
			return new double[] { lMA + (rC+100*(75-latD))*Math2.cosd(lonD),
					(rC+100*(75-latD))*Math2.sind(lonD) };
		}
		else if (lonD <= 29) { //zone i (central zone)
			return meridianPoint(mer, Math2.linInterp(latD, 73, 0, rE, meridianLength(mer)));
		}
		else if (latD >= 73) { //zone j (frigid supple zone)
			final double distTU = meridianTUIntersect(lonD, mer);
			return meridianPoint(mer, Math2.linInterp(latD, 75, 73, rC, distTU));
		}
		else if (latD <= 15) { //zone k (torrid supple zone)
			final double distCDV = meridianCDVIntersect(mer);
			return meridianPoint(mer, Math2.linInterp(latD, 15, 0, distCDV, meridianLength(mer)));
		}
		else { //zone l (middle supple zone)
			final double distTU = meridianTUIntersect(lonD, mer);
			final double distCDV = meridianCDVIntersect(mer);
			return meridianPoint(mer, Math2.linInterp(latD, 73, 15, distTU, distCDV));
		}
	}
	
	
	public static final double[] faceInverseD(double x, double y) { //convert Mary Jo's coordinates to relative lat and lon in degrees
		if (y > x-lMA || y > x/Math.sqrt(3) || y > x*(2-Math.sqrt(3))+bDE ||
				y > (lMG-x)*(2+Math.sqrt(3))+lGF || x > lMG) //this describes the footprint of the octant
			return null;
		
		double lonD = longitudeD(x, y);
		double[][] mer = meridian(lonD);
		double len = meridianDistance(mer, x, y);
		
		if (len <= rC) { //zone c (frigid zone)
			return new double[] {90 - len/104, lonD};
		}
		else if (len <= rE && lonD < 30) { //zone e ()
			return new double[] {75 - (len-rC)/100, lonD};
		}
		else if (lonD <= 29) { //zone i (central zone)
			return new double[] {Math2.linInterp(len, rE, meridianLength(mer), 73, 0), lonD};
		}
		else if (y <= yU - Math.sqrt(3)*(x-xU)) { //zone j (frigid supple zone)
			final double distTU = meridianTUIntersect(lonD, mer);
			return new double[] {Math2.linInterp(len, rC, distTU, 75, 73), lonD};
		}
		else if (Math.hypot(x-xC, y-yC) >= lCV) { //zone k (torrid supple zone)
			final double distCDV = meridianCDVIntersect(mer);
			return new double[] {Math2.linInterp(len, distCDV, meridianLength(mer), 15, 0), lonD};
		}
		else { //zone l (middle supple zone)
			final double distTU = meridianTUIntersect(lonD, mer);
			final double distCDV = meridianCDVIntersect(mer);
			return new double[] {Math2.linInterp(len, distTU, distCDV, 73, 15), lonD};
		}
	}
	
	
	private static double[][] meridian(double lonD) { //calculate the locations of the endpoints and joints of this meridian
		if (lonD == 0)
			return new double[][] {{lMA, 0}, {lMG, 0}}; //the prime meridian is special
		
		final double m3 = lonD, m2 = lonD*2/3., m1 = lonD/3.;
		final double xE, yE;
		if (lonD <= tF) { //the equatorial joint is piecewise; it depends on if the meridian hits the vertical part,
			xE = lMG;
			yE = lonD*dMEq;
		}
		else { //or the slanted part
			xE = lMG - (lonD-tF)*dMEqx;
			yE = lGF + (lonD-tF)*dMEqy;
		}
		return new double[][] {
			{lMA, 0},
			{lMA/(1-Math2.tand(m2)/Math2.tand(m3)), lMA*Math2.tand(m3)*Math2.tand(m2)/(Math2.tand(m3)-Math2.tand(m2))},
			{(yE-xE*Math2.tand(m1))/(Math2.tand(m2)-Math2.tand(m1)), (yE*Math2.cotd(m1)-xE)/(Math2.cotd(m1)-Math2.cotd(m2))},
			{xE, yE}
		};
	}
	
	
	private static double meridianLength(double[][] meridian) { //calculate the length of a Cahill-Keyes meridian from pole to equator
		double l = 0;
		for (int i = 1; i < meridian.length; i ++)
			l += Math.hypot(meridian[i][0]-meridian[i-1][0], meridian[i][1]-meridian[i-1][1]);
		return l;
	}
	
	
	private static double[] meridianPoint(double[][] meridian, double dist) { //calculate a point on a meridian a given distance from the pole
		for (int i = 1; i < meridian.length; i ++) {
			double l = Math.hypot(meridian[i-1][0]-meridian[i][0], meridian[i-1][1]-meridian[i][1]);
			if (dist < l || i+1 == meridian.length)
				return new double[] {
						Math2.linInterp(dist, 0, l, meridian[i-1][0], meridian[i][0]),
						Math2.linInterp(dist, 0, l, meridian[i-1][1], meridian[i][1]) };
			else
				dist -= l;
		}
		return null;
	}
	
	
	private static double longitudeD(double x, double y) { //calculate the longitude of a given point
		double lonD0 = Math.toDegrees(Math.atan2(y, x-lMA)); //guess 0 for longitude
		double lonD1 = Math.toDegrees(Math.atan2(y, x))*1.5; //guess 1 for longitude
		if (lonD0 >= lonD1)
			return lonD0; //the point is north of the frigid joint; return guess 0
		double[][] mer1 = meridian(lonD1); //assuming it were on the temperate segment...
		if (x <= mer1[2][0])
			return lonD1; //the point is north of the temperate joint; return guess 1
		
		if (y <= lGF + (x-lMG)*Math2.tand(tF/3)) { //the point is on the torrid segment; we require iteration
			return NumericalAnalysis.newtonRaphsonApproximation(y, y/(lNG - lENy)*45, //does the meridian strike GF?
					(l) -> (dMEq*l + (x-lMG)*Math2.tand(l/3)),
					(l) -> (dMEq + (x-lMG)*Math.toRadians(Math.pow(Math2.secd(l/3), 2))/3),
					TOLERANCE);
		}
		else {
			return NumericalAnalysis.newtonRaphsonApproximation(y, y/(lNG - lENy)*45, //then it must strike FE!*
					(l) -> (dMEqy*(l-tF) + (x-lMG+dMEqx*(l-tF))*Math2.tand(l/3) + lGF),
					(l) -> (dMEqy + dMEqx*Math2.tand(l/3) + (x-lMG+dMEqx*(l-tF))*Math.toRadians(Math.pow(Math2.secd(l/3), 2))/3),
					TOLERANCE);
		}
	}
	
	
	private static double meridianDistance(double[][] meridian, double x, double y) { //how far up the meridian is this?
		double lenCum = 0; //the cumulative length of all previous meridians
		for (int i = 1; i < meridian.length; i ++) {
			double l = Math.hypot(meridian[i][0]-meridian[i-1][0], meridian[i][1]-meridian[i-1][1]);
			if (x <= meridian[i][0])
				return Math2.linInterp(x, meridian[i-1][0], meridian[i][0], lenCum, lenCum+l);
			else
				lenCum += l;
		}
		return lenCum; //if the point is outside the meridian, make a wild guess
	}
	
	
	private static final double meridianTUIntersect(double lonD, double[][] meridian) { //how far from the pole does the meridian intersect line TU
		if (meridian[1][1] > 2*rE - Math.sqrt(3)*(meridian[1][0]-lMA)) { //frigid joint is above TU; intersection is on the frigid segment
			return rE/Math2.cosd(30-lonD);
		}
		else { //frigid joint is below TU; intersection is on the temperate segment
			final double d = Math.hypot(xU-meridian[1][0], yU-meridian[1][1]);
			final double a = Math.atan2(rE/2-meridian[1][1], xU-meridian[1][0]);
			return Math.hypot(meridian[0][0] - meridian[1][0], meridian[0][1] - meridian[1][1])
					+ d * Math2.sind(60 + Math.toDegrees(a)) / Math2.sind(120 - 2/3.*lonD);
		}
	}
	
	
	private static final double meridianCDVIntersect(double[][] meridian) { //how far from the equator does the meridian intersect arc CDV
		final int i;
		if (Math.hypot(meridian[2][0]-xC, meridian[2][1]-yC) >= lCV) //torrid joint is outside CDV; intersection is on the temperate segment
			i = 2;
		else //torrid joint is inside CDV; intersection is on the torrid segment
			i = 3;
		final double m = (meridian[i][1]-meridian[i-1][1]) / (meridian[i][0]-meridian[i-1][0]);
		final double b = meridian[i][1] - m*meridian[i][0] - yC;
		final double l0 =
				Math.hypot(meridian[i][0]-meridian[i-1][0], meridian[i][1]-meridian[i-1][1]);
		final double x = ((xC-m*b) + Math.sqrt(Math.pow(m*b-xC,2) - (m*m+1)*(b*b+xC*xC-lCV*lCV))) / (m*m+1);
		final double l = Math2.linInterp(x, meridian[i-1][0], meridian[i][0], 0, l0);
		if (i == 2) //if it is on the temperate segment, add the length of the frigid segment
			return Math.hypot(meridian[1][0]-meridian[0][0], meridian[1][1]-meridian[0][1]) + l;
		else //otherwise, add the lengths of the frigid and temperate segments
			return Math.hypot(meridian[1][0] - meridian[0][0], meridian[1][1] - meridian[0][1])
					+ Math.hypot(meridian[2][0] - meridian[1][0], meridian[2][1] - meridian[1][1])
					+ l;
	}
} //*see bottom of Misc.java
