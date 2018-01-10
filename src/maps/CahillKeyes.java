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
import utils.Math2;

/**
 * A truncated octohedral map shaped like a butterfly. The projection is Cahill-Keyes, because it
 * was the best-documented one I found, and it's the newest one of which I know so it's presumably
 * best (not that Tobler World in a Square wasn't invented after Lambert EAC).
 * http://www.genekeyes.com/CKOG-OOo/7-CKOG-illus-&-coastline.html
 * 
 * @author jkunimune
 */
public class CahillKeyes {
	
	private static final double lMA = 940; //the distance from triangle vertex M to octant vertex A (the red length)
	private static final double lMG = 10000; //the altitude of the triangle
	private static final double lNG = lMG/Math.sqrt(3); //the height of the triangle
	private static final double lENy = lMA*(Math.sqrt(3)+1)/Math.sqrt(8); //the height difference between the triangle and the octant
	private static final double lMB = lMA*2/(Math.sqrt(3)-1); //the distance from triangle vertex M to octant vertex B (the blue length)
	private static final double lGF = lMG/Math.sqrt(3) - lMB; //the distance from triangle vertex G to octant vertex F
	private static final double lFE = lMA*Math.sqrt(2)/(Math.sqrt(3)-1); //the distance from octant vertex A to octant vertex B (the green length)
	private static final double lCV = 5760.8557; //the radius of arc CDV
	private static final double rC = 1560; //15 deg * 104 mm/deg = radius of zone c
	private static final double rE = 1760; //rC + 2 deg * 100 mm/deg = radius of zone e
	private static final double tF = 45*lGF/(lGF+lFE); //the longitude of point F
	private static final double dMEq = lGF/tF; //the number of mm per degree of longitude on the equator
	private static final double xU = lMA + Math.sqrt(3)/2*rE;
	private static final double xC = 2786.8887; //the x coordinate of the centre of arc CDV
	private static final double yC = 1609.0110; //the y coordinate of the centre of arc CDV
	
	
	public static final Projection BUTTERFLY =
			new Projection(
					"Cahill-Keyes Butterfly", "An aesthetically pleasing octohedral map arrangement.",
					4*lMG-2*lMA, 4*lNG-2*lENy, 0b1110, Type.POLYHEDRAL, Property.COMPROMISE) {
		
		private final double OFFSET_Y = 10000/Math.sqrt(3);
		
		public double[] project(double lat, double lon) {
			final double centralLon = Math.floor(lon/(Math.PI/2))*Math.PI/2 + Math.PI/4; //get the octant longitude
			final double[] mjCoords = innerProject(
					Math.toDegrees(Math.abs(lat)), Math.toDegrees(Math.abs(lon - centralLon)));
			double mjX = mjCoords[0];
			double mjY = mjCoords[1];
			
			if (lat < 0) //if it's in the souther hemisphere
				mjX = 2*lMG - mjX; //flip it around
			if (lon < centralLon) //if the relative longitude is negative
				mjY = -mjY; //flip it the other way
			return new double[] {
					mjX*Math.sin(centralLon*2/3) + mjY*Math.cos(centralLon*2/3),
					-mjX*Math.cos(centralLon*2/3) + mjY*Math.sin(centralLon*2/3) + OFFSET_Y };
		}
		
		public double[] inverse(double x, double y) {
			y -= OFFSET_Y;
			double quadrAngle = (Math.floor(Math.atan2(x, -y)/(Math.PI/3))+2)*Math.PI/3; //the angle of the centre of the quadrant, measured widdershins from -y
			double centralLon = quadrAngle*1.5 - 3*Math.PI/4; //the central meridian of this quadrant
			double mjX = -x*Math.cos(quadrAngle) - y*Math.sin(quadrAngle);
			double mjY =  x*Math.sin(quadrAngle) - y*Math.cos(quadrAngle);
			
			double[] relCoords = innerInverse(Math.min(mjX, 2*lMG-mjX), Math.abs(mjY));
			if (relCoords == null)
				return null;
			
			if (mjY < 0) 	relCoords[1] *= -1; //the left half of the octant gets shifted west
			if (mjX > lMG) 	relCoords[0] *= -1; //the outer rim of the map is the southern hemisphere
			return new double[] { Math.toRadians(relCoords[0]),
					Math.toRadians(relCoords[1]) + centralLon };
		}
	};
	
	
	public static final Projection M_MAP =
			new Projection(
					"Cahill-Keyes M-Profile", "A simple pleasing octohedral map arrangement.",
					4*lMG, 3*lNG-2*lENy, 0b1110, Type.POLYHEDRAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			final int octantNum = (int)Math.floor(lon/(Math.PI/2));
			final double centralLon = octantNum*Math.PI/2 + Math.PI/4; //get the octant longitude
			final double[] mjCoords = innerProject(
					Math.toDegrees(Math.abs(lat)), Math.toDegrees(Math.abs(lon - centralLon)));
			double mjX = 10000 - mjCoords[0];
			double mjY = mjCoords[1];
			
			if (lat < 0) //if it's in the southern hemisphere
				mjX = -mjX; //flip it around
			if (lon < centralLon) //if the relative longitude is negative
				mjY = -mjY; //flip it the other way
			double offsetX = centralLon*20000/Math.PI;
			double rotDirec = (octantNum%2 == 0) ? 1. : -1.;
			double sinRot = Math.sqrt(0.75);
			double cosRot = 0.5;
			return new double[] {
					rotDirec*cosRot*mjX + sinRot*mjY + offsetX,
					sinRot*mjX - rotDirec*cosRot*mjY };
		}
		
		public double[] inverse(double x, double y) {
			int quadrantNum = (int)Math.floor(x/lMG) + 2;
			double offsetX = (quadrantNum - 1.5)*lMG; //quadrant center coordinates
			double rotDirec = (quadrantNum%2 == 0) ? 1. : -1.;
			double sinRot = Math.sqrt(0.75);
			double cosRot = 0.5;
			double diagX = rotDirec*cosRot*(x-offsetX) + sinRot*y;
			double diagY = sinRot*(x-offsetX) - rotDirec*cosRot*y;
			
			double[] relCoords = innerInverse(lMG - Math.abs(diagX), Math.abs(diagY));
			if (relCoords == null)
				return null;
			
			if (diagY < 0) 	relCoords[1] *= -1; //the left half of the octant gets shifted west
			if (diagX < 0) 	relCoords[0] *= -1; //the bottom half of the map is the southern hemisphere
			return new double[] {Math.toRadians(relCoords[0]),
					Math.toRadians(relCoords[1] + (quadrantNum-1.5)*90) };
		}
	};
	
	
	public static final Projection OCTANT = 
			new Projection(
					"Cahill-Keyes Octant", "A singular octant from a Cahill-Keyes projection.",
					lMG, lNG, 0b0111, Type.POLYHEDRAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			lat = Math.max(0, lat);
			lon = Math.max(0, Math.min(Math.PI/4, lon));
			final double[] mjCoords = innerProject(
					Math2.round(Math.toDegrees(Math.abs(lat)),6), //I apply the round so that graticule lines exactly on the thing get done properly
					Math2.round(Math.toDegrees(Math.abs(lon)),6));
			return new double[] { mjCoords[0] - width/2, mjCoords[1] - height/2 };
		}
		
		public double[] inverse(double x, double y) {
			double[] coordsDeg = innerInverse(x+width/2, y+height/2);
			return new double[] { Math.toRadians(coordsDeg[0]), Math.toRadians(coordsDeg[1]) };
		}
	};
	
	
	
	private static final double[] innerProject(double lat, double lon) { //convert adjusted lat and lon in degrees to Mary Jo's coordinates
		final double[][] mer = meridian(lon);
		if (lat >= 75) { //zone c (frigid zone)
			return new double[] {
					lMA + 104*(90-lat)*Math2.cosd(lon), 104*(90-lat)*Math2.sind(lon) };
		}
		else if (lat >= 73 && lon <= 30) { //zone e
			return new double[] {
					lMA + (rC+100*(75-lat))*Math2.cosd(lon), (rC+100*(75-lat))*Math2.sind(lon) };
		}
		else if (lon <= 29) { //zone i (central zone)
			return meridianPoint(mer, Math2.linInterp(lat, 73, 0, rE, meridianLength(mer)));
		}
		else if (lat >= 73) { //zone j (frigid supple zone)
			final double distTU = meridianTUIntersect(lon, mer);
			return meridianPoint(mer, Math2.linInterp(lat, 75, 73, rC, distTU));
		}
		else if (lat <= 15) { //zone k (equatorial supple zone)
			final double distCDV = meridianCDVIntersect(mer);
			return meridianPoint(mer, Math2.linInterp(lat, 15, 0, distCDV, meridianLength(mer)));
		}
		else { //zone l (middle supple zone)
			final double distTU = meridianTUIntersect(lon, mer);
			final double distCDV = meridianCDVIntersect(mer);
			return meridianPoint(mer, Math2.linInterp(lat, 73, 15, distTU, distCDV));
		}
	}
	
	
	private static final double[] innerInverse(double x, double y) { //convert Mary Jo's coordinates to relative lat and lon in degrees
		if (Math.atan2(y, x) > Math.PI/6) return null;
		return new double[] { 90 - Math.hypot(x, y)/110, Math.toDegrees(Math.atan2(y, x)*1.5) };
	}
	
	
	private static double[][] meridian(double lon) { //calculate the locations of the endpoints and joints of this meridian
		if (lon == 0)
			return new double[][] {{lMA, 0}, /*{lMA+lFE, 0}, {lMG-lFE, 0}, */{lMG, 0}}; //I'll add more vertices here if I must
		final double m3 = lon, m2 = lon*2/3., m1 = lon/3.;
		final double xE, yE;
		if (lon <= tF) { //the equatorial joint is piecewise; it depends on if the meridian hits the vertical part,
			xE = lMG;
			yE = lon*dMEq;
		}
		else { //or the slanted supple zone part
			xE = lMG - (lon-tF)*dMEq*(Math.sqrt(3)-1)/Math.sqrt(8);
			yE = lGF + (lon-tF)*dMEq*(Math.sqrt(3)+1)/Math.sqrt(8);
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
	
	
	private static final double meridianTUIntersect(double lon, double[][] meridian) { //how far from the pole does the meridian intersect line TU
		if (meridian[1][1] > 2*rE - Math.sqrt(3)*(meridian[1][0]-lMA)) { //frigid joint is above TU; intersection is on the frigid segment
			return rE/Math2.cosd(30-lon);
		}
		else { //frigid joint is below TU; intersection is on the temperate segment
			final double d = Math.hypot(xU-meridian[1][0], rE/2-meridian[1][1]);
			final double a = Math.atan2(rE/2-meridian[1][1], xU-meridian[1][0]);
			return Math.hypot(meridian[0][0] - meridian[1][0], meridian[0][1] - meridian[1][1])
					+ d * Math2.sind(60 + Math.toDegrees(a)) / Math2.sind(120 - 2/3.*lon);
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
}
