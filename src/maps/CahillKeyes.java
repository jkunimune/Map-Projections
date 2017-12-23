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
 * A polyhedral map shaped like a butterfly. The projection is Cahill-Keyes, because it was the
 * best-documented one I found.
 * http://www.genekeyes.com/CKOG-OOo/7-CKOG-illus-&-coastline.html
 * 
 * @author jkunimune
 */
public class CahillKeyes {
	
	public static final Projection BUTTERFLY =
			new Projection(
					"Cahill-Keyes Butterfly", "An aesthetically pleasing octohedral map arrangement.",
					40000, 40000/Math.sqrt(3), 0b1110, Type.POLYHEDRAL, Property.COMPROMISE) {
		
		private final double OFFSET_Y = 10000/Math.sqrt(3);
		
		public double[] project(double lat, double lon) {
			final double centralLon = Math.floor(lon/(Math.PI/2))*Math.PI/2 + Math.PI/4; //get the octant longitude
			final double[] mjCoords = innerProject(
					Math.toDegrees(Math.abs(lat)), Math.toDegrees(Math.abs(lon - centralLon)));
			double mjX = mjCoords[0];
			double mjY = mjCoords[1];
			
			if (lat < 0) //if it's in the souther hemisphere
				mjX = 20000 - mjX; //flip it around
			if (lon < centralLon) //if the relative longitude is negative
				mjY = -mjY; //flip it the other way
			return new double[] {
					mjX*Math.sin(centralLon*2/3) + mjY*Math.cos(centralLon*2/3),
					-mjX*Math.cos(centralLon*2/3) + mjY*Math.sin(centralLon*2/3) + OFFSET_Y };
		}
		
		public double[] inverse(double x, double y) {
			return null; //TODO: projection wishlist
		}
	};
	
	
	public static final Projection M_MAP =
			new Projection(
					"Cahill-Keyes M-Profile", "A simple pleasing octohedral map arrangement.",
					40000, 30000/Math.sqrt(3), 0b1110, Type.POLYHEDRAL, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			final int octantNum = (int)Math.floor(lon/(Math.PI/2));
			final double centralLon = octantNum*Math.PI/2 + Math.PI/4; //get the octant longitude
			final double[] mjCoords = innerProject(
					Math.toDegrees(Math.abs(lat)), Math.toDegrees(Math.abs(lon - centralLon)));
			double mjX = 10000 - mjCoords[0];
			double mjY = mjCoords[1];
			
			if (lat < 0) //if it's in the souther hemisphere
				mjX = -mjX; //flip it around
			if (lon < centralLon) //if the relative longitude is negative
				mjY = -mjY; //flip it the other way
			double offsetX = centralLon*20000/Math.PI;
			double rotDirec = octantNum%2==0 ? 1. : -1.;
			return new double[] {
					rotDirec/2*mjX + Math.sqrt(3)/2*mjY + offsetX,
					Math.sqrt(3)/2*mjX - rotDirec/2*mjY };
		}
		
		public double[] inverse(double x, double y) {
			return null; //TODO: projection wishlist
		}
	};
	
	
	
	private static final double[] innerProject(double lat, double lon) { //convert adjusted lat and lon in degrees to Mary Jo's coordinates
		return new double[] { (90-lat)*111, (90-lat)*111*Math.tan(Math.toRadians(lon*2/3)) };
	}
	
	
	private static final double[] innerInverse(double x, double y) { //convert Mary Jo's coordinates to relative lat and lon in degrees
		return null;
	}
}
