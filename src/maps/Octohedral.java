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

import de.jtem.mfc.field.Complex;
import maps.Projection.Property;

/**
 * A class of maps that use octohedral octants. Very similar to Polyhedral, but much faster since
 * it takes advantage of the fact that everything is orthogonal.
 * 
 * @author jkunimune
 */
public class Octohedral {
	
	public static final Projection WATERMAN = new OctohedralProjection(
			"Waterman Butterfly", "A simple Cahill-esque octohedral map arrangement, with Antarctica left on.",
			2*Math.sqrt(3), (Math.sqrt(3)-1)/2, 0b1010, Property.COMPROMISE, 3,
			Configuration.BUTTERFLY) {
		
		protected double[] faceProject(double lat, double lon) {
			return Waterman.faceProject(lat, lon);
		}
		
		protected double[] faceInverse(double x, double y) {
			return Waterman.faceInverse(x, y);
		}
	};
	
	
	public static final Projection KEYES_BUTTERFLY = new OctohedralProjection(
			"Cahill-Keyes Butterfly", "A simple Cahill-esque octohedral map arrangement, with Antarctica left on.",
			CahillKeyes.lMG, CahillKeyes.lMA, 0b1010, Property.COMPROMISE, 4,
			Configuration.BUTTERFLY) {
		
		protected double[] faceProject(double lat, double lon) {
			return CahillKeyes.faceProjectD(Math.toDegrees(lat), Math.toDegrees(lon));
		}
		
		protected double[] faceInverse(double x, double y) {
			double[] coords = CahillKeyes.faceInverseD(x, y);
			return (coords == null) ? null :
				new double[] {Math.toRadians(coords[0]), Math.toRadians(coords[1])};
		}
	};
	
	
	public static final Projection KEYES_BASIC_M = new OctohedralProjection(
			"Cahill-Keyes Basic", "A simple M-shaped octohedral projection, with Antarctica broken into three pieces.",
			CahillKeyes.lMG, CahillKeyes.lMA, 0b1010, Property.COMPROMISE, 3,
			Configuration.M_PROFILE) {
		
		protected double[] faceProject(double lat, double lon) {
			return CahillKeyes.faceProjectD(Math.toDegrees(lat), Math.toDegrees(lon));
		}
		
		protected double[] faceInverse(double x, double y) {
			double[] coords = CahillKeyes.faceInverseD(x, y);
			return (coords == null) ? null :
				new double[] {Math.toRadians(coords[0]), Math.toRadians(coords[1])};
		}
	};
	
	
	public static final Projection KEYES_CONCIALDI = new OctohedralProjection(
			"Cahill-Keyes Bat", "I doubt I'll include this mashup in the final product; I just kind of want to see like what it looks.",
			CahillKeyes.lMG, CahillKeyes.lMA, 0b1010, Property.COMPROMISE, 4,
			Configuration.BAT_SHAPE) {
		
		protected double[] faceProject(double lat, double lon) {
			return CahillKeyes.faceProjectD(Math.toDegrees(lat), Math.toDegrees(lon));
		}
		
		protected double[] faceInverse(double x, double y) {
			double[] coords = CahillKeyes.faceInverseD(x, y);
			return (coords == null) ? null :
				new double[] {Math.toRadians(coords[0]), Math.toRadians(coords[1])};
		}
	};
	
	
	public static final Projection CAHILL_CONCIALDI = new OctohedralProjection(
			"Cahill-Concialdi Bat", "A conformal octohedral projection with no cuts and a unique arrangement.",
			Math.sqrt(3)/2, 0, 0b1000, Property.CONFORMAL, 4, Configuration.BAT_SHAPE) {
		
		private final double HEXAGON_SCALE = 1.112913; //this is 2^(2/3)/6*\int_0^\pi sin^(-1/3) x dx
		private final double TOLERANCE = 1e-3;
		private final double[] VERTEX = {0, Math.PI/4, -3*Math.PI/4};
		
		protected double[] faceProject(double lat, double lon) {
			double[] poleCoords = {lat, lon};
			double[] vertCoords = obliquifySphc(lat, lon, VERTEX); //look at an oblique aspect from the nearest vertex
			if (poleCoords[0] > vertCoords[0]) { //if this point is closer to the pole
				Complex w = Complex.fromPolar(Math.pow(Math.tan(Math.PI/4-lat/2), 2/3.), lon*2/3.);
				Complex z = polynomial(w); //project it as normal
				return new double[] {z.getRe(), z.getIm()};
			}
			else { //if it is closer to the vertex
				Complex w = Complex.fromPolar(
						Math.pow(Math.tan(Math.PI/4-vertCoords[0]/2), 2/3.), vertCoords[1]*2/3.);
				Complex zSkew = polynomial(w); //use the maclaurin series centred there
				return new double[] {
						-1/2.*zSkew.getRe() + Math.sqrt(3)/2*zSkew.getIm() + Math.sqrt(3)/2,
						-Math.sqrt(3)/2*zSkew.getRe() - 1/2.*zSkew.getIm() + 1/2. };
			}
		}
		
		protected double[] faceInverse(double x, double y) {
			Complex z;
			if (x < (1-y)/Math.sqrt(3)) //do the Newton Raphson from whichever vertex to which it is closest
				z = new Complex(x, y);
			else
				z = new Complex(-1/2.*(x-Math.sqrt(3)/2) - Math.sqrt(3)/2*(y-1/2.),
						Math.sqrt(3)/2*(x-Math.sqrt(3)/2) - 1/2.*(y-1/2.));
			Complex w = z.divide(HEXAGON_SCALE);
			Complex error = polynomial(w).minus(z);
			for (int i = 0; i < 8 && error.abs() > TOLERANCE; i ++) {
				Complex dzdw = derivative(w);
				w = w.minus(error.divide(dzdw));
				error = polynomial(w).minus(z);
			}
			double[] latLon = { Math.PI/2 - 2*Math.atan(Math.pow(w.abs(), 3/2.)), w.arg()*3/2. }; //inverse conic it back to spherical coordinates
			if (x < (1-y)/Math.sqrt(3)) //if it was closest to that vertex, the result is easy
				return latLon;
			else //if it was closer to the other vertex, do some obliquifying
				return obliquifyPlnr(latLon, VERTEX);
		}
		
		private Complex polynomial(Complex w) { //an approximation of the true conformal mapping function
			w = w.times(Complex.fromPolar(1, -Math.PI/6));
			Complex z = w.plus(w.pow(7).divide(21))
					.plus(w.pow(11).divide(99)).plus(w.pow(13).divide(1287/16.));
			return z.divide(Complex.fromPolar(HEXAGON_SCALE, -Math.PI/6));
		}
		
		private Complex derivative(Complex w) { //the derivative of polynomial()
			w = w.times(Complex.fromPolar(1, -Math.PI/6));
			Complex z = new Complex(1).plus(w.pow(6).divide(3))
					.plus(w.pow(10).divide(9)).plus(w.pow(12).divide(99/16.));
			return z.divide(Complex.fromPolar(HEXAGON_SCALE, -Math.PI/6));
		}
	};
	
	
	
	private static abstract class OctohedralProjection extends Projection {
		
		private final double size;
		private Configuration config;
		
		
		public OctohedralProjection(String name, String desc, double altitude, double cutSize,
				int fisc, Property property, int rating, Configuration config) {
			super(name, desc,
					config.fullWidth*altitude-config.cutWidth*cutSize,
					config.fullHeight*altitude-config.cutHeight*cutSize, fisc,
					(cutSize == 0) ? Type.OCTOHEDRAL : Type.TETRADECAHEDRAL, property, rating);
			this.size = altitude;
			this.config = config;
		}
		
		
		protected abstract double[] faceProject(double lat, double lon);
		
		protected abstract double[] faceInverse(double x, double y);
		
		
		public double[] project(double lat, double lon) {
			double[] octant = config.project(lat, lon); //octant properties
			double x0 = octant[0]*size, y0 = octant[1]*size, tht0 = octant[2], lon0 = octant[3];
			
			double[] coords = this.faceProject(Math.abs(lat), Math.abs(lon-lon0));
			double xMj = coords[0], yMj = coords[1]; //relative octant coordinates (Mj stands for "Mary Jo Graca")
			
			if (lat < 0) //reflect the southern hemisphere over the equator
				xMj = 2*size - xMj;
			if (lon-lon0 < 0)
				yMj = -yMj;
			
			return new double[] {
					x0 + Math.sin(tht0)*xMj + Math.cos(tht0)*yMj,
					y0 - Math.cos(tht0)*xMj + Math.sin(tht0)*yMj + config.fullHeight*size/2 };
		}
		
		
		public double[] inverse(double x, double y) {
			y = y - config.fullHeight*size/2; //measure from extrapolated top of map, not centre
			double[] octant = config.inverse(x/size, y/size);
			if (octant == null) 	return null;
			double x0 = size*octant[0], y0 = size*octant[1], tht0 = octant[2], lon0 = octant[3];
			
			double xMj = Math.sin(tht0)*(x-x0) - Math.cos(tht0)*(y-y0);
			double yMj = Math.cos(tht0)*(x-x0) + Math.sin(tht0)*(y-y0);
			
			if (Math.abs(yMj) > Math.min(xMj, 2*size-xMj)/Math.sqrt(3)) 	return null; //restrict to one rhombus
			double[] coords = this.faceInverse(Math.min(xMj, 2*size-xMj), Math.abs(yMj));
			if (coords == null) 	return null;
			double lat = coords[0], lon = coords[1];
			lat *= Math.signum(size-xMj);
			lon *= Math.signum(yMj);
			
			if (octant.length > 5 && (lat < octant[4] || lat > octant[5])) 	return null; //optional latitude limits
			return new double[] { lat, lon + lon0 };
		}
		
	}
	
	
	
	private enum Configuration {
		
		BUTTERFLY(4, 2, 4/Math.sqrt(3), Math.sqrt(3)) { //the classic four octants splayed out in a nice butterfly shape, with Antarctica divided and attached
			
			private final double Y_OFFSET = -1/Math.sqrt(3);
			
			public double[] project(double lat, double lon) {
				if (Math.abs(lon) >= Math.PI && lat < 0) {
					double sign = Math.signum(lon);
					return new double[] {sign, 2/Math.sqrt(3), sign*Math.PI/6, 5*sign*Math.PI/4};
				}
				double centralMerid = Math.floor(lon/(Math.PI/2))*Math.PI/2 + Math.PI/4;
				return new double[] { 0, Y_OFFSET, centralMerid*2/3., centralMerid };
			}
			
			public double[] inverse(double x, double y) {
				if (y > (1-Math.abs(x))/Math.sqrt(3)) {
					double sign = Math.signum(x);
					return new double[] { sign, 2/Math.sqrt(3), sign*Math.PI/6, sign*5*Math.PI/4 };
				}
				double tht = Math.atan2(x, -y+Y_OFFSET);
				if (Math.abs(tht) > 5*Math.PI/6)
					return null;
				double centralAngle = Math.floor(tht/(Math.PI/3))*Math.PI/3 + Math.PI/6;
				return new double[] { 0, Y_OFFSET, centralAngle, centralAngle*3/2. };
			}
		},
		
		
		M_PROFILE(4, 0, Math.sqrt(3), Math.sqrt(3)) { //The more compact zigzag configuration with Antarctica divided and attached
			
			public double[] project(double lat, double lon) {
				double centralMerid = Math.floor(lon/(Math.PI/2))*Math.PI/2 + Math.PI/4;
				double sign = Math.signum(centralMerid);
				return new double[] {
						sign, 0, sign*(Math.abs(centralMerid)*2/3.-Math.PI/3), centralMerid };
			}
			
			public double[] inverse(double x, double y) {
				double tht = Math.atan2(Math.abs(x)-1, -y);
				if (tht < -Math.PI/3) 	return null;
				double centralAngle = Math.floor(tht/(Math.PI/3))*Math.PI/3 + Math.PI/6;
				double sign = Math.signum(x);
				return new double[] {
						sign, 0, sign*centralAngle, sign*(centralAngle*3/2.+Math.PI/2) };
			}
		
		},
		
		
		M_W_S_POLE(4, 0, 3.5/Math.sqrt(3), 1.5*Math.sqrt(3)) { //Keyes's current configuration, with Antarctica reassembled in the center
			
			public double[] project(double lat, double lon) {
				// TODO: Implement this
				return null;
			}
			
			public double[] inverse(double x, double y) {
				// TODO: Implement this
				return null;
			}
		},
		
		
		BAT_SHAPE(2*Math.sqrt(3), 0, 2, 0) { //Luca Concialdi's obscure "Bat" arrangement that I liked. I don't think it's the best map possible as Luca does, but I do think it's quite neat
			
			public double[] project(double lat, double lon) {
				double centralMerid = Math.floor((lon+Math.PI/4)/(Math.PI/2))*Math.PI/2;
				if (Math.abs(centralMerid) == Math.PI && lat < 0) //the outer wings
					return new double[] { Math.signum(lon)*Math.sqrt(3), .5, 0, centralMerid };
				else if (centralMerid == 0 && lat < -Math.PI/4) //the bottoms of the wings
					return new double[] { 0, -2.5, Math.signum(lon)*2*Math.PI/3, centralMerid };
				else //the bulk of the map
					return new double[] { 0, -.5, centralMerid*2/3., centralMerid };
			}
			
			public double[] inverse(double x, double y) {
						double sign = Math.signum(x);
				if (y+.5 > Math.sqrt(3)*Math.abs(x))
					return null; //the empty top
				else if (y >= Math.sqrt(3)*(Math.sqrt(3)/2-Math.abs(x))) {
					if (y > -.5) 	return null; //more empty space
					return new double[] { sign*Math.sqrt(3), .5, 0, sign*Math.PI }; //the outer wings
				}
				else if (y <= -1.5 && Math.abs(x) > 1) {
					if (y+2.5 >= Math.abs(x)/Math.sqrt(3))
						return new double[] { 0, -2.5, Math.signum(x)*2*Math.PI/3, 0,
								-Math.PI/2, -Math.PI/4 }; //the bottoms of the wings
					else if (y+4.5 >= Math.abs(x)*Math.sqrt(3))
						return new double[] { 0, -2.5, sign*2*Math.PI/3, sign*2*Math.PI,
								-Math.PI/2, -Math.PI/4 }; //some more wing bottom
					else
						return new double[] { sign*Math.sqrt(3), -3.5, Math.PI, sign*3*Math.PI/2,
								-Math.PI/2, -Math.PI/4 }; //the bottoms of the bottoms of the wings
				}
				else {
					double tht = Math.atan2(x, -y-.5);
					double centralAngle = Math.floor((tht+Math.PI/6)/(Math.PI/3))*Math.PI/3;
					if (centralAngle != 0)
						return new double[] { 0, -.5, centralAngle, centralAngle*3/2.}; //the bulk of the map
					else
						return new double[] { 0, -.5, centralAngle, centralAngle*3/2.,
								-Math.PI/4, Math.PI/2 }; //the part with Antarctica cut off
				}
			}
		};
		
		
		public final double fullWidth, cutWidth, fullHeight, cutHeight;
		
		private Configuration(double fullWidth, double cutWidth,
				double fullHeight, double cutHeight) {
			this.fullWidth = fullWidth;
			this.cutWidth = cutWidth;
			this.fullHeight = fullHeight;
			this.cutHeight = cutHeight;
		}
		
		public abstract double[] project(double lat, double lon); //calculate the x, y, rotation, and central meridian for this quadrant
		public abstract double[] inverse(double x, double y); //calculate the x, y, rotation, central meridian, and min and max latitude for this quadrant
	}
	
}
