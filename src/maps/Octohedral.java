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
			Math.sqrt(3)/2, 0, 0b1000, Property.CONFORMAL, 3, Configuration.BAT_SHAPE) {
		
		protected double[] faceProject(double lat, double lon) {
			double[] poleCoords = {lat, lon};
			double[] vertCoords = obliquifySphc(lat, lon, new double[] {0, Math.PI/4, -3*Math.PI/4}); //look at an oblique aspect from the nearest vertex
			if (poleCoords[0] > vertCoords[0]) { //if this point is closer to the pole
				return fortyoctantProject(lat, lon); //project it as normal
			}
			else { //if it is closer to the vertex
				double[] skewCoords = fortyoctantProject(vertCoords[0], vertCoords[1]); //use the maclaurin series centered there
				return new double[] {
						-1/2.*skewCoords[0] + Math.sqrt(3)/2*skewCoords[1] + Math.sqrt(3)/2,
						-Math.sqrt(3)/2*skewCoords[0] - 1/2.*skewCoords[1] + 1/2. };
			}
		}
		
		private double[] fortyoctantProject(double lat, double lon) { //an even finer projection than the half-octant
			final Complex w = Complex.fromPolar(
					Math.pow(Math.tan(Math.PI/4-lat/2), 2/3.), lon*2/3.-Math.PI/6);
			Complex z = w.plus(w.pow(7).divide(21)).plus(w.pow(11).divide(99)).plus(w.pow(13).divide(1287/16.));
			if (z.isInfinite() || z.isNaN())	z = new Complex(0);
			z = z.times(new Complex(Math.sqrt(3)/2, 1/2.)).divide(1.112913); //this number is 2^(2/3)/6*\int_0^\pi sin^(-1/3) x dx
			return new double[] { z.getRe(), z.getIm() };
		}
		
		protected double[] faceInverse(double x, double y) {
			return null;
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
			double lon0 = octant[0], x0 = size*octant[1], y0 = size*octant[2], tht0 = octant[3];
			
			double xMj = Math.sin(tht0)*(x-x0) - Math.cos(tht0)*(y-y0);
			double yMj = Math.cos(tht0)*(x-x0) + Math.sin(tht0)*(y-y0);
			
			double[] coords = this.faceInverse(Math.min(xMj, 2*size-xMj), Math.abs(yMj));
			if (coords == null) 	return null;
			double lat = coords[0], lon = coords[1];
			
			return new double[] { Math.signum(size-xMj)*lat, Math.signum(yMj)*lon + lon0 };
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
					return new double[] { sign*5*Math.PI/4, sign, 2/Math.sqrt(3), sign*Math.PI/6 };
				}
				double tht = Math.atan2(x, -y+Y_OFFSET);
				if (Math.abs(tht) > 5*Math.PI/6)
					return null;
				double centralAngle = Math.floor(tht/(Math.PI/3))*Math.PI/3 + Math.PI/6;
				return new double[] { centralAngle*3/2., 0, Y_OFFSET, centralAngle };
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
						sign*(centralAngle*3/2.+Math.PI/2), sign, 0, sign*centralAngle };
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
				// TODO: Implement this
				return null;
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
		public abstract double[] inverse(double x, double y); //calculate the central meridian, x, y, and rotation for this quadrant
	}
	
}
