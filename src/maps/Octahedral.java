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
import image.Path;
import maps.Projection.Property;
import utils.Shape;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.NaN;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.toRadians;
import static utils.Math2.cosd;
import static utils.Math2.cotd;
import static utils.Math2.floorMod;
import static utils.Math2.max;
import static utils.Math2.sind;

/**
 * A class of maps that use octahedral octants. Very similar to Polyhedral, but much faster since
 * it takes advantage of the fact that everything is orthogonal.
 * 
 * @author jkunimune
 */
public class Octahedral {
	
	public static final Projection CONFORMAL_CAHILL_FACE = new Projection(
			"Cahill Conformal (face)", "The conformal projection from an octant to an equilateral triangle",
			Shape.polygon(new double[][] {{0., 0.}, {0., -sqrt(3)/2.}, {1/2., -sqrt(3)/2.}}),
			true, true, false, false, Projection.Type.OCTAHEDRAL, Property.CONFORMAL, 3) {
		
		private final double HEXAGON_SCALE = 1.112913; //this is 2^(2/3)/6*\int_0^\pi sin^(-1/3) x dx
		private final double TOLERANCE = 1e-3;
		private final double[] VERTEX = {0, PI/4, -3*PI/4};
		
		public double[] project(double lat, double lon) {
			double[] poleCoords = {lat, lon};
			double[] vertCoords = transformFromOblique(lat, lon, VERTEX); //look at an oblique aspect from the nearest vertex
			if (poleCoords[0] > vertCoords[0]) { //if this point is closer to the pole
				Complex w = Complex.fromPolar(pow(tan(PI/4-lat/2), 2/3.), lon*2/3.);
				Complex z = polynomial(w); //project it as normal
				return new double[] {z.getIm(), -z.getRe()}; //rotate 90°
			}
			else { //if it is closer to the vertex
				Complex w = Complex.fromPolar(
						pow(tan(PI/4-vertCoords[0]/2), 2/3.), vertCoords[1]*2/3.);
				Complex zSkew = polynomial(w); //use the maclaurin series centred there
				return new double[] { //rotate and translate in the plane appropriately
						-sqrt(3)/2*zSkew.getRe() - 1/2.*zSkew.getIm() + 1/2.,
						 1/2.*zSkew.getRe() - sqrt(3)/2*zSkew.getIm() - sqrt(3)/2 };
			}
		}
		
		public double[] inverse(double x, double y) {
			Complex z;
			if (y > (x-1)/sqrt(3)) //do the Newton Raphson from whichever vertex to which it is closest
				z = new Complex(-y, x); //applying a transformation in the plane as appropriate
			else
				z = new Complex(-sqrt(3)/2*(x-1/2.) + 1/2.*(y+sqrt(3)/2),
				                -1/2.*(x-1/2.) - sqrt(3)/2*(y+sqrt(3)/2));
			Complex w = z.divide(HEXAGON_SCALE);
			Complex error = polynomial(w).minus(z);
			for (int i = 0; i < 8 && error.abs() > TOLERANCE; i ++) {
				Complex dzdw = derivative(w);
				w = w.minus(error.divide(dzdw));
				error = polynomial(w).minus(z);
			}
			double[] latLon = { PI/2 - 2*atan(pow(w.abs(), 3/2.)), w.arg()*3/2. }; //inverse conic it back to spherical coordinates
			if (y > (x-1)/sqrt(3)) //if it was closest to that vertex, the result is easy
				return latLon;
			else //if it was closer to the other vertex, do some obliquifying
				return transformToOblique(latLon, VERTEX);
		}
		
		private Complex polynomial(Complex w) { //an approximation of the true conformal mapping function
			w = w.times(Complex.fromPolar(1, -PI/6));
			Complex z = w.plus(w.pow(7).divide(21))
					.plus(w.pow(11).divide(99)).plus(w.pow(13).divide(1287/16.));
			return z.divide(Complex.fromPolar(HEXAGON_SCALE, -PI/6));
		}
		
		private Complex derivative(Complex w) { //the derivative of polynomial()
			w = w.times(Complex.fromPolar(1, -PI/6));
			Complex z = new Complex(1).plus(w.pow(6).divide(3))
					.plus(w.pow(10).divide(9)).plus(w.pow(12).divide(99/16.));
			return z.divide(Complex.fromPolar(HEXAGON_SCALE, -PI/6));
		}
	};
	
	
	public static final OctahedralProjection CONFORMAL_CAHILL_BUTTERFLY = new OctahedralProjection(
			"Cahill Conformal", "The conformal and only reproducible variant of Cahill's original map",
			0, Property.CONFORMAL, 3, CONFORMAL_CAHILL_FACE, Configuration.BUTTERFLY);
	
	
	public static final Projection CAHILL_CONCIALDI = new OctahedralProjection(
			"Cahill\u2013Concialdi", "A conformal octahedral projection with a unique arrangement",
			0, Property.CONFORMAL, 4, CONFORMAL_CAHILL_FACE, Configuration.BAT_SHAPE);
	
	
	public static final OctahedralProjection CONFORMAL_CAHILL_OCTANT = new OctahedralProjection(
			"Cahill Conformal (single octant)", "A single octant of Cahill's conformal projection (for memory economization in the case of very large maps)",
			0, Property.CONFORMAL, 3, CONFORMAL_CAHILL_FACE, Configuration.SINGLE_OCTANT);
	
	
	public static final Projection WATERMAN_BUTTERFLY = new OctahedralProjection(
			"Waterman Butterfly", "A Cahill-esque octahedral map arrangement, with Antarctica broken off and reassembled",
			(sqrt(3)-1)/8, Property.COMPROMISE, 3,
			Waterman.FACE, Configuration.BUTTERFLY_WITH_SOUTH_POLE);
	
	
	public static final Projection WATERMAN_SIMPLIFIED = new OctahedralProjection(
			"Waterman Butterfly (simplified)", "A simple Cahill-esque octahedral map arrangement, with Antarctica broken into four pieces",
			(sqrt(3)-1)/8, Property.COMPROMISE, 3,  // TODO: stating compromise here is unnecessary
			Waterman.FACE, Configuration.BUTTERFLY);
	
	
	public static final Projection WATERMAN_OCTANT = new OctahedralProjection(
			"Waterman (single octant)", "A single octant of Waterman's octahedral projection (for memory economization in the case of very large maps)",
			(sqrt(3)-1)/8, Property.COMPROMISE, 3,
			Waterman.FACE, Configuration.SINGLE_OCTANT);
	
	
	public static final Projection KEYES_STANDARD = new OctahedralProjection(
			"Cahill\u2013Keyes", "An M-shaped octahedral projection with Antarctica assembled in the center",
			CahillKeyes.POLE_OFFSET, Property.COMPROMISE, 4,
			CahillKeyes.FACE, Configuration.M_PROFILE_WITH_SOUTH_POLE);
	
	
	public static final Projection KEYES_SIMPLIFIED = new OctahedralProjection(
			"Cahill\u2013Keyes (simplified)", "A simple M-shaped octahedral projection, with Antarctica broken into three pieces",
			CahillKeyes.POLE_OFFSET, Property.COMPROMISE, 3,
			CahillKeyes.FACE, Configuration.M_PROFILE);
	
	
	public static final Projection KEYES_BUTTERFLY = new OctahedralProjection(
			"Cahill\u2013Keyes (butterfly)", "Keyes's octahedral projection with Cahill's original butterfly arrangement",
			CahillKeyes.POLE_OFFSET, Property.COMPROMISE, 3,
			CahillKeyes.FACE, Configuration.BUTTERFLY);
	
	
	public static final Projection KEYES_OCTANT = new OctahedralProjection(
			"Cahill\u2013Keyes (single octant)", "A single octant of the Cahill\u2013Keyes projection (for memory economization in the case of very large maps)",
			CahillKeyes.POLE_OFFSET, Property.COMPROMISE, 3,
			CahillKeyes.FACE, Configuration.SINGLE_OCTANT);
	
	
	
	private static class OctahedralProjection extends Projection {
		
		private final Projection faceProj;
		private final Octant[] octants;
		
		
		public OctahedralProjection(String name, String desc, double tipOffset,
		                            Property property, int rating,
		                            Projection faceProj, Configuration config) {
			super(name, desc, null, false, config.finite, faceProj.isSolveable(), faceProj.isInvertable(),
			      (tipOffset == 0) ? Type.OCTAHEDRAL : Type.TETRADECAHEDRAL, property, rating,
			      new String[] {}, new double[][] {}, config.hasAspect);
			this.octants = config.placeOctants(tipOffset);
			this.faceProj = faceProj;
			this.shape = config.drawShape(tipOffset, faceProj);
		}
		
		
		public double[] project(double lat, double lon) {
			for (Octant octant: this.octants) { // try each octant
				double lonr = floorMod(lon - octant.centralLongitude + PI, 2*PI) - PI; // relative longitude
				if (abs(lonr) > PI/4 || lat < octant.minLatitude || lat > octant.maxLatitude) // if it doesn't fit...
					continue; // check the next one
				if (lonr < octant.minLongitude || lonr > octant.maxLongitude) // also check the longitude restriction
					continue;
				double th = octant.planeRotation; // rotation angle
				
				double[] coords = this.faceProj.project(abs(lat), abs(lonr));
				double xMj = coords[0], yMj = coords[1]; //relative octant coordinates (Mj stands for "Mary Jo Graca")
				
				if (lat < 0) //reflect the southern hemisphere over the equator
					yMj = -sqrt(3) - yMj;
				if (lonr < 0) //and reflect the western hemisphere over the central meridian
					xMj = -xMj;
				
				return new double[] {
						octant.x + cos(th)*xMj - sin(th)*yMj,
						octant.y + sin(th)*xMj + cos(th)*yMj };
			}
			return new double[] {NaN, NaN}; // if none of the octants fit, return null
		}
		
		
		public double[] inverse(double x, double y) {
			for (Octant octant: this.octants) { // try each octant
				double th = octant.planeRotation; // rotation angle
				
				double xMj =  cos(th)*(x - octant.x) + sin(th)*(y - octant.y); // do the coordinate change
				double yMj = -sin(th)*(x - octant.x) + cos(th)*(y - octant.y);
				if (abs(xMj) > -max(yMj, -sqrt(3) - yMj)/sqrt(3) + 1e-12) // if the angle is wrong,
					continue; // check the next one
				
				double[] coords = this.faceProj.inverse(abs(xMj), max(yMj, -sqrt(3) - yMj));
				if (coords == null)
					continue; // if you got nothing, keep looking
				double lat = coords[0], lon = coords[1]; // project
				
				lat *= signum(yMj + sqrt(3)/2); // undo the reflections
				lon *= signum(xMj);
				if (lat < octant.minLatitude - 1e-6 || lat > octant.maxLatitude + 1e-6) // if the resulting coordinates are wrong
					continue; // move on
				if (lon < octant.minLongitude - 1e-6 || lon > octant.maxLongitude + 1e-6) // also check the longitude restriction
					continue;
				
				lon = floorMod(lon + octant.centralLongitude + PI, 2*PI)- PI;
				
				return new double[] { lat, lon };
			}
			return null;
		}
		
	}
	
	
	
	private enum Configuration {
		
		/** the classic four quadrants splayed out in a nice butterfly shape, with Antarctica divided and attached */
		BUTTERFLY(true, true) {
			Octant[] placeOctants(double tipOffset) {
				return new Octant[] {
						new Octant(0, 0, -PI/2, -PI/2, PI/2, -3*PI/4),
						new Octant(0, 0, -PI/6, -PI/2, PI/2,   -PI/4),
						new Octant(0, 0,  PI/6, -PI/2, PI/2,    PI/4),
						new Octant(0, 0,  PI/2, -PI/2, PI/2,  3*PI/4),
						};
			}
			Shape drawShape(double tipOffset, Projection projection) {
				List<Path.Command> shape = new ArrayList<>();
				shape.addAll(starShape( 0.0        ,  0.0,  2*PI/3, -2*PI/3, 4, tipOffset));
				shape.addAll(starShape(-0.5*sqrt(3),  0.5,    PI/3,   -PI/3, 2, tipOffset));
				shape.addAll(starShape(-1.0*sqrt(3),  0.0,  2*PI/3,    PI/3, 1, tipOffset));
				shape.addAll(starShape(-0.5*sqrt(3), -0.5,  4*PI/3,     0  , 4, tipOffset));
				shape.addAll(starShape(-0.5*sqrt(3), -1.5,    PI  ,  2*PI/3, 1, tipOffset));
				shape.addAll(starShape( 0.0        , -1.0,  5*PI/3,    PI/3, 4, tipOffset));
				shape.addAll(starShape( 0.5*sqrt(3), -1.5, -2*PI/3,   -PI  , 1, tipOffset));
				shape.addAll(starShape( 0.5*sqrt(3), -0.5,     0  , -4*PI/3, 4, tipOffset));
				shape.addAll(starShape( 1.0*sqrt(3),  0.0,   -PI/3, -2*PI/3, 1, tipOffset));
				shape.addAll(starShape( 0.5*sqrt(3),  0.5,    PI/3,   -PI/3, 2, tipOffset));
				return Shape.polygon(Path.asArray(shape));
			}
		},
		
		/** the classic four quadrants splayed out in a nice butterfly shape, with Antarctica divided and attached */
		BUTTERFLY_WITH_SOUTH_POLE(true, false) {
			private final double ySouthPole = -1.32;
			Octant[] placeOctants(double tipOffset) {
				double quadrantLength = sqrt(3) - tipOffset;
				return new Octant[] {
						new Octant(0, 0, -PI/2, toRadians(-63), toRadians( 90), toRadians(-155)),
						new Octant(0, 0, -PI/6, toRadians(-63), toRadians( 90), toRadians( -65)),
						new Octant(0, 0,  PI/6, toRadians(-63), toRadians( 90), toRadians(  25)),
						new Octant(0, 0,  PI/2, toRadians(-63), toRadians( 90), toRadians( 115)),
						new Octant(-quadrantLength/sqrt(2), ySouthPole - quadrantLength/sqrt(2),  3*PI/4,
						           toRadians(-90), toRadians(-63), toRadians(-155)),
						new Octant(-quadrantLength/sqrt(2), ySouthPole + quadrantLength/sqrt(2),    PI/4,
						           toRadians(-90), toRadians(-63), toRadians( -65)),
						new Octant( quadrantLength/sqrt(2), ySouthPole + quadrantLength/sqrt(2),   -PI/4,
						           toRadians(-90), toRadians(-63), toRadians(  25)),
						new Octant( quadrantLength/sqrt(2), ySouthPole - quadrantLength/sqrt(2), -3*PI/4,
						           toRadians(-90), toRadians(-63), toRadians( 115)),
						};
			}
			Shape drawShape(double tipOffset, Projection projection) {
				double shortEdge = tipOffset/(2*sind(15));
				List<Path.Command> parallelSegment = projection.drawLoxodrome(
						toRadians(63), toRadians(45), toRadians(63), toRadians(0), .1);
				parallelSegment.addAll(Path.scaled(-1, 1, Path.reversed(parallelSegment)));  // expand the parallel segment to cover 90°
				// compose the shape of the main butterfly body
				List<Path.Command> mainShape = new ArrayList<>();
				mainShape.addAll(starShape( 0.0        ,  0.0,  2*PI/3, -2*PI/3, 4, tipOffset));
				mainShape.addAll(starShape(-0.5*sqrt(3),  0.5,    PI/3,   -PI/3, 2, tipOffset));
				mainShape.addAll(Path.rotated(PI/2, Path.translated(0, sqrt(3), parallelSegment)));
				mainShape.addAll(starShape(-0.5*sqrt(3), -0.5,  4*PI/3,     0  , 4, tipOffset));
				mainShape.addAll(Path.rotated(5*PI/6, Path.translated(0, sqrt(3), parallelSegment)));
				mainShape.addAll(starShape( 0.0        , -1.0,  5*PI/3,    PI/3, 4, tipOffset));
				mainShape.addAll(Path.rotated(-5*PI/6, Path.translated(0, sqrt(3), parallelSegment)));
				mainShape.addAll(starShape( 0.5*sqrt(3), -0.5,     0  , -4*PI/3, 4, tipOffset));
				mainShape.addAll(Path.rotated(-PI/2, Path.translated(0, sqrt(3), parallelSegment)));
				mainShape.addAll(starShape( 0.5*sqrt(3),  0.5,    PI/3,   -PI/3, 2, tipOffset));
				// compose the shape of the Antarctica island
				List<Path.Command> sideShape = new ArrayList<>();
				sideShape.add(new Path.Command('M', 0, ySouthPole + shortEdge));
				sideShape.addAll(Path.translated(tipOffset/sqrt(2), -tipOffset/sqrt(2) + ySouthPole,
				                                 Path.rotated(-3*PI/4, Path.reversed(parallelSegment))));
				sideShape.add(new Path.Command('L', -shortEdge, ySouthPole));
				sideShape.addAll(Path.translated(tipOffset/sqrt(2), tipOffset/sqrt(2) + ySouthPole,
				                                 Path.rotated(-PI/4, Path.reversed(parallelSegment))));
				sideShape.add(new Path.Command('L', 0, ySouthPole - shortEdge));
				sideShape.addAll(Path.translated(-tipOffset/sqrt(2), tipOffset/sqrt(2) + ySouthPole,
				                                 Path.rotated(PI/4, Path.reversed(parallelSegment))));
				sideShape.add(new Path.Command('L', shortEdge, ySouthPole));
				sideShape.addAll(Path.translated(-tipOffset/sqrt(2), -tipOffset/sqrt(2) + ySouthPole,
				                                 Path.rotated(3*PI/4, Path.reversed(parallelSegment))));
				// determine the bounds and set the path 'M's and 'L's correctly
				double xMax = NEGATIVE_INFINITY, yMin = POSITIVE_INFINITY, yMax = NEGATIVE_INFINITY;
				for (int i = 0; i < mainShape.size(); i ++) {
					if (mainShape.get(i).args[0] > xMax)
						xMax = mainShape.get(i).args[0];
					if (mainShape.get(i).args[1] > yMax)
						yMax = mainShape.get(i).args[1];
					if (i > 0 && mainShape.get(i).type == 'M')
						mainShape.set(i, new Path.Command('L', mainShape.get(i).args));
				}
				for (int i = 0; i < sideShape.size(); i ++) {
					if (sideShape.get(i).args[1] < yMin)
						yMin = sideShape.get(i).args[1];
					if (i > 0 && sideShape.get(i).type == 'M')
						sideShape.set(i, new Path.Command('L', sideShape.get(i).args));
				}
				// puth the paths together and construct the shape directly
				mainShape.addAll(sideShape);
				return new Shape(-xMax, xMax, yMin, yMax, mainShape);
			}
		},
		
		/** The more compact zigzag configuration with Antarctica divided and attached */
		M_PROFILE(true, true) {
			Octant[] placeOctants(double tipOffset) {
				return new Octant[] {
					new Octant(-sqrt(3)/2, 0, -PI/6, -PI/2, PI/2, -3*PI/4),
					new Octant(-sqrt(3)/2, 0,  PI/6, -PI/2, PI/2,   -PI/4),
					new Octant( sqrt(3)/2, 0, -PI/6, -PI/2, PI/2,    PI/4),
					new Octant( sqrt(3)/2, 0,  PI/6, -PI/2, PI/2,  3*PI/4),
				};
			}
			Shape drawShape(double tipOffset, Projection projection) {
				List<Path.Command> shape = new ArrayList<>();
				shape.addAll(starShape( 0.0        , -0.5,  2*PI/3, -2*PI/3, 4, tipOffset));
				shape.addAll(starShape(-0.5*sqrt(3),  0.0,    PI/3,   -PI/3, 2, tipOffset));
				shape.addAll(starShape(-1.0*sqrt(3), -0.5,  2*PI/3,     0  , 2, tipOffset));
				shape.addAll(starShape(-1.0*sqrt(3), -1.5,    PI  ,  2*PI/3, 1, tipOffset));
				shape.addAll(starShape(-0.5*sqrt(3), -1.0,  5*PI/3,    PI/3, 4, tipOffset));
				shape.addAll(starShape( 0.0        , -1.5,  4*PI/3,  2*PI/3, 2, tipOffset));
				shape.addAll(starShape( 0.5*sqrt(3), -1.0,  5*PI/3,    PI/3, 4, tipOffset));
				shape.addAll(starShape( 1.0*sqrt(3), -1.5, -2*PI/3,   -PI  , 1, tipOffset));
				shape.addAll(starShape( 1.0*sqrt(3), -0.5,     0  , -2*PI/3, 2, tipOffset));
				shape.addAll(starShape( 0.5*sqrt(3),  0.0,    PI/3,   -PI/3, 2, tipOffset));
				return Shape.polygon(Path.asArray(shape));
			}
		},

		/** Gene Keyes's current configuration, with Antarctica reassembled in the center */
		M_PROFILE_WITH_SOUTH_POLE(true, false) {
			Octant[] placeOctants(double tipOffset) {
				double xSouthPole = -tipOffset/2.;
				double ySouthPole = -1.5 + tipOffset*sqrt(3)/2.;
				double quadrantLength = sqrt(3) - tipOffset;
				return new Octant[] {
					new Octant(-sqrt(3)/2, 0., -PI/6, toRadians(-65), PI/2, toRadians(-155)),
					new Octant(-sqrt(3)/2, 0.,  PI/6, toRadians(-90), PI/2, toRadians( -65)),
					new Octant( sqrt(3)/2, 0., -PI/6, toRadians(-65), PI/2, toRadians(  25)),
					new Octant( sqrt(3)/2, 0.,  PI/6, toRadians(-65), PI/2, toRadians( 115)),
					new Octant(xSouthPole - quadrantLength*sqrt(3)/2, ySouthPole - quadrantLength/2, 2*PI/3,
					           toRadians(-90), toRadians(-65), toRadians(-155)),
					new Octant(xSouthPole + quadrantLength*sqrt(3)/2, ySouthPole + quadrantLength/2, -PI/3,
					           toRadians(-90), toRadians(-65), toRadians(25)),
					new Octant(xSouthPole + quadrantLength/2, ySouthPole - quadrantLength*sqrt(3)/2, -5*PI/6,
					           toRadians(-90), toRadians(-65), toRadians(115)),
				};
			}
			Shape drawShape(double tipOffset, Projection projection) {
				double xSouthPole = -tipOffset/2.;
				double ySouthPole = -1.5 + tipOffset*sqrt(3)/2.;
				double shortEdge = tipOffset/(2*sind(15));
				List<Path.Command> parallelSegment = projection.drawLoxodrome(
						toRadians(65), toRadians(45), toRadians(65), toRadians(0), .1);
				parallelSegment.addAll(Path.scaled(-1, 1, Path.reversed(parallelSegment)));  // expand the parallel segment to cover 90°
				List<Path.Command> shape = new ArrayList<>();
				shape.addAll(starShape( 0.0        , -0.5,  2*PI/3, -2*PI/3, 4, tipOffset));
				shape.addAll(starShape(-0.5*sqrt(3),  0.0,    PI/3,   -PI/3, 2, tipOffset));
				shape.addAll(starShape(-1.0*sqrt(3), -0.5,  2*PI/3,     0  , 2, tipOffset));
				shape.addAll(Path.translated(-sqrt(3), -1.5, Path.rotated(5*PI/6, parallelSegment)));
				shape.addAll(starShape(-0.5*sqrt(3), -1.0,  5*PI/3,    PI/3, 4, tipOffset));
				shape.add(new Path.Command('L', xSouthPole - shortEdge*cosd(15), ySouthPole + shortEdge*sind(15)));
				shape.addAll(Path.translated(xSouthPole + tipOffset*sqrt(3)/2, ySouthPole + tipOffset/2,
				                             Path.rotated(-PI/3, Path.reversed(parallelSegment))));
				shape.add(new Path.Command('L', xSouthPole - shortEdge*sind(15), ySouthPole - shortEdge*cosd(15)));
				shape.addAll(Path.translated(xSouthPole - tipOffset/2, ySouthPole + tipOffset*sqrt(3)/2,
				                             Path.rotated(PI/6, Path.reversed(parallelSegment))));
				shape.add(new Path.Command('L', xSouthPole + shortEdge*cosd(15), ySouthPole - shortEdge*sind(15)));
				shape.addAll(Path.translated(xSouthPole - tipOffset*sqrt(3)/2, ySouthPole - tipOffset/2,
				                             Path.rotated(2*PI/3, Path.reversed(parallelSegment))));
				shape.add(new Path.Command('L', xSouthPole + shortEdge*sind(15), ySouthPole + shortEdge*cosd(15)));
				shape.addAll(Path.translated(0, -1.5, Path.rotated(5*PI/6, parallelSegment)));
				shape.addAll(starShape( 0.5*sqrt(3), -1.0,  5*PI/3,    PI/3, 4, tipOffset));
				shape.addAll(Path.translated(sqrt(3), -1.5, Path.rotated(-5*PI/6, parallelSegment)));
				shape.addAll(starShape( 1.0*sqrt(3), -0.5,     0  , -2*PI/3, 2, tipOffset));
				shape.addAll(starShape( 0.5*sqrt(3),  0.0,    PI/3,   -PI/3, 2, tipOffset));
				return Shape.polygon(Path.asArray(shape));
			}
		},

		/** Luca Concialdi's "Bat" arrangement */
		BAT_SHAPE(true, false) {
			Octant[] placeOctants(double tipOffset) {
				return rotateOctants(toRadians(5), new Octant[]{
					new Octant( 0.0,  0.0        , -2*PI/3,   0  ,  PI/2, toRadians(-160), toRadians(-9), 1), // Alaska
					new Octant(-1.5,  0.5*sqrt(3),     0  , -PI/2,   0  , toRadians(-160), toRadians(9), 1),  // West Antarctica
					new Octant( 0.0,  0.0        ,   -PI/3, -PI/2,  PI/2, toRadians( -70)),                   // America
					new Octant( 0.0,  0.0        ,     0  , -PI/4,  PI/2, toRadians(  20)),                   // Africa
					new Octant( 0.0, -1.0*sqrt(3), -2*PI/3, -PI/2, -PI/4, toRadians(  20), -1, 0),            // Coats Land
					new Octant( 0.0, -1.0*sqrt(3),  2*PI/3, -PI/2, -PI/4, toRadians(  20),  0, 1),            // Dronning Maund Land
					new Octant( 0.0,  0.0        ,    PI/3, -PI/2,  PI/2, toRadians( 110)),                   // Asia
					new Octant( 0.0,  0.0        ,  2*PI/3,   0  ,  PI/2, toRadians( 200), -1, toRadians(-9)),// Chukotka
					new Octant( 1.5,  0.5*sqrt(3),     0  , -PI/2,   0  , toRadians( 200), -1, toRadians(9)), // New Zealand
				});
			}
			Shape drawShape(double tipOffset, Projection projection) {
				List<Path.Command> meridian = projection.drawLoxodrome(0, toRadians(-9), PI/2, toRadians(-9), .1);
				List<Path.Command> parallel = projection.drawLoxodrome(PI/4, toRadians(0), PI/4, toRadians(45), .1);
				List<Path.Command> shape = new ArrayList<>();
				shape.addAll(starShape( 0.0,  0.0        ,    PI/2,   -PI/2, 3, tipOffset));
				shape.addAll(Path.rotated(-2*PI/3, Path.reversed(meridian)));
				shape.addAll(starShape(-1.0,  0.0        ,  5*PI/6,   -PI/2, 4, tipOffset));
				shape.addAll(Path.translated(-1.5, -0.5*sqrt(3), Path.rotated(PI, meridian)));
				shape.addAll(starShape(-1.5, -0.5*sqrt(3),  5*PI/6,    PI/2, 1, tipOffset));
				shape.addAll(Path.translated(-1.5, -0.5*sqrt(3), Path.rotated(PI/3, parallel)));
				shape.addAll(starShape(-0.5, -0.5*sqrt(3),  3*PI/2,    PI/6, 4, tipOffset));
				shape.addAll(Path.translated(0.0, -1.0*sqrt(3), Path.rotated(PI, Path.reversed(parallel))));
				shape.addAll(Path.translated(0.0, -1.0*sqrt(3), Path.rotated(PI, Path.scaled(-1, 1, parallel))));
				shape.addAll(starShape( 0.5, -0.5*sqrt(3),   -PI/6, -3*PI/2, 4, tipOffset));
				shape.addAll(Path.translated(1.5, -0.5*sqrt(3), Path.rotated(-PI/3, Path.scaled(-1, 1, Path.reversed(parallel)))));
				shape.addAll(starShape( 1.5, -0.5*sqrt(3),   -PI/2, -5*PI/6, 1, tipOffset));
				shape.addAll(Path.translated(1.5, -0.5*sqrt(3), Path.rotated(PI, Path.reversed(meridian))));
				shape.addAll(starShape( 1.0,  0.0        ,    PI/2, -5*PI/6, 4, tipOffset));
				shape.addAll(Path.rotated(2*PI/3, meridian));
				return Shape.polygon(Path.asArray(Path.rotated(toRadians(5), shape)));
			}
		},

		/** an octahedron that actually only covers the positive octant, in case you want to do each octant separately */
        SINGLE_OCTANT(false, true) {
			Octant[] placeOctants(double tipOffset) {
				return new Octant[]{
						new Octant(0.0, 0.0, PI/6, 0, PI/2, PI/4),
						};
			}
			Shape drawShape(double tipOffset, Projection projection) {
				double cutDepth = tipOffset*(cotd(15) + sqrt(3))/2;
				return Shape.polygon(new double[][] {
						{cutDepth*sqrt(3)/2., -cutDepth/2.},
						{tipOffset/2., -tipOffset*sqrt(3)/2.},
						{0, -cutDepth},
						{0, -1 + cutDepth},
						{tipOffset/2., -1 + tipOffset*sqrt(3)/2.},
						{cutDepth*sqrt(3)/2., -1 + cutDepth/2.},
						{sqrt(3)/2. - cutDepth*sqrt(3)/2, -1/2. - cutDepth/2.},
						{sqrt(3)/2. - tipOffset, -1/2.},
						{sqrt(3)/2. - cutDepth*sqrt(3)/2., -1/2. + cutDepth/2.},
				});
			}
		};
		
		public final boolean finite;
		public final boolean hasAspect;

		/**
		 * @param finite whether there are enuff octants to project every part of the earth
		 * @param hasAspect whether it would make any sense to change the aspect
		 */
		Configuration(boolean finite, boolean hasAspect) {
			this.finite = finite;
			this.hasAspect = hasAspect;
		}

		/**
		 * construct the array of octant specifications.  each octant must specify the x and y of the pole,
		 * the central meridian, and bounding latitudes and longitudes of the face.
		 * @param tipOffset if the projection is for a truncated octahedron, this is the distance from where the
		 *                  tip of the octant would be for a true octahedron to where it actually is
		 */
		abstract Octant[] placeOctants(double tipOffset);

		/**
		 * construct the bounding polygon of a projection using this configuration.
		 * @param tipOffset if the projection is for a truncated octahedron, this is the distance from where the
		 *                  tip of the octant would be for a true octahedron to where it actually is
		 * @param projection the projection to be used for each face, in case meridians and parallels are involved
		 *                   in the shape TODO: right now I'm using the full projection but ideally this should be a projection for just one octant
		 */
		abstract Shape drawShape(double tipOffset, Projection projection);

		/**
		 * generate a path representing a partial star shape, comprising some number of points, separated by
		 * notches, through some set of angles.  the notches will always be right angles.
		 * @param x0 the x coordinate of the center of the star
		 * @param y0 the y coordinate of the center of the star
		 * @param th0 the bearing from the center to the first point (0 is down, increasing goes rite)
		 * @param th1 the bearing from the center to the last point (0 is down, increasing goes rite)
		 * @param numPoints the number of points (the first and last ones are tecnicly only half points so depending on how you count it may also be the number of points minus one)
		 * @param innerRadius the distance from the center to the nearest vertex
		 * @return an open path tracing from the first point to the last point
		 */
		private static List<Path.Command> starShape(double x0, double y0, double th0, double th1, int numPoints, double innerRadius) {
			if (innerRadius > 0) {
				double outerRadius = innerRadius*(cotd(15) + sqrt(3))/2;
				List<Path.Command> shape = new ArrayList<>(2*numPoints + 1);
				shape.add(new Path.Command('M', x0 + outerRadius*sin(th0), y0 - outerRadius*cos(th0)));
				for (int i = 0; i < numPoints; i ++) {
					double th = th0 + (i + 0.5)*(th1 - th0)/numPoints;
					shape.add(new Path.Command('L', x0 + innerRadius*sin(th), y0 - innerRadius*cos(th)));
					th = th0 + (i + 1)*(th1 - th0)/numPoints;
					shape.add(new Path.Command('L', x0 + outerRadius*sin(th), y0 - outerRadius*cos(th)));
				}
				return shape;
			}
			else {
				return List.of(new Path.Command('M', x0, y0));
			}
		}

		/**
		 * take some Octants defining a map layout and figure out what octants would correspond to that may
		 * rotated counterclockwise by some angle.
		 * @param th the angle in radians, with positive meaning widdershins and negative meaning clockwise
		 * @param in the unrotated octants
		 */
		private static Octant[] rotateOctants(double th, Octant[] in) { // apply that rotation to the bat
			Octant[] out = new Octant[in.length];
			for (int i = 0; i < in.length; i ++) {
				out[i] = new Octant(
						in[i].x*cos(th) - in[i].y*sin(th),
						in[i].x*sin(th) + in[i].y*cos(th),
						in[i].planeRotation + th,
						in[i].minLatitude, in[i].maxLatitude,
						in[i].centralLongitude,
						in[i].minLongitude, in[i].maxLongitude);
			}
			return out;
		}
	}


	/**
	 * the information that defines a contiguusly projecting region of the globe and where/how it's placed/oriented in
	 * the plane.  in some situations this may be used to represent a full quadrant of the globe (bounded in longitude
	 * but not in latitude).
	 */
	private static class Octant {
		public final double x;
		public final double y;
		public final double planeRotation;
		public final double minLatitude;
		public final double maxLatitude;
		public final double centralLongitude;
		public final double minLongitude;
		public final double maxLongitude;

		Octant(double x, double y, double planeRotation, double minLatitude, double maxLatitude, double centralLongitude) {
			this(x, y, planeRotation, minLatitude, maxLatitude, centralLongitude, -PI, PI);
		}

		Octant(double x, double y, double planeRotation, double minLatitude, double maxLatitude,
		       double centralLongitude, double minLongitude, double maxLongitude) {
			this.x = x;
			this.y = y;
			this.planeRotation = planeRotation;
			this.minLatitude = minLatitude;
			this.maxLatitude = maxLatitude;
			this.centralLongitude = centralLongitude;
			this.minLongitude = minLongitude;
			this.maxLongitude = maxLongitude;
		}
	}
	
}
