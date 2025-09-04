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
import maps.Projection.Type;
import utils.NumericalAnalysis;
import utils.Shape;

import java.util.Arrays;

import static java.lang.Double.NaN;
import static java.lang.Double.isNaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.hypot;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.lang.Math.sinh;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import static java.lang.Math.toRadians;
import static utils.Math2.atan;
import static utils.NumericalAnalysis.polynomialInterpolate;

/**
 * A class containing map projections with four-way symmetry and curved parallels.
 * 
 * @author jkunimune
 */
public class Lenticular {
	
	public static final Projection AITOFF = new Projection(
			"Aitoff", "David A. Aitoff", "A compromise projection shaped like an ellipse",
			Shape.ellipse(Math.PI, Math.PI/2), true, true, true, true,
			Type.PSEUDOAZIMUTHAL, Property.COMPROMISE, 2) {

		public double[] project(double lat, double lon) {
			final double a = acos(cos(lat)*cos(lon/2));
			if (a == 0) 	return new double[] {0, 0};
			return new double[] {
					2*cos(lat)*sin(lon/2)*a/sin(a),
					sin(lat)*a/sin(a)};
		}
		
		public double[] inverse(double x, double y) {
			final double[] intermediate = Azimuthal.EQUIDISTANT.inverse(x/2, y);
			double[] transverse = transformToOblique(intermediate, new double[] {0, 0, 0});
			if (transverse != null) 	transverse[1] *= 2;
			return transverse;
		}
	};
	
	
	public static final Projection HAMMER = new Projection(
			"Hammer", "Ernst H. H. von Hammer", "An equal-area projection shaped like an ellipse",
			Shape.ellipse(2, 1), true, true, true, true, Type.PSEUDOAZIMUTHAL, Property.EQUAL_AREA, 1) {

		public double[] project(double lat, double lon) {
			final double z = sqrt(1+cos(lat)*cos(lon/2));
			return new double[] {2*cos(lat)*sin(lon/2)/z, sin(lat)/z};
		}
		
		public double[] inverse(double x, double y) {
			final double z = sqrt(1 - x*x/8 - y*y/2);
			final double shift = (hypot(x/2, y) > 1) ? 2*PI*signum(x) : 0;
			return new double[] {
					asin(z*y*sqrt(2)),
					2*atan(sqrt(.5)*z*x / (2*z*z - 1)) + shift};
		}
	};
	
	
	public static final Projection VAN_DER_GRINTEN = new Projection(
			"Van der Grinten", "Alphons J. van der Grinten", "A circular compromise map that is popular for some reason",
			Shape.circle(1), true, true, true, true, Type.OTHER, Property.COMPROMISE, 0) {

		public double[] project(double lat, double lon) {
			if (lat == 0) //special case 1: equator
				return new double[] {lon/PI, 0};
			if (lon == 0 || lat >= PI/2 || lat <= -PI/2) //special case 3: prime meridian
				return new double[] {0, tan(asin(2*lat/PI)/2)};
			
			final double t = abs(asin(2*lat/PI));
			final double A = abs(PI/lon - lon/PI)/2;
			final double G = cos(t)/(sin(t)+cos(t)-1);
			final double P = G*(2/sin(t) - 1);
			final double Q = A*A + G;
			return new double[] {
					signum(lon)*(A*(G-P*P)+sqrt(A*A*(G-P*P)*(G-P*P)-(P*P+A*A)*(G*G-P*P)))/(P*P+A*A),
					signum(lat)*(P*Q-A*sqrt((A*A+1)*(P*P+A*A)-Q*Q))/(P*P+A*A)};
		}
		
		public double[] inverse(double x, double y) {
			if (y == 0) // special case 1: equator
				return new double[] {0, x*PI};
			if (x == 0) // special case 3: prime meridian
				return new double[] {PI/2 * sin(2*atan(y)), 0};

			double r4 = pow(x*x + y*y, 2);
			double c1 = -abs(y) * (1 + x*x + y*y);
			double c2 = c1 - 2*y*y + x*x;
			double c3 = -2 * c1 + 1 + 2*y*y + r4;
			double d = y*y / c3 + 1 / 27.0 * (2*pow(c2 / c3, 3) - 9*c1*c2 / (c3*c3));
			double a1 = 1 / c3*(c1 - c2*c2 / (3*c3));
			double m1 = 2 * sqrt(-a1 / 3);
			double t1 = acos(3*d / (a1 * m1)) / 3;
			return new double[] {
					signum(y) * PI * (-m1 * cos(t1 + PI/3) - c2 / (3*c3)),
					PI*(x*x + y*y - 1 + sqrt(1 + 2*(x*x - y*y) + r4))
							/ (2*x)};
		}
	};
	
	
	public static final Projection STREBE_95 = new Projection(
			"Strebe 1995", "Daniel Strebe", "An equal-area map with curvy poles that pushes distortion to the edges",
			null, true, true, false, true, Type.STREBE, Property.COMPROMISE, 2,
			new String[] {"Scale Factor"},
			new double[][] {{sqrt(2*PI/(4+PI)), sqrt((4+PI)/PI*2), 1.35}}) {
		
		private double factor;
		
		public void initialize(double... params) {
			factor = params[0];
			shape = Shape.meridianEnvelope(this);
		}
		
		public double[] project(double lat, double lon) {
			double[] xy1 = Pseudocylindrical.ECKERT_IV.project(lat, lon);
			xy1[0] *= 2*sqrt(PI/(4+PI))*factor/sqrt(2);
			xy1[1] *= 2*sqrt(PI/(4+PI))/factor/sqrt(2);
			double[] ll2 = Pseudocylindrical.MOLLWEIDE.inverse(xy1);
			double[] xy3 = Lenticular.HAMMER.project(ll2);
			xy3[0] *= 1/factor;
			xy3[1] *= factor;
			return xy3;
		}
		
		public double[] inverse(double x, double y) {
			double[] ll2 = Lenticular.HAMMER.inverse(x*factor, y/factor);
			double[] xy1 = Pseudocylindrical.MOLLWEIDE.project(ll2);
			xy1[0] /= 2*sqrt(PI/(4+PI))*factor/sqrt(2);
			xy1[1] /= 2*sqrt(PI/(4+PI))/factor/sqrt(2);
			double[] ll0 = Pseudocylindrical.ECKERT_IV.inverse(xy1);
			
			if (isNaN(ll0[0]))
				return null;
			else
				return ll0;
		}
		
	};
	
	
	public static final Projection BERTIN = new Projection(
			"Bertin", "Jacques Bertin", "An artistically conceived oblique map projection",
			Shape.ellipse(1.68, 1), true, true, true, false, Type.OTHER, Property.COMPROMISE, 3) {
		
		private final double[] POLE = {toRadians(42), toRadians(-163.5), toRadians(180)};

		public double[] project(double lat, double lon) {
			double[] oblique = transformFromOblique(lat, lon, POLE); // start with a slightly oblique globe
			lat = oblique[0];
			lon = oblique[1];
			if (lat + lon < -1.4) { // apply controlled smooshing to the resulting coordinates
				double u = (lon - lat + 1.6)*(lat + lon + 1.4) / 8;
				lon += u;
				lat -= 0.8 * u * cos(lat);
			}
			
			double[] coords = HAMMER.project(oblique); // apply a Hammer projection
			double x = coords[0], y = coords[1];
			x *= 1.68/2; // change the aspect ratio
			double d = (1 - cos(lat * lon)) / 12; // apply controlled smooshing to the resulting coordinates
			if (y < 0)
				x *= 1 + d; // depending on whether it is in the top
			else
				y *= 1 + d / 1.5 * x*x; // or bottom half
			return new double[] {x, y};
		}
		
		
		public double[] inverse(double x, double y) {
//			return obliquifyPlnr(HAMMER.inverse(x/1.63*2, y), POLE); // TODO
			return null;
		}
	};
	
	
	public static final Projection LAGRANGE = new Projection(
			"Lagrange", "Johann H. Lambert", "A circular conformal map",
			null, true, true, true, true, Type.OTHER, Property.CONFORMAL, 2,
			new String[] {"Central parallel", "Longitudinal scale"},
			new double[][] {{-89., 89., 0.}, {0.01, 1.5, 0.5}}) {
		
		private double prefactor, n, lonMax;
		
		public double[] project(double lat, double lon) {
			if (abs(lat) == PI/2)
				return new double[] { 0, signum(lat) };
			lon = max(-lonMax, min(lonMax, lon)); // don't let longitude exceed ± lonMax lest x -> inf
			double v = pow(prefactor*(1 + sin(lat))/(1 - sin(lat)), n/2);
			double c = (v + 1/v)/2 + cos(n*lon);
			double x = sin(n*lon)/c;
			double y = (v - 1/v)/(2*c);
			return new double[] {x, y};
		}
		
		public double[] inverse(double x, double y) {
			double r2 = x*x + y*y;
			double th = 2*y/(1 + r2);
			double t = pow((1 + th)/(1 - th), 1/n);
			double lat = asin((t - prefactor)/(t + prefactor));
			double lon = 1/n*(atan(2*x/(1 - r2)) + ((r2 > 1) ? PI*signum(x) : 0));
			return new double[] {lat, lon};
		}
		
		public void initialize(double... params) {
			double lat0 = Math.toRadians(params[0]);
			this.prefactor = (1 - sin(lat0))/(1 + sin(lat0));
			this.n = params[1];
			// put a limit on longitudes so you don't get self-intersection when n > 1
			this.lonMax = 0.999*PI/n;
			// we need to take care with the shape to avoid it getting too big
			double maxSize = PI*n;
			double lonEnvelope = min(PI, lonMax);
			double xIntercept = sin(n*lonEnvelope)/(1 + cos(n*lonEnvelope));
			double envelopeRadius = (pow(xIntercept, 2) + 1)/(2*xIntercept);
			double envelopeHeight = (n < 1/2.) ? 1 : envelopeRadius;
			double xCenter = xIntercept - envelopeRadius;
			double xOffset = sqrt(pow(envelopeRadius, 2) - pow(maxSize, 2));
			double yIntersect = sqrt(pow(envelopeRadius, 2) - pow(maxSize - xCenter, 2));
			// if the poles are too pointy
			if (n < 1/2. && envelopeHeight > maxSize) {
				this.shape = new Shape( // cut them off
						-xIntercept, xIntercept, -maxSize, maxSize,
						Arrays.asList(
								new Path.Command('M', xCenter + xOffset, -maxSize),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, xCenter + xOffset, maxSize),
								new Path.Command('H', -xCenter - xOffset),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, -xCenter - xOffset, -maxSize),
								new Path.Command('Z')));
			}
			// if the whole thing fits within the square
			else if (xIntercept <= maxSize) {
				// just use the meridian envelope
				this.shape = new Shape(
						-xIntercept, xIntercept, -envelopeHeight, envelopeHeight,
						Arrays.asList(
								new Path.Command('M', 0, -1),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, (n > 1/2.) ? 1 : 0, 1, 0, 1),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, (n > 1/2.) ? 1 : 0, 1, 0, -1),
								new Path.Command('Z')));
			}
			// if it breaches the side of the square but not the top
			else if (envelopeRadius <= maxSize) {
				// cut the sides off
				this.shape = new Shape(
						-maxSize, maxSize, -envelopeHeight, envelopeHeight,
						Arrays.asList(
								new Path.Command('M', maxSize, yIntersect),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, 0, 1),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, -maxSize, yIntersect),
								new Path.Command('V', -yIntersect),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, 0, -1),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, maxSize, -yIntersect),
								new Path.Command('Z')));
			}
			// NOTE: there's another case here where it breaches the top but doesn't encompass the corners.
			// however, this happens for such a small part of parameter space that it's currently not possible to save
			// a map in this case, and even if it was the corner artifact would be almost imperceptible.
			// so I'm leaving that case out.
			// if it breaches the top as well and encompasses the corners
			else {
				// it becomes a rectangle with just a few cuts out of it
				this.shape = new Shape(
						-maxSize, maxSize, -maxSize, maxSize,
						Arrays.asList(
								new Path.Command('M', maxSize, maxSize),
								new Path.Command('H', xCenter - xOffset),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, 0, 1),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, -xCenter + xOffset, maxSize),
								new Path.Command('H', -maxSize),
								new Path.Command('V', -maxSize),
								new Path.Command('H', -xCenter + xOffset),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, 0, -1),
								new Path.Command('A', envelopeRadius, envelopeRadius, 0, 0, 1, xCenter - xOffset, -maxSize),
								new Path.Command('H', maxSize),
								new Path.Command('Z')));
			}
		}
	};
	
	
	public static final Projection EISENLOHR = new Projection(
			"Eisenlohr", "Friedrich Eisenlohr", "The optimal conventional conformal map",
			null, true, true, true, false, Type.OTHER, Property.CONFORMAL, 2) {
		
		public double[] project(double lat, double lon) {
			if (abs(lat) == PI/2)
				return new double[] { 0, signum(lat)*(1 - PI/4) };
			Complex w = new Complex(lon, log(tan(PI/4+lat/2)));
			Complex v = w.divide(4).minus(PI/8).tan().minus(1).divide(-sqrt(2));
			Complex z = v.log().plus(v.invert().minus(v).divide(sqrt(2)));
			return new double[] { z.getRe(), z.getIm() };
		}
		
		public double[] inverse(double x, double y) {
			if (x > 0) { // It converges on the right half, but not the left
				double[] res = inverse(-x, -y); // I'm not sure why it does this,
				return new double[] {-res[0], -res[1]}; // but the fix is easy.
			}
			
			Complex z = new Complex(x, y);
			Complex v = NumericalAnalysis.newtonRaphsonApproximation(z, z.exp(),
					(t)->(t.log().plus(t.invert().minus(t).divide(sqrt(2)))),
					(t)->(t.invert().plus(t.pow(-2).neg().minus(1).divide(sqrt(2)))),
					1e-4);
			Complex w = atan(v.times(sqrt(2)).minus(1)).minus(PI/8).times(-4);
			return new double[] { atan(sinh(w.getIm())), w.getRe() };
		}

		public void initialize(double... params) throws IllegalArgumentException {
			this.shape = Shape.meridianEnvelope(this);
		}
	};
	
	
	public static final Projection WAGNER_VIII = new Projection(
			"Wagner VIII", "Karlheinz Wagner", "A compromise projection with pseudoazimuthal energy",
			null, true, true, true, true, Type.OTHER, Property.COMPROMISE, 3) {

		private final double m1 = 0.92118, m2 = 0.8855, n = 3.,
				cX = 5.62290, cY = 2.61626;
		
		public double[] project(double lat, double lon) {
			double psi = asin(m1*sin(m2*lat));
			if (lon == 0)
				return new double[] { 0, cY*sin(psi/2) };
			double del = acos(cos(lon/n)*cos(psi));
			double alp = signum(lon)*acos(sin(psi)/sin(del));
			return new double[] {
					cX*sin(del/2)*sin(alp),
					cY*sin(del/2)*cos(alp) };
		}
		
		public double[] inverse(double x, double y) {
			double alp = atan2(x/cX, y/cY);
			double del = 2*asin(hypot(x/cX, y/cY));
			double psi = asin(cos(alp)*sin(del));
			double lat = asin(sin(psi)/m1)/m2;
			double lon = signum(x)*acos(cos(del)/cos(psi))*n;
			if (isNaN(lat) || abs(lat) > PI/2)
				return null;
			else
				return new double[] {lat, lon};
		}

		public void initialize(double... params) throws IllegalArgumentException {
			this.shape = Shape.meridianEnvelope(this);
		}
	};
	
	
	public static final Projection POLYCONIC = new Projection(
			"American polyconic", "Ferdinand R. Hassler", "A map made for narrow strips of longitude that was really popular with the USGS for a while",
			null, true, true, true, false, Type.OTHER, Property.EQUIDISTANT, 3) {
		
		public double[] project(double lat, double lon) {
			if (lat == 0)
				return new double[] {lon, 0};
			double E = lon*sin(lat);
			return new double[] { sin(E)/tan(lat), lat + (1 - cos(E))/tan(lat) };
		}
		
		public double[] inverse(double x, double y) {
			if (y == 0)
				return new double[] {0, x};
			
			double lat;
			try {
				lat = NumericalAnalysis.bisectionFind(
						(ph) -> {
							if (ph == 0)
								return y;
							else
								return 1/tan(ph) - signum(ph)*hypot(x, y - ph - 1/tan(ph));
						},
						-PI/2, PI/2, 1e-4);
			} catch (NumericalAnalysis.AlgorithmFailedException e) {
				lat = NaN;
			}
			if (isNaN(lat))
				return null;
			return new double[] { lat, atan2(x, -signum(lat)*(y - lat - 1/tan(lat)))/abs(sin(lat)) };
		}

		public void initialize(double... params) throws IllegalArgumentException {
			this.shape = Shape.meridianEnvelope(this);
		}
	};
	
	
	public static final Projection CHINA = new Projection(
			"Equal difference latitude polyconic", "Chinese Bureau of Surveying and Mapping",
			"Also known as the 等差分纬线多圆锥 projection – a compromise projection that's extremely common in China.  This implementation uses the equations published by Yè Yuǎnzhì in Cèhuì Tōngbào (2012 Nov), p. 68.",
			null, true, true, true, false, Type.OTHER, Property.COMPROMISE, 3) {
		
		private final double SCALE = 6371116;
		private final double[] TABLE_LATITUDE = {10, 15, 20, 23.44, 30, 40, 45, 50, 60, 66.56, 70, 75, 80, 90};
		private final double[] TABLE_RADIUS = {2399.943, 1643.348, 1201.174, 1024.159, 807.634, 618.060, 556.629, 507.561, 430.763, 388.633, 375.478, 372.030, 388.525, 508.062};
		private final double[] TABLE_ARC = {3.9029, 5.6360, 7.5961, 8.7978, 10.8572, 13.4857, 14.5329, 15.4043, 16.6511, 17.2088, 17.0908, 16.1835, 14.4826, 9.3451};
		private double[] referenceLatitude;
		private double[] referenceCurvature;
		private double[] referenceParallelLength;
		
		public double[] project(double lat, double lon) {
			double y0 = (0.9953537*lat + 0.01476138*lat*lat*lat);
			double ρ = 1/polynomialInterpolate(lat, referenceLatitude, referenceCurvature, 3);
			if (Double.isFinite(ρ)) {
				double δn = polynomialInterpolate(lat, referenceLatitude, referenceParallelLength, 3)/ρ;
				double δ = δn*1.1*(1 - 1/11.*abs(lon/PI))*lon/PI;
				return new double[] {ρ*sin(δ), y0 + ρ*(1 - cos(δ))};
			}
			else {
				double xn = polynomialInterpolate(lat, referenceLatitude, referenceParallelLength, 3);
				double x = xn*1.1*(1 - 1/11.*abs(lon/PI))*lon/PI;
				return new double[] {x, y0};
			}
		}
		
		public double[] inverse(double x, double y) {
			double lat;
			try {
				lat = NumericalAnalysis.bisectionFind(
						(ph) -> {
							if (ph == 0)
								return y;
							double y0 = (0.9953537*ph + 0.01476138*ph*ph*ph);
							double ρ = 1/polynomialInterpolate(ph, referenceLatitude, referenceCurvature, 3);
							return ρ - signum(ρ)*hypot(x, y - y0 - ρ);
						},
						-PI/2, PI/2, shape.height*1e-4);
			} catch (NumericalAnalysis.AlgorithmFailedException e) {
				lat = NaN;
			}
			if (isNaN(lat))
				return null;
			else if (lat == 0) {
				double xn = polynomialInterpolate(0, referenceLatitude, referenceParallelLength, 3);
				double lon = signum(x)*PI*(5.5 - sqrt(30.25 - 10*abs(x/xn)));
				return new double[] {0, lon};
			}

			double y0 = (0.9953537*lat + 0.01476138*lat*lat*lat);
			double ρ = 1/polynomialInterpolate(lat, referenceLatitude, referenceCurvature, 3);
			double δn = polynomialInterpolate(lat, referenceLatitude, referenceParallelLength, 3)/ρ;
			double δ = atan2(x, -signum(ρ)*(y - y0 - ρ));
			double lon = signum(δ)*PI*(5.5 - sqrt(30.25 - 10*abs(δ/δn)));
			return new double[] { lat, lon };
		}
		
		public void initialize(double... params) throws IllegalArgumentException {
			referenceLatitude = new double[2*TABLE_LATITUDE.length + 1];
			referenceCurvature = new double[2*TABLE_LATITUDE.length + 1];
			referenceParallelLength = new double[2*TABLE_LATITUDE.length + 1];
			for (int i = 0; i < referenceLatitude.length; i ++) {
				if (i != TABLE_LATITUDE.length) {
					int j, sign;
					if (i < TABLE_LATITUDE.length) {
						j = TABLE_LATITUDE.length - 1 - i;
						sign = -1;
					}
					else {
						j = i - TABLE_LATITUDE.length - 1;
						sign = 1;
					}
					referenceLatitude[i] = sign*toRadians(TABLE_LATITUDE[j]);
					referenceCurvature[i] = sign/(TABLE_RADIUS[j]*1e5/SCALE);
					referenceParallelLength[i] = TABLE_RADIUS[j]*1e5/SCALE*toRadians(TABLE_ARC[j]);
				}
				else {
					referenceLatitude[i] = 0;
					referenceCurvature[i] = 0;
					referenceParallelLength[i] = 16500000/SCALE;
				}
			}
			this.shape = Shape.meridianEnvelope(this);
		}
	};

}
