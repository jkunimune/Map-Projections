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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.function.DoubleConsumer;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import dialogs.ProgressBarDialog;
import ellipticFunctions.Jacobi;
import util.Dixon;
import util.Elliptic;
import util.Math2;
import util.NumericalAnalysis;
import util.Vector;

/**
 * Map projections!
 * 
 * @author jkunimune
 */
public enum Projection {

	NULL("null", "", 1., 0b0000, "", "", new String[0], new double[0][]) {
		public double[] project(double lat, double lon, double[] params) {
			return null;
		}
		public double[] inverse(double x, double y, double[] params) {
			return null;
		}
	},
	
	MERCATOR("Mercator", 1., 0b0111, "cylindrical", "conformal", "very popular") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {lon, Math.log(Math.tan(Math.PI/4+lat/2))};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] {Math.atan(Math.sinh(y*Math.PI)), x*Math.PI};
		}
	},
	
	PLATE_CARREE("Plate Carrée", 2., 0b1111, "cylindrical", "equidistant", null, "focused on the equator") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {lon, lat};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] {y*Math.PI/2, x*Math.PI};
		}
	},
	
	EQUIRECTANGULAR("Equirectangular", "A linear mapping from longitude and latitude to x and y",
			2, 0b1111, "cylindrical", "equidistant",
			new String[]{"Std. parallel"}, new double[][]{{0, 85, 0}}) {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {lon, lat/Math.cos(Math.toRadians(params[0]))};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] {y*Math.PI/2, x*Math.PI};
		}
		public double getAspectRatio(double[] params) {
			return Math.PI*Math.pow(Math.cos(Math.toRadians(params[0])), 2);
		}
	},
	
	GALL_PETERS("Gall-Peters", 1.571, 0b1111, "cylindrical", "equal-area", "somewhat controversial") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {lon, Math.sin(lat)*Math.PI/1.571};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
	},
	
	HOBO_DYER("Hobo-Dyer", 1.977, 0b1111, "cylindrical", "equal-area",
			null, "with least distortion at 37.5°") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {lon, Math.sin(lat)*Math.PI/1.977};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
	},
	
	BEHRMANN("Behrmann", 2.356, 0b1111, "cylindrical", "equal-area",
			null, "with least distortion at 30°") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {lon, Math.sin(lat)*Math.PI/2.356};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
	},
	
	LAMBERT_CYLIND("Lambert cylindrical", Math.PI, 0b1111, "cylindrical", "equal-area",
			null, "with least distortion along the equator") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {lon, Math.sin(lat)};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
	},
	
	E_A_CYLIND("Equal-area cylindrical", "A generalized equal-area cylindrical projection",
			1.977, 0b1111, "cylindrical", "equal-area",
			new String[]{"Std. parallel"}, new double[][]{{0, 85, 37.5}}) {
		public double[] project(double lat, double lon, double[] params) {
			final double a = Math.pow(Math.cos(Math.toRadians(params[0])), 2);
			return new double[] {lon, Math.sin(lat)/a};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { Math.asin(y), x*Math.PI };
		}
		public double getAspectRatio(double[] params) {
			return Math.PI*Math.pow(Math.cos(Math.toRadians(params[0])), 2);
		}
	},
	
	GALL("Gall Stereographic", 4/3., 0b1111, "cylindrical", "compromise") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {lon, Math.tan(lat/2)*(1+Math.sqrt(2))};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { 2*Math.atan(y), x*Math.PI };
		}
	},
	
	MILLER("Miller", 4*Math.PI/5/Math.log(Math.tan(9*Math.PI/20)), 0b1111, "cylindrical", "compromise") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {lon, 1.25*Math.log(Math.tan(Math.PI/4+.8*lat/2))};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] {1.25*Math.atan(Math.sinh(y*Math.log(Math.tan(9*Math.PI/20)))), x*Math.PI};
		}
	},
	
	STEREOGRAPHIC("Stereographic", 1., 0b0111, "azimuthal", "conformal", "mathematically important") {
		public double[] project(double lat, double lon, double[] params) {
			final double r = 1.5/(Math.tan(lat/2 + Math.PI/4));
			return new double[] {r*Math.sin(lon), -r*Math.cos(lon)};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { Math.PI/2 - 2*Math.atan(2*Math.hypot(x, y)),
					Math.atan2(y, x) + Math.PI/2 };
		}
	},
	
	POLAR("Polar", 1., 0b0111, "azimuthal", "equidistant") {
		public double[] project(double lat, double lon, double[] params) {
			final double r = Math.PI/2 - lat;
			return new double[] {r*Math.sin(lon), -r*Math.cos(lon)};
		}
		public double[] inverse(double x, double y, double[] params) {
			double phi = Math.PI/2 - Math.PI * Math.hypot(x, y);
			if (phi > -Math.PI/2)
				return new double[] {phi, Math.atan2(y, x) + Math.PI/2};
			else
				return null;
		}
	},
	
	E_A_AZIMUTH("Azimuthalal Equal-Area", 1., 0b0111, "azimuthal", "equal-area") {
		public double[] project(double lat, double lon, double[] params) {
			final double r = Math.PI*Math.cos((Math.PI/2+lat)/2);
			return new double[] {r*Math.sin(lon), -r*Math.cos(lon)};
		}
		public double[] inverse(double x, double y, double[] params) {
			double R = Math.hypot(x, y);
			if (R <= 1)
				return new double[] {Math.asin(1-2*R*R), Math.atan2(y,x)+Math.PI/2};
			else
				return null;
		}
	},
	
	ORTHOGRAPHIC("Orthographic", "A projection that mimics the Earth viewed from a great distance",
			1., 0b1110, "azimuthal", "orthographic") {
		public double[] project(double lat, double lon, double[] params) {
			if (lat < 0)	lat = 0;
			final double r = Math.PI*Math.cos(lat);
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon) };
		}
		public double[] inverse(double x, double y, double[] params) {
			double R = Math.hypot(x, y);
			if (R <= 1)
				return new double[] { Math.acos(R), Math.atan2(y, x) + Math.PI/2 };
			else
				return null;
		}
	},
	
	GNOMONIC("Gnomonic", "A projection that draws all great circles as straight lines",
			1., 0b0110, "azimuthal", "gnomonic") {
		public double[] project(double lat, double lon, double[] params) {
			if (lat <= 0)	lat = 1e-5;
			final double r = Math.tan(Math.PI/2 - lat);
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon)};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { Math.PI/2 - Math.atan(2*Math.hypot(x, y)),
					Math.atan2(y, x) + Math.PI/2 };
		}
	},
	
	LAMBERT_CONIC("Conformal Conic", 0b0111, "conformal") {
		public double[] project(double lat, double lon, double[] params) {
			final double lat1 = Math.toRadians(params[0]);
			final double lat2 = Math.toRadians(params[1]);
			final double n;
			if (Math.abs(lat1) != Math.abs(lat2))
				n = Math.log(Math.cos(lat1)/Math.cos(lat2))/Math.log(Math.tan(Math.PI/4+lat2/2)/Math.tan(Math.PI/4+lat1/2));
			else if (lat1 == lat2 && lat1 != 0)
				n = Math.sin(lat1);
			else
				return MERCATOR.project(lat, lon, params);
			final double F = Math.cos(lat1)*Math.pow(Math.tan(Math.PI/4+lat1/2), n)/n;
			final double r0 = F*Math.pow(Math.tan(Math.PI/4+(lat1+lat2)/4), -n);
			final double r = F*Math.pow(Math.tan(Math.PI/4+lat/2), -n);
			return new double[] { r*Math.sin(n*lon), r0 - r*Math.cos(n*lon) };
		}
		public double[] inverse(double x, double y, double[] params) {
			final double lat1 = Math.toRadians(params[0]);
			final double lat2 = Math.toRadians(params[1]);
			final double n;
			if (lat1 + lat2 > 0)
				n = Math.log(Math.cos(lat1)/Math.cos(lat2))/Math.log(Math.tan(Math.PI/4+lat2/2)/Math.tan(Math.PI/4+lat1/2));
			else if (Math.abs(lat1) != Math.abs(lat2))
				return Math2.invArr(this.inverse(-x, -y, new double[] {-params[0],-params[1]}));
			else if (lat1 == lat2 && lat1 != 0)
				n = Math.sin(lat1);
			else
				return MERCATOR.inverse(x, y, params);
			final double F = Math.cos(lat1)*Math.pow(Math.tan(Math.PI/4+lat1/2), n)/n;
			final double r0 = F*Math.pow(Math.tan(Math.PI/4+(lat1+lat2)/4), -n);
			final double r = Math.hypot(x, r0/Math.PI-y)*Math.PI;
			final double phi = 2*Math.atan(Math.pow(F/r, 1/n)) - Math.PI/2;
			final double lam = Math.atan2(x, r0/Math.PI-y)/n;
			if (Math.abs(lam) > Math.PI) 	return null;
			else 							return new double[] {phi, lam};
		}
		public double getAspectRatio(double[] params) {
			return 1; //TODO: implement this
		}
	},
	
	E_D_CONIC("Equidistant Conic", 0b1111, "equidistant") {
		public double[] project(double lat, double lon, double[] params) {
			final double r = 3*Math.PI/5 - 4/5.*lat;
			final double tht = lon/2 - Math.PI/2;
			return new double[] { r*Math.cos(tht), r*Math.sin(tht) + Math.PI/2 };
		}
		public double[] inverse(double x, double y, double[] params) {
			y = (y-1)/2;
			final double r = Math.hypot(x, y);
			if (r < 0.2 || r > 1)	return null;
			return new double[] {
					3*Math.PI/4 - 5*Math.PI/4*r, 2*Math.atan2(y, x)+Math.PI };
		}
		public double getAspectRatio(double[] params) {
			return 2; //TODO: implement this
		}
	},
	
	ALBERS("Albers", 0b1111, "equal-area") {
		public double[] project(double lat, double lon, double[] params) {
			final double r = 2*Math.sqrt(1.2 - Math.sin(lat));
			final double tht = lon/2 - Math.PI/2;
			return new double[] { r*Math.cos(tht), r*Math.sin(tht) + Math.PI/2 };
		}
		public double[] inverse(double x, double y, double[] params) {
			y = (y-1)/2;
			final double r = Math.hypot(x, y);
			if (r < Math.sqrt(1/11.) || r > 1)	return null;
			return new double[] {
					Math.asin(1.2-2.2*Math.pow(r, 2)), 2*Math.atan2(y, x)+Math.PI };
		}
		public double getAspectRatio(double[] params) {
			return 2; //TODO: implement this
		}
	},
	
	LEE("Lee", Math.sqrt(3), 0b1001, "tetrahedral", "conformal", null, "that really deserves more attention") {
		public double[] project(double lat, double lon, double[] params) {
			return tetrahedralProjectionForward(lat, lon, (coordR) -> {
				final mfc.field.Complex z = mfc.field.Complex.fromPolar(
						Math.pow(2, 5/6.)*Math.tan(Math.PI/4-coordR[0]/2), coordR[1]);
				final mfc.field.Complex w = Dixon.invFunc(z);
				return new double[] { w.abs()*1.186, w.arg() };
			});
		}
		public double[] inverse(double x, double y, double[] params) {
			final double[] doubles = tetrahedralProjectionInverse(x,y);
			final double[] faceCenter = { doubles[0], doubles[1], doubles[2] };
			final double tht = doubles[3], xp = doubles[4], yp = doubles[5];
			
			final mfc.field.Complex w = mfc.field.Complex.fromPolar(
					Math.hypot(xp, yp)*1.53,
					Math.atan2(yp, xp)+tht - Math.PI/2);
			final mfc.field.Complex ans = Dixon.leeFunc(w).times(Math.pow(2, -5/6.));
			final double[] triCoords = {
					Math.PI/2 - 2*Math.atan(ans.abs()),
					ans.arg() + Math.PI };
			
			return obliquifyPlnr(triCoords, faceCenter);
		}
	},
	
	AUTHAGRAPH("AuthaGraph", "A hip new Japanese map that is almost authagraphic (this is an approximation; they won't give me their actual equations)",
			4/Math.sqrt(3), 0b1001, "tetrahedral", "compromise") {
		public double[] project(double lat, double lon, double[] params) {
			return null;
		}
		public double[] inverse(double x, double y, double[] params) {
			final double[] faceCenter = new double[3];
			final double rot, localX, localY;
			if (y-1 < 4*x && y-1 < -4*x) {
				faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
				faceCenter[1] = 0;
				rot = 0;
				localX = 4/Math.sqrt(3)*x;
				localY = y+1/3.0;
			}
			else if (y-1 < -4*(x+1)) {
				faceCenter[0] = -Math.PI/2;
				faceCenter[1] = Math.PI;
				rot = 0;
				localX = 4/Math.sqrt(3)*(x+1);
				localY = y+1/3.0;
			}
			else if (y-1 < 4*(x-1)) {
				faceCenter[0] = -Math.PI/2;
				faceCenter[1] = Math.PI;
				rot = 0;
				localX = 4/Math.sqrt(3)*(x-1);
				localY = y+1/3.0;
			}
			else if (x < 0) {
				faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
				faceCenter[1] = 4*Math.PI/3;
				rot = Math.PI/3;
				localX = 4/Math.sqrt(3)*(x+0.5);
				localY = y-1/3.0;
			}
			else {
				faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
				faceCenter[1] = 2*Math.PI/3;
				rot = -Math.PI/3;
				localX = 4/Math.sqrt(3)*(x-0.5);
				localY = y-1/3.0;
			}
			faceCenter[2] = 0;
			
			final double t = Math.atan2(localY, localX) + rot;
			final double t0 = Math.floor((t+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double dt = t-t0;
			final double z = 2.49*Math.hypot(localX, localY)*Math.cos(dt);
			final double g = 0.03575*z*z*z + 0.0219*z*z + 0.4441*z;
			double[] triCoords = {
					Math.PI/2 - Math.atan(Math.tan(g)/Math.cos(dt)),
					Math.PI/2 + t0 + dt };
			return obliquifyPlnr(triCoords, faceCenter);
		}
	},
	
	SINUSOIDAL("Sinusoidal", "An equal-area map shaped like a sine-wave",
			2., 0b1111, "pseudocylindrical", "equal-area") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] { Math.cos(lat)*lon, lat };
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { y*Math.PI/2, x*Math.PI/Math.cos(y*Math.PI/2) };
		}
	},
	
	MOLLWEIDE("Mollweide", "An equal-area projection shaped like an ellipse",
			2., 0b1101, "pseudocylindrical", "equal-area") {
		public double[] project(double lat, double lon, double[] params) {
			double tht = lat;
			for (int i = 0; i < 10; i ++)
				tht -= (2*tht+Math.sin(2*tht)-Math.PI*Math.sin(lat))/
						(2+2*Math.cos(2*tht));
			return new double[] { lon*Math.cos(tht), Math.PI/2*Math.sin(tht) };
		}
		public double[] inverse(double x, double y, double[] params) {
			double tht = Math.asin(y);
			return new double[] {
					Math.asin((2*tht + Math.sin(2*tht)) / Math.PI),
					Math.PI * x / Math.cos(tht)};
		}
	},
	
	AITOFF("Aitoff", "A compromise projection shaped like an ellipse",
			2., 0b1011, "pseudoazimuthal", "equal-area") {
		public double[] project(double lat, double lon, double[] params) {
			final double a = Math.acos(Math.cos(lat)*Math.cos(lon/2));
			return new double[] {
					2*Math.cos(lat)*Math.sin(lon/2)*a/Math.sin(a),
					Math.sin(lat)*a/Math.sin(a)};
		}
		public double[] inverse(double x, double y, double[] params) {
			final double[] intermediate = POLAR.inverse(x/2, y/2, params);
			double[] transverse = obliquifyPlnr(intermediate, new double[] {0,0,0});
			if (transverse != null) 	transverse[1] *= 2;
			return transverse;
		}
	},
	
	HAMMER("Hammer", "An equal-area projection shaped like an ellipse",
			2., 0b1111, "pseudoazimuthal", "equal-area") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] {
					Math.PI*Math.cos(lat)*Math.sin(lon/2)/Math.sqrt(1+Math.cos(lat)*Math.cos(lon/2)),
					Math.PI/2*Math.sin(lat)/Math.sqrt(1+Math.cos(lat)*Math.cos(lon/2)) };
		}
		public double[] inverse(double x, double y, double[] params) {
			final double X = x * Math.sqrt(8);
			final double Y = y * Math.sqrt(2);
			final double z = Math.sqrt(1 - Math.pow(X/4, 2) - Math.pow(Y/2, 2));
			return new double[] {
					Math.asin(z * Y), 2*Math.atan(0.5*z*X / (2*z*z - 1))};
		}
	},
	
	TOBLER("Tobler", "An equal-area projection shaped like a hyperellipse (in case you're wondering about gamma, it's calculated automatically)",
			2., 0b1001, "pseudocylindrical", "equal-area",new String[]{"Std. Parallel","alpha","K"},
			new double[][] {{0,89,37.5}, {0,1,0}, {1,8,2.5}}) {
		private double[] Z; //Z[i] = sin(phi) when y = i/(Z.length-1)
		private double[] lastParams;
		private void generateZ(int n, double[] params) {
			final double e = NumericalAnalysis.simpsonIntegrate(0, 1,
					Tobler::hyperEllipse, .0001, params);
			Z = NumericalAnalysis.simpsonODESolve(1, n,
					Tobler::dZdY, Math.min(1./n,.0001),
					e, params[1], params[2]);
			lastParams = params.clone();
		}
		public double[] project(double lat, double lon, double[] params) {
			if (Z == null || !Arrays.equals(params, lastParams))
				generateZ(10000, params);
			final double z0 = Math.abs(Math.sin(lat));
			final int i = Arrays.binarySearch(Z, z0);
			final double y;
			if (i >= 0)
				y = i/(Z.length-1.);
			else if (-i-1 >= Z.length)
				y = Z[Z.length-1];
			else
				y = Math2.linInterp(z0, Z[-i-2], Z[-i-1], -i-2, -i-1)/
						(Z.length-1.);
			final double ar = Math.pow(Math.cos(Math.toRadians(params[0])),2);
			return new double[] {
					Tobler.X(y, lon, params), y*Math.signum(lat)/ar };
		}
		public double[] inverse(double x, double y, double[] params) {
			if (Z == null || !Arrays.equals(params, lastParams))
				generateZ(10000, params);
			return new double[] {
					Math.asin(Z[(int)Math.round(Math.abs(y)*(Z.length-1))])*Math.signum(y),
					Tobler.lam(x,y,params) };
		}
		public double getAspectRatio(double[] params) {
			return Math.pow(Math.cos(Math.toRadians(params[0])),2) * Math.PI;
		}
	},
	
	VAN_DER_GRINTEN("Van der Grinten", "A circular compromise map that is popular for some reason",
			1., 0b1111, "other", "compromise") {
		public double[] project(double lat, double lon, double[] params) {
			final double t = Math.asin(Math.abs(2*lat/Math.PI));
			if (lat == 0) // special case 1: equator
				return new double[] {lon, 0};
			if (lon == 0 || lat >= Math.PI/2 || lat <= -Math.PI/2) // special case 3: meridian
				return new double[] {0, Math.signum(lat)*Math.PI*Math.tan(t/2)};
			final double A = Math.abs(Math.PI/lon - lon/Math.PI)/2;
			final double G = Math.cos(t)/(Math.sin(t)+Math.cos(t)-1);
			final double P = G*(2/Math.sin(t) - 1);
			final double Q = A*A + G;
			return new double[] {
					Math.PI*Math.signum(lon)*(A*(G-P*P)+Math.sqrt(A*A*(G-P*P)*(G-P*P)-(P*P+A*A)*(G*G-P*P)))/(P*P+A*A),
					Math.PI*Math.signum(lat)*(P*Q-A*Math.sqrt((A*A+1)*(P*P+A*A)-Q*Q))/(P*P+A*A)};
		}
		public double[] inverse(double x, double y, double[] params) {
			if (y == 0) // special case 1: equator
				return new double[] {0, x*Math.PI};
			if (x == 0) // special case 3: meridian
				return new double[] {Math.PI/2 * Math.sin(2*Math.atan(y)), 0};
			
			double c1 = -Math.abs(y) * (1 + x*x + y*y);
			double c2 = c1 - 2*y*y + x*x;
			double c3 = -2 * c1 + 1 + 2*y*y + Math.pow(x*x + y*y, 2);
			double d = y*y / c3 + 1 / 27.0 * (2*Math.pow(c2 / c3, 3) - 9*c1*c2 / (c3*c3));
			double a1 = 1 / c3*(c1 - c2*c2 / (3*c3));
			double m1 = 2 * Math.sqrt(-a1 / 3);
			double tht1 = Math.acos(3*d / (a1 * m1)) / 3;
			return new double[] {
					Math.signum(y) * Math.PI * (-m1 * Math.cos(tht1 + Math.PI/3) - c2 / (3*c3)),
					Math.PI*(x*x + y*y - 1 + Math.sqrt(1 + 2*(x*x - y*y) + Math.pow(x*x + y*y, 2)))
							/ (2*x)};
		}
	},
	
	ROBINSON("Robinson", "A visually pleasing piecewise compromise map",
			1.9716, 0b1111, "pseudocylindrical", "compromise") {
		public double[] project(double lat, double lon, double[] params) {
			return new double[] { Robinson.plenFromLat(lat)*lon,
								Robinson.pdfeFromLat(lat)*Math.PI/2};
		}
		public double[] inverse(double x, double y, double[] params) {
			return new double[] { Robinson.latFromPdfe(y),
								x/Robinson.plenFromPdfe(y)*Math.PI };
		}
	},
	
	WINKEL_TRIPEL("Winkel Tripel", "National Geographic's compromise projection of choice",
			0, 0b1011, "other", "compromise",
			new String[] {"Std. Parallel"}, new double[][] {{0,90,50.4598}}) {
		public double[] project(double lat, double lon, double[] params) {
			final double cosphi0 = Math.cos(Math.toRadians(params[0]));
			return new double[] { WinkelTripel.f1pX(lat,lon,cosphi0)/(1+cosphi0),
					WinkelTripel.f2pY(lat,lon,cosphi0)/(1+cosphi0) };
		}
		public double[] inverse(double x, double y, double[] params) {
			final double cosphi0 = Math.cos(Math.toRadians(params[0]));
			return NumericalAnalysis.newtonRaphsonApproximation(
					x*Math.PI*(1 + cosphi0), y*Math.PI,
					y*Math.PI/2, x*Math.PI*(1 + Math.cos(y*Math.PI/2))/2,
					WinkelTripel::f1pX, WinkelTripel::f2pY,
					WinkelTripel::df1dphi, WinkelTripel::df1dlam,
					WinkelTripel::df2dphi, WinkelTripel::df2dlam, .002, cosphi0);
		}
		public double getAspectRatio(double[] params) {
			return 1 + Math.cos(Math.toRadians(params[0]));
		}
	},
	
	PEIRCE_QUINCUNCIAL("Peirce Quincuncial", "A conformal projection that uses complex elliptic functions",
			1., 0b1001, "other", "conformal") {
		public double[] project(double lat, double lon, double[] params) {
			final double alat = Math.abs(lat);
			final double wMag = Math.tan(Math.PI/4-alat/2);
			final Complex w = new Complex(wMag*Math.sin(lon), -wMag*Math.cos(lon));
			final Complex k = new Complex(Math.sqrt(0.5));
			Complex z = Elliptic.F(w.acos(),k).multiply(Math.PI/1.854).subtract(Math.PI).negate();
			if (z.isInfinite() || z.isNaN())	z = new Complex(0);
			double x = z.getReal(), y = z.getImaginary();
			
			if (lat < 0) {
				if (x >= 0 && y >= 0)
					z = new Complex(Math.PI-y, Math.PI-x);
				else if (x >= 0 && y < 0)
					z = new Complex(Math.PI+y, -Math.PI+x);
				else if (y >= 0)
					z = new Complex(-Math.PI+y, Math.PI+x);
				else
					z = new Complex(-Math.PI-y, -Math.PI-x);
			}
			return new double[] {z.getReal(), z.getImaginary()};
		}
		public double[] inverse(double x, double y, double[] params) {
			mfc.field.Complex u = new mfc.field.Complex(1.854*(x+1), 1.854*y); // 1.854 is approx K(sqrt(1/2)
			mfc.field.Complex k = new mfc.field.Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus
			mfc.field.Complex ans = Jacobi.cn(u, k);
			double p = 2 * Math.atan(ans.abs());
			double theta = ans.arg() - Math.PI/2;
			double lambda = Math.PI/2 - p;
			return new double[] {lambda, theta};
		}
	},
	
	GUYOU("Guyou", "Peirce Quincuncial, rearranged a bit", 2., 0b1001, "other", "conformal") {
		private final double[] POLE = {0, -Math.PI/2, Math.PI/4};
		
		public double[] project(double lat, double lon, double[] params) {
			final double[] coords = obliquifySphc(lat,lon, POLE);
			final double alat = Math.abs(coords[0]);
			final double wMag = Math.tan(Math.PI/4-alat/2);
			final Complex w = new Complex(wMag*Math.sin(coords[1]), -wMag*Math.cos(coords[1]));
			final Complex k = new Complex(Math.sqrt(0.5));
			Complex z = Elliptic.F(w.acos(),k).multiply(new Complex(Math.PI/3.708,Math.PI/3.708)).subtract(new Complex(0,Math.PI/2));
			if (z.isInfinite() || z.isNaN()) 	z = new Complex(0);
			if (coords[0] < 0) 	z = z.conjugate().negate();
			return new double[] {z.getReal(), z.getImaginary()};
		}
		public double[] inverse(double x, double y, double[] params) {
			mfc.field.Complex u = new mfc.field.Complex(1.8558*(x - y/2 - 0.5), 1.8558*(x + y/2 + 0.5)); // don't ask me where 3.7116 comes from
			mfc.field.Complex k = new mfc.field.Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus
			mfc.field.Complex ans = Jacobi.cn(u, k);
			double p = 2 * Math.atan(ans.abs());
			double theta = ans.arg();
			double lambda = Math.PI/2 - p;
			return obliquifyPlnr(new double[] {lambda,theta}, POLE);
		}
	},
	
	LEMONS("Lemons", "BURN LIFE'S HOUSE DOWN!!!", 2., 0b1110, "pseudocylindrical", "?") {
		public double[] project(double lat, double lon, double[] params) {
			return null;
		}
		public double[] inverse(double x, double y, double[] params) {
			x = x+2;
			final double lemWdt = 1/6.0;
			if (Math.abs(x % lemWdt - lemWdt / 2.0) <= Math.cos(y*Math.PI/2) * lemWdt/2.0) // if it is in
				return new double[] { y*Math.PI/2,	// a sine curve
						Math.PI * (x%lemWdt - lemWdt/2.0) / (Math.cos(y*Math.PI/2))
								+ (int)(x/lemWdt) * Math.PI/6 };
			else
				return null;
		}
	},
	
	MAGNIFIER("Magnifier", "A novelty map projection that blows up the center way out of proportion",
			1., 0b1011, "azimuthal", "pointless") {
		public double[] project(double lat, double lon, double[] params) {
			final double p = 1/2.0+lat/Math.PI;
			final double fp = 1 - 0.1*p - 0.9*Math.pow(p,7);
			final double r = Math.PI*fp;
			return new double[] { r*Math.sin(lon), -r*Math.cos(lon) };
		}
		public double[] inverse(double x, double y, double[] params) {
			double R = Math.hypot(x, y);
			if (R <= 1)
				return new double[] {
						Math.PI/2 * (1 - R*.2 - R*R*R*1.8),
						Math.atan2(y, x) + Math.PI/2};
			else
				return null;
		}
	},
	
	EXPERIMENT("Experiment", "What happens when you apply a complex differentiable function to a stereographic projection?",
			1., 0b0000, "?", "conformal") {
		public double[] project(double lat, double lon, double[] params) {
			final double wMag = Math.tan(Math.PI/4-lat/2);
			final Complex w = new Complex(wMag*Math.sin(lon), -wMag*Math.cos(lon));
			Complex z = w.asin();
			if (z.isInfinite() || z.isNaN())	z = new Complex(0);
			return new double[] { z.getReal(), z.getImaginary() };
		}
		public double[] inverse(double x, double y, double[] params) {
			Complex z = new Complex(x*3, y*3);
			Complex ans = z.sin();
			double p = 2 * Math.atan(ans.abs());
			double theta = ans.getArgument();
			double lambda = Math.PI/2 - p;
			return new double[] {lambda, theta};
		}
	},
	
	TETRAGRAPH("TetraGraph", Math.sqrt(3), 0b1111, "tetrahedral", "equidistant", null, "that I invented") {
		public double[] project(double lat, double lon, double[] params) {
			return tetrahedralProjectionForward(lat, lon, (coordR) -> {
				final double tht = coordR[1] - Math.floor(coordR[1]/(2*Math.PI/3))*(2*Math.PI/3) - Math.PI/3;
				return new double[] {
						Math.atan(1/Math.tan(coordR[0])*Math.cos(tht))/Math.cos(tht)*Math.PI/3/Math.atan(Math.sqrt(2)),
						coordR[1]
				};
			});
		}
		public double[] inverse(double x, double y, double[] params) {
			final double[] doubles = tetrahedralProjectionInverse(x,y);
			final double[] faceCenter = { doubles[0], doubles[1], doubles[2] };
			final double tht = doubles[3], xp = doubles[4], yp = doubles[5];
			final double t = Math.atan2(yp, xp) + tht;
			final double t0 = Math.floor((t+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double dt = t-t0;
			double[] triCoords = {
					Math.PI/2 - Math.atan(Math.tan(Math.atan(Math.sqrt(2))/(Math.sqrt(3)/3)*Math.hypot(xp,yp)*Math.cos(dt))/Math.cos(dt)),
					Math.PI/2 + t0 + dt};
			return obliquifyPlnr(triCoords, faceCenter);
		}
	},
	
	HYPERELLIPOWER("Hyperellipower", "A parametric projection that I'm still testing",
			2., 0b1111, "pseudocylindrical", "compromise", new String[] {"k","n","a"},
			new double[][] {{1,5,3.7308},{.5,2.,1.2027},{.5,2.,1.1443}}) {
		public double[] project(double lat, double lon, double[] params) {
			final double k = params[0], n = params[1], a = params[2];
			final double ynorm = (1-Math.pow(1-Math.abs(lat/(Math.PI/2)), n));
			return new double[] {
					Math.pow(1 - Math.pow(ynorm, k),1/k)*lon,
					ynorm*Math.PI/2/Math.sqrt(n)*a*Math.signum(lat)
				};
		}
		public double[] inverse(double x, double y, double[] params) {
			final double k = params[0], n = params[1];
			return new double[] {
					(1 - Math.pow(1-Math.abs(y), 1/n))*Math.PI/2*Math.signum(y),
					x/Math.pow(1 - Math.pow(Math.abs(y),k),1/k)*Math.PI };
		}
		public double getAspectRatio(double[] params) {
			final double n = params[1], a = params[2];
			return 2*Math.sqrt(n)/a;
		}
	},
	
	TETRAPOWER("Tetrapower", "A parametric projection that I'm still testing",
			Math.sqrt(3), 0b1111, "tetrahedral", "compromise", new String[] {"k1","k2","k3"},
			new double[][] {{.25,4.,1.4586},{.25,4.,1.2891},{.25,4.,1.9583}}) {
		public double[] project(double lat, double lon, double[] params) {
			final double k1 = params[0], k2 = params[1], k3 = params[2];
			return tetrahedralProjectionForward(lat, lon, (coordR) -> {
				final double t0 = Math.floor(coordR[1]/(2*Math.PI/3))*(2*Math.PI/3) + Math.PI/3;
				final double tht = coordR[1] - t0;
				final double thtP = Math.PI/3*(1 - Math.pow(1-Math.abs(tht)/(Math.PI/2),k1))/(1 - 1/Math.pow(3,k1))*Math.signum(tht);
				final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
				final double rmax = .5/Math.cos(thtP); //the max normalized radius of this triangle (in the plane)
				final double rtgf = Math.atan(1/Math.tan(coordR[0])*Math.cos(tht))/Math.atan(Math.sqrt(2))*rmax;
				return new double[] {
						(1 - Math.pow(1-rtgf,kRad))/(1 - Math.pow(1-rmax,kRad))*rmax*2*Math.PI/3,
						thtP + t0 };
			});
		}
		public double[] inverse(double x, double y, double[] params) {
			final double k1 = params[0], k2 = params[1], k3 = params[2];
			final double[] doubles = tetrahedralProjectionInverse(x,y);
			final double[] faceCenter = { doubles[0], doubles[1], doubles[2] };
			final double tht = doubles[3], xp = doubles[4], yp = doubles[5];
			final double R = Math.hypot(xp, yp)*Math.sqrt(3)/2;
			final double t = Math.atan2(yp, xp) + tht;
			final double t0 = Math.floor((t+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double thtP = t-t0;
			final double lamS = (1-Math.pow(1-Math.abs(thtP)*(1-1/Math.pow(3,k1))/(Math.PI/3), 1/k1))*Math.PI/2*Math.signum(thtP);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax = .5/Math.cos(thtP); //the max normalized radius of this triangle (in the plane)
			final double rtgf = 1-Math.pow(1-R/rmax*(1-Math.pow(Math.abs(1-rmax), kRad)), 1/kRad); //normalized tetragraph radius
			double[] triCoords = {
					Math.atan(Math.cos(lamS)/Math.tan(rtgf/rmax*Math.atan(Math.sqrt(2)))),
					Math.PI/2 + t0 + lamS };
			return obliquifyPlnr(triCoords, faceCenter);
		}
	},
	
	TETRAFILLET("Tetrafillet", "A parametric projection that I'm still testing",
			Math.sqrt(3), 0b1111, "other", "compromise", new String[] {"k1","k2","k3"},
			new double[][] {{.25,4.,1.1598},{.25,4.,.36295},{.25,4.,1.9553}}) {
		public double[] project(double lat, double lon, double[] params) {
			final double k1 = params[0], k2 = params[1], k3 = params[2];
			return tetrahedralProjectionForward(lat, lon, (coordR) -> {
				final double t0 = Math.floor(coordR[1]/(2*Math.PI/3))*(2*Math.PI/3) + Math.PI/3;
				final double tht = coordR[1] - t0;
				final double thtP = Math.PI/3*(1 - Math.pow(1-Math.abs(tht)/(Math.PI/2),k1))/(1 - 1/Math.pow(3,k1))*Math.signum(tht);
				final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
				final double rmax; //the max normalized radius of this triangle (in the plane)
				if (Math.abs(thtP) < .70123892) 	rmax = .5/Math.cos(thtP);
				else 	rmax = .75 - 1.5972774*Math.pow(Math.PI/3-Math.abs(thtP),2)/2;
				final double rtgf = Math.atan(1/Math.tan(coordR[0])*Math.cos(tht))/Math.atan(Math.sqrt(2))*rmax; //normalized tetragraph radius
				return new double[] {
						(1 - Math.pow(1-rtgf,kRad))/(1 - Math.pow(1-rmax,kRad))*rmax*2*Math.PI/3,
						thtP + t0
				};
			});
		}
		public double[] inverse(double x, double y, double[] params) {
			final double k1 = params[0], k2 = params[1], k3 = params[2];
			final double[] doubles = tetrahedralProjectionInverse(x,y);
			final double[] faceCenter = { doubles[0], doubles[1], doubles[2] };
			final double tht = doubles[3], xp = doubles[4], yp = doubles[5];
			final double R = Math.hypot(xp, yp)*Math.sqrt(3)/2;
			final double t = Math.atan2(yp, xp) + tht;
			final double t0 = Math.floor((t+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double thtP = t-t0;
			final double lamS = (1-Math.pow(1-Math.abs(thtP)*(1-1/Math.pow(3,k1))/(Math.PI/3), 1/k1))*Math.PI/2*Math.signum(thtP);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax; //the max normalized radius of this triangle (in the plane)
			if (Math.abs(thtP) < .70123892) 	rmax = .5/Math.cos(thtP);
			else 	rmax = .75 - 1.5972774*Math.pow(Math.PI/3-Math.abs(thtP),2)/2;
			final double rtgf = 1-Math.pow(1-R/rmax*(1-Math.pow(Math.abs(1-rmax), kRad)), 1/kRad); //normalized tetragraph radius
			if (R > rmax) 	return null;
			double[] triCoords = {
					Math.atan(Math.cos(lamS)/Math.tan(rtgf/rmax*Math.atan(Math.sqrt(2)))),
					Math.PI/2 + t0 + lamS };
			return obliquifyPlnr(triCoords, faceCenter);
		}
	},
	
	TETRACHAMFER("Tetrachamfer", "A parametric projection that I'm still testing",
			Math.sqrt(3), 0b1111, "other", "compromise", new String[] {"k1","k2","k3"},
			new double[][] {{.25,4.,1.1598},{.25,4.,.36295},{.25,4.,1.9553}}) {
		public double[] project(double lat, double lon, double[] params) {
			final double k1 = params[0], k2 = params[1], k3 = params[2];
			return tetrahedralProjectionForward(lat, lon, (coordR) -> {
				final double t0 = Math.floor(coordR[1]/(2*Math.PI/3))*(2*Math.PI/3) + Math.PI/3;
				final double tht = coordR[1] - t0;
				final double thtP = Math.PI/3*(1 - Math.pow(1-Math.abs(tht)/(Math.PI/2),k1))/(1 - 1/Math.pow(3,k1))*Math.signum(tht);
				final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
				final double rmax = Math.min(.5/Math.cos(thtP), .75/Math.cos(Math.PI/3-Math.abs(thtP))); //the max normalized radius of this triangle (in the plane)
				final double rtgf = Math.atan(1/Math.tan(coordR[0])*Math.cos(tht))/Math.atan(Math.sqrt(2))*rmax; //normalized tetragraph radius
				return new double[] {
						(1 - Math.pow(1-rtgf,kRad))/(1 - Math.pow(1-rmax,kRad))*rmax*2*Math.PI/3,
						thtP + t0
				};
			});
		}
		public double[] inverse(double x, double y, double[] params) {
			final double k1 = params[0], k2 = params[1], k3 = params[2];
			final double[] doubles = tetrahedralProjectionInverse(x,y);
			final double[] faceCenter = { doubles[0], doubles[1], doubles[2] };
			final double tht = doubles[3], xp = doubles[4], yp = doubles[5];
			final double R = Math.hypot(xp, yp)*Math.sqrt(3)/2;
			final double t = Math.atan2(yp, xp) + tht;
			final double t0 = Math.floor((t+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double thtP = t-t0;
			final double lamS = (1-Math.pow(1-Math.abs(thtP)*(1-1/Math.pow(3,k1))/(Math.PI/3), 1/k1))*Math.PI/2*Math.signum(thtP);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax = Math.min(.5/Math.cos(thtP), .75/Math.cos(Math.PI/3-Math.abs(thtP))); //the max normalized radius of this triangle (in the plane)
			final double rtgf = 1-Math.pow(1-R/rmax*(1-Math.pow(Math.abs(1-rmax), kRad)), 1/kRad); //normalized tetragraph radius
			if (R > rmax) 	return null;
			double[] triCoords = {
					Math.atan(Math.cos(lamS)/Math.tan(rtgf/rmax*Math.atan(Math.sqrt(2)))),
					Math.PI/2 + t0 + lamS };
			return obliquifyPlnr(triCoords, faceCenter);
		}
	};
	
	
	
	private String name;
	private String description;
	private String[] paramNames;
	private double[][] paramValues;
	
	private double aspectRatio; //default W/H for the map
	private boolean finite; //is it completely bounded?
	private boolean invertable; //is the inverse solution closed-form?
	private boolean solveable; //is the solution closed-form?
	private boolean continuous; //can you see the whole earth without inerruption?
	private String type; //cylindrical, azimuthal, etc.
	private String property; //what it is good for
	
	
	
	private Projection(String name, int fisc, String property) { //this one is just for conic maps, because they're so similar
		this(name, buildDescription("conic",property,null,null), 0., fisc, "conic", property,
				new String[] {"Std. Parallel 1", "Std. Parallel 2"},
				new double[][] {{-89,89,15}, {-89,89,45}});
	}
	
	private Projection(String name, double aspectRatio, int fisc, String type, String property) {
		this(name, buildDescription(type,property,null,null), aspectRatio, fisc, type, property, new String[0], new double[0][]);
	}
	
	private Projection(String name, double aspectRatio, int fisc, String type, String property, String adjective) {
		this(name, buildDescription(type,property,adjective,null), aspectRatio, fisc, type, property, new String[0], new double[0][]);
	}
	
	private Projection(String name, double aspectRatio, int fisc, String type, String property, String adjective, String addendum) {
		this(name, buildDescription(type,property,adjective,addendum), aspectRatio, fisc, type, property, new String[0], new double[0][]);
	}
	
	private Projection(String name, String description, double aspectRatio, int fisc, String type, String property) {
		this(name, description, aspectRatio, fisc, type, property, new String[0], new double[0][]);
	}
	
	private Projection(String name, String description, double aspectRatio,
			int fisc, String type, String property, String[] paramNames, double[][] paramValues) {
		this.name = name;
		this.description = description;
		this.paramNames = paramNames;
		this.paramValues = paramValues;
		this.aspectRatio = aspectRatio;
		this.finite = (fisc&0b1000) > 0;
		this.invertable = (fisc&0b0100) > 0;
		this.solveable = (fisc&0b0010) > 0;
		this.continuous = (fisc&0b0001) > 0;
		this.type = type;
		this.property = property;
	}
	
	
	private static String buildDescription(String type, String property, String adjective, String addendum) { //these should all be lowercase
		String description = property+" "+type+" projection";
		if (adjective != null)
			description = adjective+" "+description;
		if (addendum != null)
			description += " "+addendum;
		if (description.charAt(0) == 'a' || description.charAt(0) == 'e' || description.charAt(0) == 'i' || description.charAt(0) == 'o' || description.charAt(0) == 'u')
			return "An "+description;
		else
			return "A "+description;
	}
	
	
	
	public abstract double[] project(double lat, double lon, double[] params);
	
	public abstract double[] inverse(double x, double y, double[] params);
	
	
	public double[] project(double[] coords, double[] params) {
		return project(coords[0], coords[1], params);
	}
	
	public double[] project(double[] coords, double[] params, double[] pole) {
		return project(coords[0], coords[1], params, pole);
	}
	
	public double[] project(double lat, double lon, double[] params, double[] pole) {
		return project(obliquifySphc(lat, lon, pole), params);
	}
	
	
	public double[] inverse(double[] coords, double[] params) {
		return inverse(coords[0], coords[1], params);
	}
	
	public double[] inverse(double[] coords, double[] params, double[] pole) {
		return inverse(coords[0], coords[1], params, pole);
	}
	
	public double[] inverse(double x, double y, double[] params, double[] pole) {
		return obliquifyPlnr(inverse(x, y, params), pole);
	}
	
	
	public List<List<double[]>> transform(List<List<double[]>> curves, int step,
			double[] params, double[] pole, DoubleConsumer tracker) {
		List<List<double[]>> output = new LinkedList<List<double[]>>();
		int i = 0;
		for (List<double[]> curve0: curves) {
			if (curve0.size() < step*3)	continue;
			
			List<double[]> curve1 = new ArrayList<double[]>(curve0.size()/step);
			for (int j = 0; j < curve0.size(); j += step)
				curve1.add(project(curve0.get(j), params, pole));
			output.add(curve1);
			
			if (tracker != null) {
				i ++;
				tracker.accept((double)i/curves.size());
			}
		}
		return output;
	}
	
	
	public double[][][] map(int size, double[] params) {
		return map(size, params, new double[] {Math.PI/2,0,0});
	}
	
	public double[][][] map(int size, double[] params, double[] pole) {
		final double ratio = this.getAspectRatio(params);
		if (ratio < 1)
			return map(Math.max(Math.round(size*ratio),1), size, params, pole, null);
		else
			return map(size, Math.max(Math.round(size/ratio),1), params, pole, null);
	}
	
	public double[][][] map(double w, double h, double[] params, double[] pole, DoubleConsumer tracker) { //generate a matrix of coordinates based on a map projection
		double[][][] output = new double[(int) h][(int) w][2]; //the coordinate matrix
		
		for (int y = 0; y < output.length; y ++) {
			for (int x = 0; x < output[y].length; x ++)
				output[y][x] = inverse((2*x+1)/w-1, 1-(2*y+1)/h, params, pole);
			if (tracker != null) 	tracker.accept((double) y/output.length);
		}
		
		return output;
	}
	
	
	public static double[][][] globe(double dt) { //generate a matrix of coordinates based on the sphere
		List<double[]> points = new ArrayList<double[]>();
		for (double phi = -Math.PI/2+dt/2; phi < Math.PI/2; phi += dt) { // make sure phi is never exactly +-tau/4
			for (double lam = -Math.PI; lam < Math.PI; lam += dt/Math.cos(phi)) {
				points.add(new double[] {phi, lam});
			}
		}
		return new double[][][] {points.toArray(new double[0][])};
	}
	
	
	public double[] avgDistortion(double[][][] points, double[] params) {
		final double[][][] distDist = calculateDistortion(points, params);
		return new double[] {Math2.stdDev(distDist[0]), Math2.mean(distDist[1])};
	}
	
	
	public double[][][] calculateDistortion(double[][][] points, double[] params) {
		return calculateDistortion(points, params, null);
	}
	
	public double[][][] calculateDistortion(double[][][] points,
			double[] params, ProgressBarDialog pBar) { //calculate both kinds of distortion over the given region
		double[][][] output = new double[2][points.length][points[0].length]; //the distortion matrix
		
		for (int y = 0; y < points.length; y ++) {
			for (int x = 0; x < points[y].length; x ++) {
				if (points[y][x] != null) {
					final double[] dists = getDistortionAt(points[y][x], params);
					output[0][y][x] = dists[0]; //the output matrix has two layers:
					output[1][y][x] = dists[1]; //area and angular distortion
				}
				else {
					output[0][y][x] = Double.NaN;
					output[1][y][x] = Double.NaN; //NaN means no map here
				}
			}
			if (pBar != null)
				pBar.setProgress((double)(y+1)/points.length);
		}
		
		final double avgArea = Math2.mean(output[0]); //don't forget to normalize output[0] so the average is zero
		for (int y = 0; y < output[0].length; y ++)
			for (int x = 0; x < output[0][y].length; x ++)
				output[0][y][x] -= avgArea;
		
		return output;
	}
	
	
	public double[] getDistortionAt(double[] s0, double[] params) { //calculate both kinds of distortion at the given point
		final double[] output = new double[2];
		final double dx = 1e-6;
		
		final double[] s1 = { s0[0], s0[1]+dx/Math.cos(s0[0]) }; //consider a point slightly to the east
		final double[] s2 = { s0[0]+dx, s0[1] }; //and slightly to the north
		final double[] p0 = project(s0, params);
		final double[] p1 = project(s1, params);
		final double[] p2 = project(s2, params);
		
		final double dA = 
				(p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]);
		output[0] = Math.log(Math.abs(dA/(dx*dx))); //the zeroth output is the size (area) distortion
		if (Math.abs(output[0]) > 25)
			output[0] = Double.NaN; //discard outliers
		
		final double s1ps2 = Math.hypot((p1[0]-p0[0])+(p2[1]-p0[1]), (p1[1]-p0[1])-(p2[0]-p0[0]));
		final double s1ms2 = Math.hypot((p1[0]-p0[0])-(p2[1]-p0[1]), (p1[1]-p0[1])+(p2[0]-p0[0]));
		output[1] = Math.abs(Math.log(Math.abs((s1ps2-s1ms2)/(s1ps2+s1ms2))));
		
		return output;
	}
	
	
	private static double[] tetrahedralProjectionForward(double lat, double lon, UnaryOperator<double[]> func) { //a helper function for projections like tetragraph and lee
		final double[][] centrums = {{-Math.PI/2, 0, Math.PI/3},
				{Math.asin(1/3.0), Math.PI, Math.PI/3},
				{Math.asin(1/3.0), Math.PI/3, Math.PI/3},
				{Math.asin(1/3.0), -Math.PI/3, -Math.PI/3}};
		double latR = Double.NaN;
		double lonR = Double.NaN;
		byte poleIdx = -1;
		for (byte i = 0; i < 4; i++) {
			final double[] relCoords = obliquifySphc(lat, lon, centrums[i]);
			if (Double.isNaN(latR) || relCoords[0] > latR) {
				latR = relCoords[0]; // pick the centrum that maxes out your latitude
				lonR = relCoords[1];
				poleIdx = i;
			}
		}
		
		final double[] rtht = func.apply(new double[] { latR, lonR });
		
		switch (poleIdx) {
		case 0:
			if (Math.sin(lon) < 0)
				return new double[]
						{-2*Math.PI/3 + rtht[0]*Math.sin(rtht[1]-Math.PI/6), -Math.PI/Math.sqrt(3) - rtht[0]*Math.cos(rtht[1]-Math.PI/6)}; // lower left
			else
				return new double[]
						{2*Math.PI/3 - rtht[0]*Math.sin(rtht[1]-Math.PI/6), -Math.PI/Math.sqrt(3) + rtht[0]*Math.cos(rtht[1]-Math.PI/6)}; // lower right
		case 1:
			if (Math.sin(lon) < 0)
				return new double[]
						{-2*Math.PI/3 + rtht[0]*Math.sin(rtht[1]-Math.PI/6), Math.PI/Math.sqrt(3) - rtht[0]*Math.cos(rtht[1]-Math.PI/6)}; // upper left
			else
				return new double[]
						{2*Math.PI/3 - rtht[0]*Math.sin(rtht[1]-Math.PI/6), Math.PI/Math.sqrt(3) + rtht[0]*Math.cos(rtht[1]-Math.PI/6)}; // upper right
		case 2:
			return new double[]
					{Math.PI/3 + rtht[0]*Math.cos(rtht[1]), rtht[0]*Math.sin(rtht[1])}; // right
		case 3:
			return new double[]
					{-Math.PI/3 - rtht[0]*Math.cos(rtht[1]), -rtht[0]*Math.sin(rtht[1])}; // left
		default:
			return null;
		}
	}
	
	
	private static double[] tetrahedralProjectionInverse(double x, double y) { // a function to help with tetrahedral projections
		if (y < x-1)
			return new double[] {
					-Math.PI/2, 0, 0,
					-Math.PI/2, Math.sqrt(3)*(x-2/3.), y+1 };
		else if (y < -x-1)
			return new double[] {
					-Math.PI/2, 0, 0,
					Math.PI/2, Math.sqrt(3)*(x+2/3.), y+1 };
		else if (y > -x+1)
			return new double[] {
					Math.PI/2-Math.asin(Math.sqrt(8)/3), Math.PI, 0,
					-Math.PI/2, Math.sqrt(3)*(x-2/3.), y-1 };
		else if (y > x+1)
			return new double[] {
					Math.PI/2-Math.asin(Math.sqrt(8)/3), Math.PI, 0,
					Math.PI/2, Math.sqrt(3)*(x+2/3.), y-1 };
		else if (x < 0)
			return new double[] {
					Math.PI/2-Math.asin(Math.sqrt(8)/3), -Math.PI/3, 0,
					Math.PI/6, Math.sqrt(3)*(x+1/3.), y };
		else
			return new double[] {
					Math.PI/2-Math.asin(Math.sqrt(8)/3), Math.PI/3, 0,
					-Math.PI/6, Math.sqrt(3)*(x-1/3.), y };
		
	}
	
	
	private static final double[] obliquifySphc(double latF, double lonF, double[] pole) { // go from polar coordinates to relative
		final double lat0 = pole[0];
		final double lon0 = pole[1];
		final double tht0 = pole[2];
		Vector r0 = new Vector (1, lat0, lon0);
		Vector rF = new Vector (1, latF, lonF);
		Vector r0XrF = r0.cross(rF);
		Vector r0Xk = r0.cross(Vector.K);
		
		double lat1 = Math.asin(r0.dot(rF)); // relative latitude
		double lon1;
		if (lat0 == Math.PI/2) // accounts for all the 0/0 errors at the poles
			lon1 = lonF-lon0;
		else if (lat0 == -Math.PI/2)
			lon1 = lon0-lonF+Math.PI;
		else {
			lon1 = Math.acos(r0XrF.dot(r0Xk)/(r0XrF.abs()*r0Xk.abs()))-Math.PI; // relative longitude
			if (Double.isNaN(lon1))
				lon1 = 0;
			else if (r0XrF.cross(r0Xk).dot(r0)/(r0XrF.abs()*r0Xk.abs()) > 0) // it's a plus-or-minus arccos.
				lon1 = 2*Math.PI-lon1;
		}
		lon1 = lon1-tht0;
		lon1 = Math2.mod(lon1+Math.PI, 2*Math.PI) - Math.PI;
		
		return new double[] {lat1, lon1};
	}
	
	
	private static final double[] obliquifyPlnr(double[] coords, double[] pole) { //go from relative coordinates to polar
		if (coords == null) 	return null;
		
		double lat1 = coords[0];
		double lon1 = coords[1];
		final double lat0 = pole[0];
		final double lon0 = pole[1];
		final double tht0 = pole[2];
		lon1 += tht0;
		double latf = Math.asin(Math.sin(lat0)*Math.sin(lat1) - Math.cos(lat0)*Math.cos(lon1)*Math.cos(lat1));
		double lonf;
		double innerFunc = Math.sin(lat1)/Math.cos(lat0)/Math.cos(latf) - Math.tan(lat0)*Math.tan(latf);
		if (lat0 == Math.PI / 2) // accounts for special case when lat0 = pi/2
			lonf = lon1+lon0;
		else if (lat0 == -Math.PI / 2) // accounts for special case when lat0 = -pi/2
			lonf = -lon1+lon0 + Math.PI;
		else if (Math.abs(innerFunc) > 1) { // accounts for special case when cos(lat1) = --> 0
			if ((lon1 == 0 && lat1 < -lat0) || (lon1 != 0 && lat1 < lat0))
				lonf = lon0 + Math.PI;
			else
				lonf = lon0;
		}
		else if (Math.sin(lon1) > 0)
			lonf = lon0 +
					Math.acos(innerFunc);
		else
			lonf = lon0 -
					Math.acos(innerFunc);
		
		double thtf = pole[2];
		
		double[] output = {latf, lonf, thtf};
		return output;
	}
	
	
	@Override
	public String toString() {
		return this.getName();
	}
	
	
	public String getName() {
		return this.name;
	}
	
	public String getDescription() {
		return this.description;
	}
	
	public boolean isParametrized() {
		return this.paramNames.length > 0;
	}
	
	public int getNumParameters() {
		return this.paramNames.length;
	}
	
	public String[] getParameterNames() {
		return this.paramNames;
	}
	
	public double[] getDefaultParameters() {
		final double[] params = new double[this.getNumParameters()];
		for (int i = 0; i < this.getNumParameters(); i ++)
			params[i] = this.paramValues[i][2];
		return params;
	}
	
	public double[][] getParameterValues() {
		return this.paramValues;
	}
	
	public double getAspectRatio(double[] params) {
		return this.aspectRatio;
	}
	
	public boolean isFinite() {
		return this.finite;
	}
	
	public boolean isInvertable() {
		return this.invertable;
	}
	
	public boolean isSolveable() {
		return this.solveable;
	}
	
	public boolean isContinuous() {
		return this.continuous;
	}
	
	public String getType() {
		return this.type;
	}
	
	public String getProperty() {
		return this.property;
	}

}
