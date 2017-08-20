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
import utils.Dixon;

/**
 * Projections created by projecting onto and then unfolding a regular tetrahedron
 * 
 * @author jkunimune
 */
public class Tetrahedral {
	
	public static final Projection LEE =
			new TetrahedralProjection("Lee", Math.sqrt(3), 0b1001, Property.CONFORMAL,
					null, "that really deserves more attention") {
		
		public double[] innerProject(double lat, double lon) {
			final mfc.field.Complex z = mfc.field.Complex.fromPolar(
					Math.pow(2, 5/6.)*Math.tan(Math.PI/4-lat/2), lon);
			final mfc.field.Complex w = Dixon.invFunc(z);
			return new double[] { w.abs()*1.186, w.arg() };
		}
		
		public double[] innerInverse(double r, double tht) {
			final mfc.field.Complex w = mfc.field.Complex.fromPolar(
					r*1.53, tht - Math.PI/2);
			final mfc.field.Complex ans = Dixon.leeFunc(w).times(Math.pow(2, -5/6.));
			return new double[] {
					Math.PI/2 - 2*Math.atan(ans.abs()),
					ans.arg() + Math.PI };
		}
	};
	
	
	public static final Projection TETRAGRAPH =
			new TetrahedralProjection("TetraGraph", Math.sqrt(3), 0b1111, Property.EQUIDISTANT,
					null, "that I invented") {
		
		public double[] innerProject(double lat, double lon) {
			final double tht = lon - Math.floor(lon/(2*Math.PI/3))*(2*Math.PI/3) - Math.PI/3;
			return new double[] {
					Math.atan(1/Math.tan(lat)*Math.cos(tht))/Math.cos(tht)*Math.PI/3/Math.atan(Math.sqrt(2)),
					lon };
		}
		
		public double[] innerInverse(double r, double tht) {
			final double t0 = Math.floor((tht+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double dt = tht-t0;
			return new double[] {
					Math.PI/2 - Math.atan(Math.tan(r*Math.cos(dt)*Math.atan(Math.sqrt(2))/(Math.sqrt(3)/3))/Math.cos(dt)),
					Math.PI/2 + tht};
		}
	};
	
	
	public static final Projection TETRAPOWER =
			new TetrahedralProjection(
					"Tetrapower", "A parametric projection that I'm still testing",
					Math.sqrt(3), 0b1111, Property.COMPROMISE, new String[] {"k1","k2","k3"},
					new double[][] {{.25,4.,.87},{.25,4.,1.67},{.25,4.,1.00}}) {
		
		private double k1, k2, k3;
		
		public void setParameters(double... params) {
			this.k1 = params[0];
			this.k2 = params[1];
			this.k3 = params[2];
		}
		
		public double[] innerProject(double lat, double lon) {
			final double t0 = Math.floor(lon/(2*Math.PI/3))*(2*Math.PI/3) + Math.PI/3;
			final double tht = lon - t0;
			final double thtP = Math.PI/3*(1 - Math.pow(1-Math.abs(tht)/(Math.PI/2),k1))/(1 - 1/Math.pow(3,k1))*Math.signum(tht);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax = .5/Math.cos(thtP); //the max normalized radius of this triangle (in the plane)
			final double rtgf = Math.atan(1/Math.tan(lat)*Math.cos(tht))/Math.atan(Math.sqrt(2))*rmax;
			return new double[] {
					(1 - Math.pow(1-rtgf,kRad))/(1 - Math.pow(1-rmax,kRad))*rmax*2*Math.PI/3,
					thtP + t0 };
		}
		
		public double[] innerInverse(double r, double tht) {
			final double R = r*Math.sqrt(3)/2;
			final double t0 = Math.floor((tht+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double thtP = tht-t0;
			final double lamS = (1-Math.pow(1-Math.abs(thtP)*(1-1/Math.pow(3,k1))/(Math.PI/3), 1/k1))*Math.PI/2*Math.signum(thtP);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax = .5/Math.cos(thtP); //the max normalized radius of this triangle (in the plane)
			final double rtgf = 1-Math.pow(1-R/rmax*(1-Math.pow(Math.abs(1-rmax), kRad)), 1/kRad); //normalized tetragraph radius
			return new double[] {
					Math.atan(Math.cos(lamS)/Math.tan(rtgf/rmax*Math.atan(Math.sqrt(2)))),
					Math.PI/2 + t0 + lamS };
		}
	};
	
	
	public static final Projection TETRAFILLET =
			new TetrahedralProjection("TetraFillet", "A parametric projection that I'm still testing",
					Math.sqrt(3), 0b1110, Property.COMPROMISE, new String[] {"k1","k2","k3"},
					new double[][] {{.25,4.,1.1598},{.25,4.,.36295},{.25,4.,1.9553}}) {
		
		private double k1, k2, k3;
		
		public void setParameters(double... params) {
			this.k1 = params[0];
			this.k2 = params[1];
			this.k3 = params[2];
		}
		
		public double[] innerProject(double lat, double lon) {
			final double t0 = Math.floor(lon/(2*Math.PI/3))*(2*Math.PI/3) + Math.PI/3;
			final double tht = lon - t0;
			final double thtP = Math.PI/3*(1 - Math.pow(1-Math.abs(tht)/(Math.PI/2),k1))/(1 - 1/Math.pow(3,k1))*Math.signum(tht);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax; //the max normalized radius of this triangle (in the plane)
			if (Math.abs(thtP) < .70123892) 	rmax = .5/Math.cos(thtP);
			else 	rmax = .75 - 1.5972774*Math.pow(Math.PI/3-Math.abs(thtP),2)/2;
			final double rtgf = Math.atan(1/Math.tan(lat)*Math.cos(tht))/Math.atan(Math.sqrt(2))*rmax; //normalized tetragraph radius
			return new double[] {
					(1 - Math.pow(1-rtgf,kRad))/(1 - Math.pow(1-rmax,kRad))*rmax*2*Math.PI/3,
					thtP + t0
			};
		}
		
		public double[] innerInverse(double r, double tht) {
			final double R = r*Math.sqrt(3)/2;
			final double t0 = Math.floor((tht+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double thtP = tht-t0;
			final double lamS = (1-Math.pow(1-Math.abs(thtP)*(1-1/Math.pow(3,k1))/(Math.PI/3), 1/k1))*Math.PI/2*Math.signum(thtP);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax; //the max normalized radius of this triangle (in the plane)
			if (Math.abs(thtP) < .70123892) 	rmax = .5/Math.cos(thtP);
			else 	rmax = .75 - 1.5972774*Math.pow(Math.PI/3-Math.abs(thtP),2)/2;
			final double rtgf = 1-Math.pow(1-R/rmax*(1-Math.pow(Math.abs(1-rmax), kRad)), 1/kRad); //normalized tetragraph radius
			if (R > rmax) 	return null;
			return new double[] {
					Math.atan(Math.cos(lamS)/Math.tan(rtgf/rmax*Math.atan(Math.sqrt(2)))),
					Math.PI/2 + t0 + lamS };
		}
	};
	
	
	public static final Projection TETRACHAMFER =
			new TetrahedralProjection(
					"TetraChamfer", "A parametric projection that I'm still testing", Math.sqrt(3),
					0b1110, Property.COMPROMISE, new String[] {"k1","k2","k3"},
					new double[][] {{.25,4.,1.1598},{.25,4.,.36295},{.25,4.,1.9553}}) {
		
		private double k1, k2, k3;
		
		public void setParameters(double... params) {
			this.k1 = params[0];
			this.k2 = params[1];
			this.k3 = params[2];
		}
		
		public double[] innerProject(double lat, double lon) {
			final double t0 = Math.floor(lon/(2*Math.PI/3))*(2*Math.PI/3) + Math.PI/3;
			final double tht = lon - t0;
			final double thtP = Math.PI/3*(1 - Math.pow(1-Math.abs(tht)/(Math.PI/2),k1))/(1 - 1/Math.pow(3,k1))*Math.signum(tht);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax = Math.min(.5/Math.cos(thtP), .75/Math.cos(Math.PI/3-Math.abs(thtP))); //the max normalized radius of this triangle (in the plane)
			final double rtgf = Math.atan(1/Math.tan(lat)*Math.cos(tht))/Math.atan(Math.sqrt(2))*rmax; //normalized tetragraph radius
			return new double[] {
					(1 - Math.pow(1-rtgf,kRad))/(1 - Math.pow(1-rmax,kRad))*rmax*2*Math.PI/3,
					thtP + t0 };
		}
		
		public double[] innerInverse(double r, double tht) {
			final double R = r*Math.sqrt(3)/2;
			final double t0 = Math.floor((tht+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double thtP = tht-t0;
			final double lamS = (1-Math.pow(1-Math.abs(thtP)*(1-1/Math.pow(3,k1))/(Math.PI/3), 1/k1))*Math.PI/2*Math.signum(thtP);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax = Math.min(.5/Math.cos(thtP), .75/Math.cos(Math.PI/3-Math.abs(thtP))); //the max normalized radius of this triangle (in the plane)
			final double rtgf = 1-Math.pow(1-R/rmax*(1-Math.pow(Math.abs(1-rmax), kRad)), 1/kRad); //normalized tetragraph radius
			if (R > rmax) 	return null;
			return new double[] {
					Math.atan(Math.cos(lamS)/Math.tan(rtgf/rmax*Math.atan(Math.sqrt(2)))),
					Math.PI/2 + t0 + lamS };
		}
	};
	
	
	public static final Projection AUTHAGRAPH =
			new TetrahedralProjection("AuthaGraph",
					"A hip new Japanese map that is almost authagraphic (this is an approximation; they won't give me their actual equations)",
					4/Math.sqrt(3), 0b1001, Property.COMPROMISE) {
		
		public double[] project(double lat, double lon) {
			return null; //TODO Can I live with myself if I never implement this?
		}
		
		public double[] innerProject(double lat, double lon) {
			return null;
		}
		
		public double[] inverse(double x, double y) {
			final double[] faceCenter;
			final double dt, xp, yp;
			if (y-1 < 4*x && y-1 < -4*x) {
				faceCenter = new double[] { Math.PI/2-Math.asin(Math.sqrt(8)/3), 0, 0};
				dt = 0;
				xp = 4/Math.sqrt(3)*x;
				yp = y+1/3.0;
			}
			else if (y-1 < -4*(x+1)) {
				faceCenter = new double[] { -Math.PI/2, Math.PI, 0 };
				dt = 0;
				xp = 4/Math.sqrt(3)*(x+1);
				yp = y+1/3.0;
			}
			else if (y-1 < 4*(x-1)) {
				faceCenter = new double[] { -Math.PI/2, Math.PI, 0 };
				dt = 0;
				xp = 4/Math.sqrt(3)*(x-1);
				yp = y+1/3.0;
			}
			else if (x < 0) {
				faceCenter = new double[] { Math.PI/2-Math.asin(Math.sqrt(8)/3), 4*Math.PI/3, 0 };
				dt = Math.PI/3;
				xp = 4/Math.sqrt(3)*(x+0.5);
				yp = y-1/3.0;
			}
			else {
				faceCenter = new double[] { Math.PI/2-Math.asin(Math.sqrt(8)/3), 2*Math.PI/3, 0 };
				dt = -Math.PI/3;
				xp = 4/Math.sqrt(3)*(x-0.5);
				yp = y-1/3.0;
			}
			
			return obliquifyPlnr(
					innerInverse(Math.hypot(xp, yp), Math.atan2(yp, xp)+dt), faceCenter);
		}
		
		protected double[] innerInverse(double r, double tht) {
			final double t0 = Math.floor((tht+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
			final double dt = tht-t0;
			final double z = 2.49*r*Math.cos(dt);
			final double g = 0.03575*z*z*z + 0.0219*z*z + 0.4441*z;
			return new double[] {
					Math.PI/2 - Math.atan(Math.tan(g)/Math.cos(dt)),
					Math.PI/2 + t0 + dt };
		}
	};
	
	
	
	/**
	 * A base for tetrahedral projections
	 * 
	 * @author jkunimune
	 */
	private static abstract class TetrahedralProjection extends Projection {
		
		public TetrahedralProjection(
				String name, double aspectRatio, int fisc, Property property, String adjective, String addendum) {
			super(name, aspectRatio, fisc, Type.TETRAHEDRAL, property, adjective, addendum);
		}
		
		public TetrahedralProjection(
				String name, String description, double aspectRatio, int fisc, Property property) {
			super(name, description, aspectRatio, fisc, Type.TETRAHEDRAL, property);
		}
		
		public TetrahedralProjection(
				String name, String description, double aspectRatio, int fisc, Property property,
				String[] paramNames, double[][] paramValues) {
			super(name, description, aspectRatio, fisc, Type.TETRAHEDRAL, property, paramNames,
					paramValues);
		}
		
		
		protected abstract double[] innerProject(double lon, double lat);
		
		protected abstract double[] innerInverse(double x, double y);
		
		
		public double[] project(double lat, double lon) {
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
			
			final double[] rtht = innerProject(latR, lonR);
			
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
		
		
		public double[] inverse(double x, double y) {
			final double[] faceCenter;
			final double dt, xp, yp;
			if (y < x-1) {
				faceCenter = new double[] { -Math.PI/2, 0, 0 };
				dt = -Math.PI/2;
				xp = Math.sqrt(3)*(x-2/3.);
				yp = y+1;
			}
			else if (y < -x-1) {
				faceCenter = new double[] { -Math.PI/2, 0, 0 };
				dt = Math.PI/2;
				xp = Math.sqrt(3)*(x+2/3.);
				yp = y+1;
			}
			else if (y > -x+1) {
				faceCenter = new double[] { Math.PI/2-Math.asin(Math.sqrt(8)/3), Math.PI, 0 };
				dt = -Math.PI/2;
				xp = Math.sqrt(3)*(x-2/3.);
				yp = y-1;
			}
			else if (y > x+1) {
				faceCenter = new double[] { Math.PI/2-Math.asin(Math.sqrt(8)/3), Math.PI, 0 };
				dt = Math.PI/2; xp = Math.sqrt(3)*(x+2/3.); yp = y-1;
			}
			else if (x < 0) {
				faceCenter = new double[] { Math.PI/2-Math.asin(Math.sqrt(8)/3), -Math.PI/3, 0 };
				dt = Math.PI/6;
				xp = Math.sqrt(3)*(x+1/3.);
				yp = y;
			}
			else {
				faceCenter = new double[] { Math.PI/2-Math.asin(Math.sqrt(8)/3), Math.PI/3, 0 };
				dt = -Math.PI/6;
				xp = Math.sqrt(3)*(x-1/3.);
				yp = y;
			}
			
			return obliquifyPlnr(
					innerInverse(Math.hypot(xp, yp), Math.atan2(yp, xp)+dt), faceCenter);
		}
	}
}
