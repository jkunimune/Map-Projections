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
import utils.NumericalAnalysis;

/**
 * Projections created by projecting onto and then unfolding a regular tetrahedron
 * 
 * @author jkunimune
 */
public class Tetrahedral {
	
	private static final double ASIN_ONE_THD = Math.asin(1/3.); //the complement of the angular radius of a face
	
	public static final Projection LEE =
			new TetrahedralProjection(
					"Lee", 0b1001, Configuration.WIDE_FACE, Property.CONFORMAL,
					null, "that really deserves more attention") {
		
		public double[] innerProject(double lat, double lon) {
			final de.jtem.mfc.field.Complex z = de.jtem.mfc.field.Complex.fromPolar(
					Math.pow(2, 5/6.)*Math.tan(Math.PI/4-lat/2), lon);
			final de.jtem.mfc.field.Complex w = Dixon.invFunc(z);
			return new double[] { w.abs()*1.132, w.arg() }; //I don't understand Dixon functions well enough to say whence the 1.132 comes
		}
		
		public double[] innerInverse(double r, double tht) {
			final de.jtem.mfc.field.Complex w = de.jtem.mfc.field.Complex.fromPolar(
					r/1.132, tht);
			final de.jtem.mfc.field.Complex ans = Dixon.leeFunc(w).times(Math.pow(2, -5/6.));
			return new double[] {
					Math.PI/2 - 2*Math.atan(ans.abs()),
					ans.arg() };
		}
	};
	
	
	public static final Projection TETRAGRAPH =
			new TetrahedralProjection(
					"TetraGraph", 0b1111, Configuration.WIDE_FACE, Property.EQUIDISTANT,
					null, "that I invented") {
		
		public double[] innerProject(double lat, double lon) {
			final double tht = lon - Math.floor((lon+Math.PI/3)/(2*Math.PI/3))*(2*Math.PI/3);
			return new double[] {
					Math.atan(1/Math.tan(lat)*Math.cos(tht))/Math.cos(tht) /Math.atan(Math.sqrt(2)),
					lon };
		}
		
		public double[] innerInverse(double r, double tht) {
			final double t0 = Math.floor((tht+Math.PI/3)/(2*Math.PI/3))*(2*Math.PI/3);
			final double dt = tht-t0;
			return new double[] {
					Math.PI/2 - Math.atan(Math.tan(r*Math.cos(dt)*Math.atan(Math.sqrt(2)))/Math.cos(dt)),
					tht };
		}
	};
	
	
	public static final Projection AUTHAGRAPH = 
			new TetrahedralProjection(
					"AuthaGraph", "A hip new Japanese map that is not authagraphic.",
					0b1011, Configuration.AUTHAGRAPH, Property.COMPROMISE) {
		
		public double[] innerProject(double lat, double lon) {
			final double lam = lon - (Math.floor(lon/(2*Math.PI/3))*(2*Math.PI/3)+Math.PI/3);
			final double tht = Math.atan((lam - Math.asin(Math.sin(lam)/Math.sqrt(3)))/Math.PI*Math.sqrt(12));
			final double p = (Math.PI/2 - lat) / Math.atan(Math.sqrt(2)/Math.cos(lam));
			return new double[] { Math.pow(p,.707)*Math.sqrt(3)/Math.cos(tht), (lon - lam)/2 + tht };
		}
		
		protected double[] innerInverse(double r, double th) {
			final double dt = th - (Math.floor(th/(Math.PI/3))*(Math.PI/3)+Math.PI/6);
			final double dl = NumericalAnalysis.newtonRaphsonApproximation(dt, dt*2,
					(l) -> Math.atan((l - Math.asin(Math.sin(l)/Math.sqrt(3)))/Math.PI*Math.sqrt(12)),
					(l) -> (1-1/Math.sqrt(1+2*Math.pow(Math.cos(l),-2)))/Math.sqrt(Math.pow(Math.PI,2)/12+Math.pow(l-Math.asin(Math.sin(l)/Math.sqrt(3)),2)),
					.01);
			final double R = r / (Math.sqrt(3)/Math.cos(dt));
			return new double[] {
					Math.PI/2 - Math.pow(R,1.41)*Math.atan(Math.sqrt(2)/Math.cos(dl)),
					(th - dt)*2 + dl };
		}
	};
	
	
	public static final Projection AUTHAPOWER = //yes, these names are getting a bit redonkadonk. I have a lot of map projections, okay? I'll prune this when I finish my optimiser
			new TetrahedralProjection(
					"AuthaPower", "A parametrised, rearranged version of my AuthaGraph approximation",
					0b1011, Configuration.WIDE_VERTEX, Property.COMPROMISE,
					new String[] {"Power"}, new double[][] {{0,1,.7}}) {
		
		private double k;
		
		public void setParameters(double... params) {
				this.k = params[0];
		}
		
		public double[] innerProject(double lat, double lon) {
			final double lam = lon - (Math.floor(lon/(2*Math.PI/3))*(2*Math.PI/3)+Math.PI/3);
			final double tht = Math.atan((lam - Math.asin(Math.sin(lam)/Math.sqrt(3)))/Math.PI*Math.sqrt(12));
			final double p = (Math.PI/2 - lat) / Math.atan(Math.sqrt(2)/Math.cos(lam));
			return new double[] { Math.pow(p,k)*Math.sqrt(3)/Math.cos(tht), (lon - lam)/2 + tht };
		}
		
		protected double[] innerInverse(double r, double th) {
			final double dt = th - (Math.floor(th/(Math.PI/3))*(Math.PI/3)+Math.PI/6);
			final double dl = NumericalAnalysis.newtonRaphsonApproximation(dt, dt*2,
					(l) -> Math.atan((l - Math.asin(Math.sin(l)/Math.sqrt(3)))/Math.PI*Math.sqrt(12)),
					(l) -> (1-1/Math.sqrt(1+2*Math.pow(Math.cos(l),-2)))/Math.sqrt(Math.pow(Math.PI,2)/12+Math.pow(l-Math.asin(Math.sin(l)/Math.sqrt(3)),2)),
					.01);
			final double R = r / (Math.sqrt(3)/Math.cos(dt));
			return new double[] {
					Math.PI/2 - Math.pow(R,1/k)*Math.atan(Math.sqrt(2)/Math.cos(dl)),
					(th - dt)*2 + dl };
		}
	};
	
	
	public static final Projection ACTUAUTHAGRAPH =
			new TetrahedralProjection(
					"Equahedral", 0b1011, Configuration.WIDE_VERTEX, Property.EQUAL_AREA,
					"to put AuthaGraph to shame", null) {
		
		public double[] innerProject(double lat, double lon) {
			final double lam = lon - (Math.floor(lon/(2*Math.PI/3))*(2*Math.PI/3)+Math.PI/3);
			final double tht = Math.atan((lam - Math.asin(Math.sin(lam)/Math.sqrt(3)))/Math.PI*Math.sqrt(12));
			return new double[] {
					Math.sqrt((1 - Math.sin(lat))/(1 - 1/Math.sqrt(1+2*Math.pow(Math.cos(lam),-2))))*Math.sqrt(3)/Math.cos(tht),
					(lon - lam)/2 + tht };
		}
		
		public double[] innerInverse(double r, double th) {
			final double dt = th - (Math.floor(th/(Math.PI/3))*(Math.PI/3)+Math.PI/6);
			final double dl = NumericalAnalysis.newtonRaphsonApproximation(dt, dt*2,
					(l) -> Math.atan((l - Math.asin(Math.sin(l)/Math.sqrt(3)))/Math.PI*Math.sqrt(12)),
					(l) -> (1-1/Math.sqrt(1+2*Math.pow(Math.cos(l),-2)))/Math.sqrt(Math.pow(Math.PI,2)/12+Math.pow(l-Math.asin(Math.sin(l)/Math.sqrt(3)),2)),
					.01);
			final double R = r / (Math.sqrt(3)/Math.cos(dt));
			final double p = Math.acos(1 - R*R*(1-1/Math.sqrt(1+2*Math.pow(Math.cos(dl), -2))));
			return new double[] { Math.PI/2 - p, (th - dt)*2 + dl };
		}
	};
	
	
	public static final Projection TETRAPOWER =
			new TetrahedralProjection(
					"TetraPower", "A parameterised tetrahedral projection that I invented.", 0b1111,
					Configuration.WIDE_FACE, Property.COMPROMISE, new String[] {"k1","k2","k3"},
					new double[][] {{.01,2.,.98},{.01,2.,1.2},{.01,2.,.98}}) {
		
		private double k1, k2, k3;
		
		public void setParameters(double... params) {
			this.k1 = params[0];
			this.k2 = params[1];
			this.k3 = params[2];
		}
		
		public double[] innerProject(double lat, double lon) {
			final double t0 = Math.floor((lon+Math.PI/3)/(2*Math.PI/3))*(2*Math.PI/3);
			final double tht = lon - t0;
			final double thtP = Math.PI/3*(1 - Math.pow(1-Math.abs(tht)/(Math.PI/2),k1))/(1 - 1/Math.pow(3,k1))*Math.signum(tht);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax = 0.5/Math.cos(thtP); //the max normalized radius of this triangle (in the plane)
			final double rtgf = Math.atan(1/Math.tan(lat)*Math.cos(tht))/Math.atan(Math.sqrt(2))*rmax;
			return new double[] {
					(1 - Math.pow(1-rtgf,kRad))/(1 - Math.pow(1-rmax,kRad))*rmax*2,
					thtP + t0 };
		}
		
		public double[] innerInverse(double r, double tht) {
			final double t0 = Math.floor((tht+Math.PI/3)/(2*Math.PI/3))*(2*Math.PI/3);
			final double thtP = tht-t0;
			final double lamS = (1-Math.pow(1-Math.abs(thtP)*(1-1/Math.pow(3,k1))/(Math.PI/3), 1/k1))*Math.PI/2*Math.signum(thtP);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax = 0.5/Math.cos(thtP); //the max normalized radius of this triangle (in the plane)
			final double rtgf = 1-Math.pow(1-r/2/rmax*(1-Math.pow(Math.abs(1-rmax), kRad)), 1/kRad); //normalized tetragraph radius
			return new double[] {
					Math.atan(Math.cos(lamS)/Math.tan(rtgf/rmax*Math.atan(Math.sqrt(2)))),
					t0 + lamS };
		}
	};
	
	
	public static final Projection TETRAFILLET =
			new TetrahedralProjection("TetraFillet",
					"A parameterised tetrahedral projection I invented with the corners filleted off.",
					0b1110, Configuration.WIDE_FACE, Property.COMPROMISE,
					new String[] {"k1","k2","k3"},
					new double[][] {{.01,2.,.78},{.01,2.,.99},{.01,2.,1.3}}) {
		
		private double k1, k2, k3;
		
		public void setParameters(double... params) {
			this.k1 = params[0];
			this.k2 = params[1];
			this.k3 = params[2];
		}
		
		public double[] innerProject(double lat, double lon) {
			final double t0 = Math.floor((lon+Math.PI/3)/(2*Math.PI/3))*(2*Math.PI/3);
			final double tht = lon - t0;
			final double thtP = Math.PI/3*(1 - Math.pow(1-Math.abs(tht)/(Math.PI/2),k1))/(1 - 1/Math.pow(3,k1))*Math.signum(tht);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax; //the max normalized radius of this triangle (in the plane)
			if (Math.abs(thtP) < .70123892) 	rmax = .5/Math.cos(thtP);
			else 	rmax = .75 - 1.5972774*Math.pow(Math.PI/3-Math.abs(thtP),2)/2;
			final double rtgf = Math.atan(1/Math.tan(lat)*Math.cos(tht))/Math.atan(Math.sqrt(2))*rmax; //normalized tetragraph radius
			return new double[] {
					(1 - Math.pow(1-rtgf,kRad))/(1 - Math.pow(1-rmax,kRad))*rmax*2,
					thtP + t0
			};
		}
		
		public double[] innerInverse(double r, double tht) {
			final double t0 = Math.floor((tht+Math.PI/3)/(2*Math.PI/3))*(2*Math.PI/3);
			final double thtP = tht-t0;
			final double lamS = (1-Math.pow(1-Math.abs(thtP)*(1-1/Math.pow(3,k1))/(Math.PI/3), 1/k1))*Math.PI/2*Math.signum(thtP);
			final double kRad = k3*Math.abs(thtP)/(Math.PI/3) + k2*(1-Math.abs(thtP)/(Math.PI/3));
			final double rmax; //the max normalized radius of this triangle (in the plane)
			if (Math.abs(thtP) < .70123892) 	rmax = .5/Math.cos(thtP);
			else 	rmax = .75 - 1.5972774*Math.pow(Math.PI/3-Math.abs(thtP),2)/2;
			final double rtgf = 1-Math.pow(1-r/2/rmax*(1-Math.pow(Math.abs(1-rmax), kRad)), 1/kRad); //normalized tetragraph radius
			if (r/2 > rmax) 	return null;
			return new double[] {
					Math.atan(Math.cos(lamS)/Math.tan(rtgf/rmax*Math.atan(Math.sqrt(2)))),
					t0 + lamS };
		}
	};
	
	
	
	/**
	 * A base for tetrahedral Projections
	 * 
	 * @author jkunimune
	 */
	private static abstract class TetrahedralProjection extends Projection {
		
		private final Configuration configuration;
		
		
		public TetrahedralProjection(
				String name, int fisc, Configuration config, Property property,
				String adjective, String addendum) {
			super(name, config.width, config.height, fisc, Type.TETRAHEDRAL, property,
					adjective, addendum);
			this.configuration = config;
		}
		
		public TetrahedralProjection(
				String name, String description, int fisc, Configuration config, Property property) {
			super(name, description, config.width, config.height, fisc, Type.TETRAHEDRAL, property);
			this.configuration = config;
		}
		
		public TetrahedralProjection(
				String name, String description, int fisc, Configuration config, Property property,
				String[] paramNames, double[][] paramValues) {
			super(name, description, config.width, config.height, fisc, Type.TETRAHEDRAL, property,
					paramNames, paramValues);
			this.configuration = config;
		}
		
		
		protected abstract double[] innerProject(double lat, double lon); //the projection from spherical to polar within a face
		
		protected abstract double[] innerInverse(double x, double y); //I think you can guess
		
		
		public double[] project(double lat, double lon) {
			double latR = Double.NEGATIVE_INFINITY;
			double lonR = Double.NEGATIVE_INFINITY;
			double[] centrum = null;
			for (double[] testCentrum: configuration.centrumSet) {
				final double[] relCoords = obliquifySphc(lat, lon, testCentrum);
				if (relCoords[0] > latR || //pick the centrum that maxes out your latitude
						(relCoords[0] == latR && Math.cos(relCoords[1]) > Math.cos(lonR))) { //or, in the event of a tie, the cosine of your longitude
					latR = relCoords[0];
					lonR = relCoords[1];
					centrum = testCentrum;
				}
			}
			
			final double[] rth = innerProject(latR, lonR); //apply the projection to the relative coordinates
			final double r = rth[0];
			final double th = rth[1] + centrum[3];
			final double x0 = centrum[4];
			final double y0 = centrum[5];
			
			double[] output = { r*Math.cos(th) + x0, r*Math.sin(th) + y0 };
			if (Math.abs(output[0]) > width/2 || Math.abs(output[1]) > height/2) { //rotate OOB bits around nearest singularity
				output = configuration.rotateOOB(output);
			}
			return output;
		}
		
		
		public double[] inverse(double x, double y) {
			double rM = Double.POSITIVE_INFINITY;
			double[] centrum = null;
			for (double[] testCentrum: configuration.centrumSet) {
				final double rR = Math.hypot(x-testCentrum[4], y-testCentrum[5]);
				if (rR < rM) {
					rM = rR;
					centrum = testCentrum;
				}
			}
			
			final double th0 = centrum[3];
			final double x0 = centrum[4];
			final double y0 = centrum[5];
			
			return obliquifyPlnr(
					innerInverse(Math.hypot(x-x0, y-y0), Math.atan2(y-y0, x-x0) - th0), centrum);
		}
	}
	
	
	
	/**
	 * A set of objects that determine the layouts of tetrahedral projections
	 * 
	 * @author jkunimune
	 */
	private static enum Configuration {
		
		/*		  LATITUDE,		 LONGITUDE,		 STD_PRLL,		 PLANE_ROT,		 X,	 Y */
		WIDE_FACE(6, 2*Math.sqrt(3), new double[][] { // [<|>] arrangement, face-centered
				{ ASIN_ONE_THD,	 Math.PI,		-2*Math.PI/3,	-2*Math.PI/3,	 2,	 Math.sqrt(3) },
				{ ASIN_ONE_THD,	 Math.PI,		 2*Math.PI/3,	-Math.PI/3,		-2,	 Math.sqrt(3) },
				{-Math.PI/2,	 0,				 2*Math.PI/3,	 2*Math.PI/3,	 2,	-Math.sqrt(3) },
				{-Math.PI/2,	 0,				-2*Math.PI/3,	 Math.PI/3,		-2,	-Math.sqrt(3) },
				{ ASIN_ONE_THD,	 Math.PI/3,		-2*Math.PI/3,	 Math.PI,		 1,	 0 },
				{ ASIN_ONE_THD,	-Math.PI/3,		 2*Math.PI/3,	 0,				-1,	 0 }}) {
			public double[] rotateOOB(double[] oob) {
				return new double[] { -oob[0], Math.sqrt(3)*Math.signum(oob[1]) - oob[1] };
			}
		},
		WIDE_VERTEX(6, 2*Math.sqrt(3), new double[][] { // [<|>] arrangement, vertex-centered
				{ Math.PI/2,	 0,				 Math.PI/3,		-Math.PI/3,		 0,	 Math.sqrt(3) },
				{-ASIN_ONE_THD,	 0,				 2*Math.PI/3,	 Math.PI/3,		 0,	-Math.sqrt(3) },
				{-ASIN_ONE_THD,	 2*Math.PI/3,	-2*Math.PI/3,	 Math.PI,		 3, 0 },
				{-ASIN_ONE_THD,	-2*Math.PI/3,	 2*Math.PI/3,	 0,				-3, 0 }}) {
			public double[] rotateOOB(double[] oob) {
				return new double[] { -oob[0], 2*Math.sqrt(3)*Math.signum(oob[1]) - oob[1] };
			}
		},
		AUTHAGRAPH(4*Math.sqrt(3), 3, new double[][] { // |\/\/`| arrangement, vertex-centered
				{-ASIN_ONE_THD, Math.PI,	 0,			 -Math.PI/2,-2*Math.sqrt(3)-.6096,  1.5 },
				{-ASIN_ONE_THD,-Math.PI/3,-2*Math.PI/3, Math.PI/2,-Math.sqrt(3)-.6096,   -1.5 },
				{ Math.PI/2,	0,			 0,	 -Math.PI/2, 0-.6096,			    1.5 },
				{-ASIN_ONE_THD, Math.PI/3, 2*Math.PI/3, Math.PI/2, Math.sqrt(3)-.6096,   -1.5 },
				{-ASIN_ONE_THD, Math.PI,	 0,			 -Math.PI/2, 2*Math.sqrt(3)-.6096,  1.5 },
				{-ASIN_ONE_THD,-Math.PI/3,-2*Math.PI/3, Math.PI/2, 3*Math.sqrt(3)-.6096, -1.5}}) {
			public double[] rotateOOB(double[] oob) {
				return new double[] { (oob[0]+3*width/2)%width-width/2, oob[1] };
			}
		};
		/*		  LATITUDE,	    LONGITUDE,	 STD_PRLL, PLANE_ROT,		 X,	 Y */
		
		public final double width, height; //the width and height of a map with this configuration
		public final double[][] centrumSet; //the mathematical information about this configuration
		public abstract double[] rotateOOB(double[] oob);
		
		private Configuration(double width, double height, double[][] centrumSet) {
			this.width = width;
			this.height = height;
			this.centrumSet = centrumSet;
		}
	}
}
