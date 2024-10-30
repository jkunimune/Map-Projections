package maps;

import utils.Shape;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.hypot;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.toRadians;
import static utils.Math2.floorMod;

/**
 * equations for Hammer's retroazimuthal projection
 */
public class HammerRetroazimuthal {
	
	/**
	 * the pure Hammer retroazimuthal projection, allowing overlap between the front and back hemispheres
	 * @param lat the latitude of the point to be projected in radians
	 * @param lon the longitude of the point to be projected in radians
	 * @param lat0 the latitude of the reference point in radians
	 * @param lon0 the longitude of the reference point in radians
	 * @return the x and y coordinates in radians
	 */
	public static double[] project(double lat, double lon, double lat0, double lon0) {
		final double Z1 = sin(lat0)*sin(lat) + cos(lat0)*cos(lat)*cos(lon-lon0);
		final double K = acos(Z1)/sqrt(1 - pow(Z1, 2));
		final double x = K*cos(lat0)*sin(lon-lon0);
		final double y = -K*(sin(lat0)*cos(lat) - cos(lat0)*sin(lat)*cos(lon-lon0));
		return new double[] {x, y};
	}
	
	/**
	 * the Hammer retroazimuthal inverse projection, with a boolean argument to switch between the
	 * front and back hemispheres
	 * @param x the x value to which to project in radians
	 * @param y the y value to which to project in radians
	 * @param lat0 the latitude of the reference point in radians
	 * @param lon0 the longitude of the reference point in radians
	 * @param front whether to look for points with longitudes near the reference longitude
	 * @return the resulting latitude and longitude in radians
	 */
	public static double[] inverse(double x, double y, double lat0, double lon0, boolean front) {
		double lat1 = PI/2 - hypot(x, y);
		if (lat1 < -PI/2) 	return null;
		double lon1 = atan2(x, -y);
		double X1 = cos(lat1)*sin(lon1);
		double Y1 = -cos(lat1)*cos(lon1);
		double Z1 = sin(lat1);
		double sinΔlat = sin(lat0)/hypot(Y1, Z1);
		double sinΔlon = X1/cos(lat0);
		if (abs(sinΔlat) > 1 || abs(sinΔlon) > 1)
			return null;
		double lat, lon;
		if (front) {
			lat = asin(sinΔlat) + atan2(Y1, Z1);
			lon = lon0 + asin(sinΔlon);
		}
		else {
			lat = PI - asin(sinΔlat) + atan2(Y1, Z1);
			lat = floorMod(lat + PI, 2*PI) - PI;
			lon = lon0 + PI - asin(sinΔlon);
		}
		lon = floorMod(lon + PI, 2*PI) - PI;
		if (lat < -PI/2 || lat > PI/2)
			return null;
		else
			return new double[] {lat, lon};
	}
	
	
	public static final Projection FRONT = new Projection(
			"Hammer Retroazimuthal (front)", "Ernst H. H. Hammer",
			"A map where bearing and distance to a reference point is preserved " +
			"(showing the hemisphere nearest the reference point)",
			null, false, false, true, true, Projection.Type.PSEUDOCONIC, Projection.Property.RETROAZIMUTHAL, 2,
			new String[] {"Latitude","Longitude"},
			new double[][] {{-89,89,21.4}, {-180,180,39.8}}, false) {
		
		private double lat0, lon0;
		
		public void initialize(double... params) {
			this.lat0 = toRadians(params[0]);
			this.lon0 = toRadians(params[1]);
			this.shape = Shape.meridianEnvelope(
					this, -PI/2, PI/2, lon0 - PI/2, lon0 + PI/2);
		}
		
		public double[] project(double lat, double lon) {
			if (floorMod(lon - lon0 + PI/2, 2*PI) > 3*PI/2)
				return HammerRetroazimuthal.project(lat, lon0 - PI/2, lat0, lon0);
			else if ((abs(lon - lon0) + PI/2)%(2*PI) > PI)
				return HammerRetroazimuthal.project(lat, lon0 + PI/2, lat0, lon0);
			else
				return HammerRetroazimuthal.project(lat, lon, lat0, lon0);
		}
		
		public double[] inverse(double x, double y) {
			return HammerRetroazimuthal.inverse(x, y, lat0, lon0, true);
		}
	};
	
	
	public static final Projection BACK = new Projection(
			"Hammer Retroazimuthal (back)", "Ernst H. H. Hammer",
			"A map where bearing and distance to a reference point is preserved " +
			"(showing the hemisphere farther from the reference point)",
			null, false, false, true, true, Projection.Type.PSEUDOCONIC, Projection.Property.RETROAZIMUTHAL, 2,
			new String[] {"Latitude","Longitude"},
			new double[][] {{-89,89,21.4}, {-180,180,39.8}}, false) {
		
		private double lat0, lon0;
		
		public void initialize(double... params) {
			this.lat0 = toRadians(params[0]);
			this.lon0 = toRadians(params[1]);
			Shape outerShape = Shape.circle(PI);
			Shape innerShape = Shape.meridianEnvelope(
					this, -PI/2, PI/2, lon0 + PI/2, lon0 + 3*PI/2);
			this.shape = Shape.combination(outerShape, innerShape);
		}
		
		public double[] project(double lat, double lon) {
			if (floorMod(lon - lon0 + PI/2, 2*PI) < PI/2)
				return HammerRetroazimuthal.project(lat, lon0 - PI/2, lat0, lon0);
			else if ((abs(lon - lon0) + PI/2)%(2*PI) < PI)
				return HammerRetroazimuthal.project(lat, lon0 + PI/2, lat0, lon0);
			else
				return HammerRetroazimuthal.project(lat, lon, lat0, lon0);
		}
		
		public double[] inverse(double x, double y) {
			return HammerRetroazimuthal.inverse(x, y, lat0, lon0, false);
		}
	};
	
	
	public static final Projection FULL = new Projection(
			"Hammer Retroazimuthal", "Ernst H. H. Hammer",
			"The full version of a map where bearing and distance to a reference point is " +
			"preserved (points separated by more than 90° of longitude have been rotated " +
			"180° to fit)",
			null, false, true, true, true, Projection.Type.PSEUDOCONIC, Projection.Property.RETROAZIMUTHAL, 2,
			new String[] {"Latitude","Longitude"},
			new double[][] {{-89,89,21.4}, {-180,180,39.8}}, false) {
		
		private double lat0, lon0;
		
		public void initialize(double... params) {
			this.lat0 = toRadians(params[0]);
			this.lon0 = toRadians(params[1]);
			FRONT.initialize(params);
			BACK.initialize(params);
			this.shape = Shape.combination(
					FRONT.getShape(), Shape.scaled(-1, -1, BACK.getShape()));
		}
		
		public double[] project(double lat, double lon) {
			if (floorMod(lon - lon0 + PI/2, 2*PI) <= PI)
				return HammerRetroazimuthal.project(lat, lon, lat0, lon0);
			else {
				double[] xy = HammerRetroazimuthal.project(lat, lon, lat0, lon0);
				return new double[] {-xy[0], -xy[1]};
			}
		}
		
		public double[] inverse(double x, double y) {
			double r = hypot(x, y);
			if (r <= PI/2 - abs(lat0) || (r <= PI/2 + abs(lat0) && (y > 0) != (lat0 > 0)))
				return HammerRetroazimuthal.inverse(x, y, lat0, lon0, true);
			else
				return HammerRetroazimuthal.inverse(-x, -y, lat0, lon0, false);
		}
	};
	
}
