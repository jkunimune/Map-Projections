package utils;

import image.Path;
import maps.Projection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;


public class Shape {
	public final List<Path.Command> path;
	public final double xMin;
	public final double xMax;
	public final double width;
	public final double yMin;
	public final double yMax;
	public final double height;
	public final double aspectRatio;
	
	/**
	 * instantiate a new Shape with its dimensions and underlying Path.  the x and y limits must
	 * match the Path, because I don't have a good way to verify that.  in general you should prefer
	 * to use the static methods circle(), ellipse(), rectangle(), annularSector(), polygon(), and
	 * meridianEnvelope(), because those will all garantee consistency for you.
	 * @param xMin the x-coordinate of the leftmost point on the path
	 * @param xMax the x-coordinate of the rightmost point on the path
	 * @param yMin the y-coordinate of the bottom of the shape
	 * @param yMax the y-coordinate of the top of the shape
	 * @param path the underlying Path object (an List of Commands)
	 */
	public Shape(double xMin, double xMax, double yMin, double yMax, List<Path.Command> path) {
		this.path = path;
		this.xMin = xMin;
		this.xMax = xMax;
		this.width = xMax - xMin;
		this.yMin = yMin;
		this.yMax = yMax;
		this.height = yMax - yMin;
		this.aspectRatio = width/height;
	}


	/**
	 * construct a circle of some radius, centered on the origin
	 */
	public static Shape circle(double radius) {
		return ellipse(radius, radius);
	}


	/**
	 * construct an ellipse centered on the origin and aligned with the principle axes
	 */
	public static Shape ellipse(double semiaxisX, double semiaxisY) {
		return new Shape(
				-semiaxisX, semiaxisX,
				-semiaxisY, semiaxisY,
				Arrays.asList(
						new Path.Command('M', 0, semiaxisY),
						new Path.Command('A', semiaxisX, semiaxisY, 0, 0, 1, 0, -semiaxisY),
						new Path.Command('A', semiaxisX, semiaxisY, 0, 0, 1, 0, semiaxisY)));
	}


	/**
	 * construct a rectangle of a given width and height, centered on the origin
	 */
	public static Shape rectangle(double width, double height) {
		return new Shape(
				-width/2, width/2, -height/2, height/2,
				Arrays.asList(
						new Path.Command('M', -width/2, -height/2),
						new Path.Command('h', width),
						new Path.Command('v', height),
						new Path.Command('h', -width),
						new Path.Command('Z')));
	}


	/**
	 * construct a thick arc shape symmetric about the x axis, composed of two lines coincident with the origin and
	 * two arcs centered on the origin.
	 * @param innerRadius the radius of the smaller arc segment
	 * @param outerRadius the radius of the larger arc segment
	 * @param angularWidth the angle between the two line segments
	 * @param facesDown true if the shape is above the origin and false if it is below the origin
	 */
	public static Shape annularSector(double innerRadius, double outerRadius, double angularWidth, boolean facesDown) {
		if (facesDown) {
			double angle = angularWidth/2;
			boolean obtuse = angularWidth >= PI;
			List<Path.Command> path = Arrays.asList(
					new Path.Command('M', outerRadius*sin(angle), outerRadius*cos(angle)),
					new Path.Command('A', outerRadius, outerRadius, 0, obtuse ? 1 : 0, 1, -outerRadius*sin(angle), outerRadius*cos(angle)),
					new Path.Command('L', -innerRadius*sin(angle), innerRadius*cos(angle)),
					new Path.Command('A', innerRadius, innerRadius, 0, obtuse ? 1 : 0, 0, innerRadius*sin(angle), innerRadius*cos(angle)),
					new Path.Command('Z'));

			// for obtuse annular sectors
			if (obtuse)
				return new Shape(
						-outerRadius, outerRadius,
						outerRadius*cos(angularWidth/2), outerRadius,
						path);
			else
				return new Shape(
						-outerRadius*sin(angularWidth/2), outerRadius*sin(angularWidth/2),
						innerRadius*cos(angularWidth/2), outerRadius,
						path);
		}
		// for upside-down annular sectors, generate a rightside-up one and invert it
		else {
			Shape reverse = annularSector(innerRadius, outerRadius, angularWidth, true);
			return scaled(-1, -1, reverse);
		}
	}


	/**
	 * construct a polygon from a set of vertices.
	 * @param points each row is a pair of x and y coordinates. the last row does not need to repeat the first row.
	 */
	public static Shape polygon(double[][] points) {
		double xMin = POSITIVE_INFINITY, xMax = NEGATIVE_INFINITY;
		double yMin = POSITIVE_INFINITY, yMax = NEGATIVE_INFINITY;
		List<Path.Command> commands = new ArrayList<>(points.length + 1);
		for (int i = 0; i < points.length; i ++) {
			commands.add(new Path.Command((i == 0) ? 'M' : 'L', points[i]));
			if (points[i][0] < xMin)
				xMin = points[i][0];
			if (points[i][0] > xMax)
				xMax = points[i][0];
			if (points[i][1] < yMin)
				yMin = points[i][1];
			if (points[i][1] > yMax)
				yMax = points[i][1];
		}
		commands.add(new Path.Command('Z'));
		return new Shape(xMin, xMax, yMin, yMax, commands);
	}


	/**
	 * infer the outline of a map projection by tracing the pole lines and prime meridians
	 */
	public static Shape meridianEnvelope(Projection projection) {
		return meridianEnvelope(projection, -PI/2, PI/2, -PI, PI);
	}


	/**
	 * infer the outline of a map projection by tracing certain parallels and meridians
	 */
	public static Shape meridianEnvelope(
			Projection projection, double фMin, double фMax, double λMin, double λMax) {
		List<Path.Command> east = projection.drawLoxodrome(фMax, λMin, фMin, λMin, .0005); // east meridian
		List<Path.Command> south = projection.drawLoxodrome(фMin, λMin, фMin, λMax, .0005); // south parallel
		List<Path.Command> west = projection.drawLoxodrome(фMin, λMax, фMax, λMax, .0005); // west meridian
		List<Path.Command> north = projection.drawLoxodrome(фMax, λMax, фMax, λMin, .0005); // north parallel
		List<Path.Command> envelope = new ArrayList<>(east.size() + south.size() + west.size() + north.size());
		envelope.addAll(east);
		envelope.addAll(south.subList(1, south.size()));
		envelope.addAll(west.subList(1, west.size()));
		envelope.addAll(north.subList(1, north.size()));
		return polygon(Path.asArray(envelope));
	}


	/**
	 * take an existing shape and scale it along the x- and y- axes
	 */
	public static Shape scaled(double xScale, double yScale, Shape shape) {
		List<Path.Command> path = Path.scaled(xScale, yScale, shape.path);
		double xMin = shape.xMin*xScale;
		double xMax = shape.xMax*xScale;
		if (xScale < 0) {
			double newXMin = xMax;
			xMax = xMin;
			xMin = newXMin;
		}
		double yMin = shape.yMin*yScale;
		double yMax = shape.yMax*yScale;
		if (yScale < 0) {
			double newYMin = yMax;
			yMax = yMin;
			yMin = newYMin;
		}
		return new Shape(xMin, xMax, yMin, yMax, path);
	}
	
	
	/**
	 * combine two shapes.  whether this is a union or subtraction depends on the polarity of the components' paths;
	 * all this really does is overlay their paths and union their bounding boxen.
	 */
	public static Shape combination(Shape... components) {
		double xMin = POSITIVE_INFINITY, xMax = NEGATIVE_INFINITY;
		double yMin = POSITIVE_INFINITY, yMax = NEGATIVE_INFINITY;
		List<Path.Command> path = new ArrayList<>();
		for (Shape component: components) {
			if (component.xMin < xMin)
				xMin = component.xMin;
			if (component.xMax > xMax)
				xMax = component.xMax;
			if (component.yMin < yMin)
				yMin = component.yMin;
			if (component.yMax > yMax)
				yMax = component.yMax;
			path.addAll(component.path);
		}
		return new Shape(xMin, xMax, yMin, yMax, path);
	}
}
