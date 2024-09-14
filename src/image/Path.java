package image;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static java.lang.Double.isFinite;
import static java.lang.Double.parseDouble;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.sin;
import static java.lang.Math.toDegrees;
import static java.lang.String.format;
import static java.util.Locale.US;

public class Path {
	/**
	 * convert an SVG path string to an editable Object
	 */
	public static List<Command> parse(String d) {
		List<Command> path = new ArrayList<>();

		// go forward thru the String
		int i = 0;
		while (i < d.length()) {
			// get the command type
			char type = d.charAt(i);
			i ++;
			// everything up to the next commandletter is part of the arguments
			int start = i;
			while (i < d.length() && !isCommandTypeLetter(d.charAt(i)))
				i ++;
			String argString = d.substring(start, i);
			argString = argString.replaceAll("([0-9.])-", "$1,-"); //this is necessary because some Adobe products leave out delimiters between negative numbers
			argString = argString.replaceAll("\\.([0-9]+)\\.", ".$1,."); //this is necessary because some Adobe products also leave out delimiters between a number containing a decimal point and a number starting with a decimal point
			argString = argString.replaceAll("\\.([0-9]+)\\.", ".$1,."); //you have to call it twice in case there are multiple point-separated numbers in a row

			// separate the arguments
			String[] argStrings;
			if (argString.trim().isEmpty())
				argStrings = new String[0];
			else
				argStrings = argString.trim().split("[\\s,]+");

			// parse the arguments
			final double[] args = new double[argStrings.length];
			for (int j = 0; j < args.length; j++) {
				args[j] = parseDouble(argStrings[j]); //parse the coordinate
				if (!isFinite(args[j]))
					throw new IllegalArgumentException("uhh... "+type+argString);
			}

			path.add(new Command(type, args));
		}
		return path;
	}


	/**
	 * take a path and apply a scaling in the x and y directions and then shift by some amount in some direction
	 */
	public static List<Command> transformed(
			double xScale, double yScale, double xShift, double yShift, List<Command> path) {
		return translated(xShift, yShift, scaled(xScale, yScale, path));
	}


	/**
	 * shift a path by some amount in some direction
	 */
	public static List<Command> translated(double xShift, double yShift, List<Command> path) {
		List<Command> newPath = new ArrayList<>();
		for (Command old: path) {
			double[] newArgs = Arrays.copyOf(old.args, old.args.length);
			switch (old.type) {
				case 'H':
					newArgs[0] += xShift;
					break;
				case 'V':
					newArgs[0] += yShift;
					break;
				case 'L': case 'M': case 'T':
					newArgs[0] += xShift;
					newArgs[1] += yShift;
					break;
				case 'Q': case 'S':
					newArgs[0] += xShift;
					newArgs[1] += yShift;
					newArgs[2] += xShift;
					newArgs[3] += yShift;
					break;
				case 'C':
					newArgs[0] += xShift;
					newArgs[1] += yShift;
					newArgs[2] += xShift;
					newArgs[3] += yShift;
					newArgs[4] += xShift;
					newArgs[5] += yShift;
					break;
				case 'A':
					newArgs[5] += xShift;
					newArgs[6] += yShift;
					break;
				default:  // anything other than these stays the same in a translation
					break;
			}
			newPath.add(new Command(old.type, newArgs));
		}
		return newPath;
	}


	/**
	 * amplify a path about the origin in the x and y direction. the scales may be different, and
	 * either or both may be negative.
	 */
	public static List<Command> scaled(double xScale, double yScale, List<Command> path) {
		List<Command> newPath = new ArrayList<>();
		for (Command old: path) {
			double[] newArgs = Arrays.copyOf(old.args, old.args.length);
			switch (old.type) {
				case 'H': case 'h':
					newArgs[0] *= xScale;
					break;
				case 'V': case 'v':
					newArgs[0] *= yScale;
					break;
				case 'L': case 'l': case 'M': case 'm': case 'T': case 't':
					newArgs[0] *= xScale;
					newArgs[1] *= yScale;
					break;
				case 'Q': case 'q': case 'S': case 's':
					newArgs[0] *= xScale;
					newArgs[1] *= yScale;
					newArgs[2] *= xScale;
					newArgs[3] *= yScale;
					break;
				case 'C': case 'c':
					newArgs[0] *= xScale;
					newArgs[1] *= yScale;
					newArgs[2] *= xScale;
					newArgs[3] *= yScale;
					newArgs[4] *= xScale;
					newArgs[5] *= yScale;
					break;
				case 'A': case 'a':
					newArgs[0] *= abs(xScale);
					newArgs[1] *= abs(yScale);
					if (xScale*yScale < 0)
						newArgs[4] = 1 - newArgs[4];
					newArgs[5] *= xScale;
					newArgs[6] *= yScale;
					break;
				case 'Z': case 'z':
					break;
				default:
					throw new IllegalArgumentException("unrecognized path command encountered in rotated(): " + old.type);
			}
			newPath.add(new Command(old.type, newArgs));
		}
		return newPath;
	}


	/**
	 * take a path and rotate it about the origin
	 * @param rotation the angle through which to rotate it, in radians.
	 *                 positive is widdershins and negative is clockwise.
	 */
	public static List<Command> rotated(double rotation, List<Command> path) {
		List<Command> newPath = new ArrayList<>();
		for (Command old: path) {
			char newType = old.type;
			double[] newArgs = Arrays.copyOf(old.args, old.args.length);
			switch (old.type) {
				case 'H': case 'h':
					newType = (char)(old.type + 4);
					newArgs = new double[] {
							old.args[0]*cos(rotation),
							old.args[0]*sin(rotation)};
					break;
				case 'V': case 'v':
					newType = (char)(old.type - 10);
					newArgs = new double[] {
							-old.args[0]*sin(rotation),
							old.args[0]*cos(rotation)};
					break;
				case 'L': case 'l': case 'M': case 'm': case 'Q': case 'q': case 'C': case 'c': case 'T': case 't': case 'S': case 's':
					for (int i = 0; i < newArgs.length; i += 2) {
						newArgs[i    ] = old.args[i]*cos(rotation) - old.args[i + 1]*sin(rotation);
						newArgs[i + 1] = old.args[i]*sin(rotation) + old.args[i + 1]*cos(rotation);
					}
					break;
				case 'A': case 'a':
					newArgs[2] -= toDegrees(rotation);
					newArgs[5] = old.args[5]*cos(rotation) - old.args[6]*sin(rotation);
					newArgs[6] = old.args[5]*sin(rotation) + old.args[6]*cos(rotation);
					break;
				case 'Z': case 'z':
					break;
				default:
					throw new IllegalArgumentException("unrecognized path command encountered in rotated(): " + old.type);
			}
			newPath.add(new Command(newType, newArgs));
		}
		return newPath;
	}


	/**
	 * return a copy of a path that is the same but with its commands in the opposite order.  this will not do anything
	 * with the command types, so be careful, as the results will not really be valid paths.  only use this when you're
	 * planning to postprocess them somewhat afterward.
	 */
	public static List<Command> reversed(List<Command> path) {
		List<Command> reversed = new ArrayList<>(path);
		Collections.reverse(reversed);
		return reversed;
	}


	/**
	 * convert a path to an array where each row is a pair of x and y coordinates. this will lose all information
	 * regarding the separation of sections with movetos, the control points of bezier curves, and cetera.
	 * please don't pass 'H' or 'V' commands because they will be ignored.
	 */
	public static double[][] asArray(List<Command> commands) {
		double[][] points = new double[commands.size()][];
		for (int i = 0; i < commands.size(); i ++) {
			Command command = commands.get(i);
			if (command.args.length >= 2)  // just ignore 'H', 'V', and 'Z' commands
				points[i] = new double[] {
						command.args[command.args.length - 2],
						command.args[command.args.length - 1]};  // use the last two arguments as the coordinates
		}
		return points;
	}


	private static boolean isCommandTypeLetter(char c) {
		return (c >= 'A' && c <= 'Z' && c != 'E') || (c >= 'a' && c <= 'z' && c != 'e');
	}


	private static String formatDouble(double d) { //format a number with the minimum number of digits, but at most 3
		String str = format(US, "%.3f", d);
		if (str.contains("."))
			str = str.replaceFirst("0+$", "");
		if (str.endsWith("."))
			str = str.substring(0, str.length() - 1);
		return str;
	}

	/**
	 * express the path as a d string as would be used in an SVG &lt;path&gt; element
	 */
	public static String toString(List<Command> path) {
		StringBuilder d = new StringBuilder();
		for (Path.Command c: path)
			d.append(c.toString()).append(" ");
		return d.toString();
	}


	/**
	 * An SVG path command, like line or bezier curve or whatever
	 * @author jkunimune
	 */
	public static class Command {
		/** 'M' for moveto, 'L' for line, 'C' for cubic, etc. */
		final public char type;
		/** the coordinates and parameters according to the type of command */
		final public double[] args;

		public Command(char type, double... args) {
			this.type = type;
			this.args = args;
		}

		public String toString() {
			StringBuilder s = new StringBuilder(Character.toString(type));
			for (double arg : args)
				s.append(formatDouble(arg)).append(",");
			return s.substring(0, max(1, s.length()-1));
		}
	}
}
