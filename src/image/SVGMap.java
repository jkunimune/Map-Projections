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
package image;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;
import org.xml.sax.helpers.DefaultHandler;
import utils.BoundingBox;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import static java.lang.Double.NaN;
import static java.lang.Double.isFinite;
import static java.lang.Double.isNaN;
import static java.lang.Double.parseDouble;
import static java.lang.Math.PI;
import static java.lang.Math.hypot;
import static java.lang.Math.max;
import static java.lang.String.format;
import static utils.Math2.linInterp;

/**
 * An input equirectangular map based on an SVG file
 * 
 * @author jkunimune
 */
public class SVGMap implements Iterable<SVGMap.SVGElement>, SavableImage {
	
	public static final double[] NULL_TRANSFORM = {1, 1, 0, 0};
	
	private static final double MAX_EDGE_LENGTH = 1/20.; // cut lines that are more than this far across the map

	private SVGHeader header; // the initial element, which contains dimension information
	private final List<SVGElement> elements; //all of the parsed file contents
	private final int numVertices; //the total number of path commands, for optimization purposes


	/**
	 * construct an SVGMap object from an SVG file. this will parse the XML syntax, break the file
	 * structure up into SVGElement objects containing all of the image data, and linearly map that
	 * image data from whatever viewBox onto the [-PI, PI], [-PI/2, PI/2] domain.
	 * @throws IOException if there's trouble accessing the ifle
	 * @throws SAXException if the XML syntax is invalid
	 * @throws ParserConfigurationException if the SAXParser object can't be constructed for some reason
	 */
	public SVGMap(File file) throws IOException, SAXException, ParserConfigurationException {
		elements = new LinkedList<>(); // the list elements

		final SAXParser parser = SAXParserFactory.newInstance().newSAXParser();
		
		final DefaultHandler handler = new DefaultHandler() {
			private final Deque<double[]> transformStack =
					new ArrayDeque<double[]>(Collections.singleton(NULL_TRANSFORM));

			@Override
			public InputSource resolveEntity(String publicId, String systemId) {
				return new InputSource(new StringReader("")); //ignore all external references - we don't need to validate
			}

			@Override
			public void startElement(
					String uri, String localName, String tagName, Attributes immutableAttributes
			) {
				AttributesImpl attributes = new AttributesImpl(immutableAttributes);

				// first parse any transformation embedded in this element
				if (attributes.getIndex("transform") >= 0) {
					transformStack.push(parseTransform(attributes.getValue("transform")));
					attributes.removeAttribute(attributes.getIndex("transform"));
				}
				else {
					transformStack.push(transformStack.peek());
				}

				// pull out width, height, and viewBox from <svg>
				if (tagName.equals("svg"))
					elements.add(parseSVGHeader(attributes));
				// pull out the d from <path>
				else if (tagName.equals("path"))
					elements.add(parsePath(attributes));
				// convert certain rectangles to Boundary objects
				else if (tagName.equals("rect")) {
					Background parsedBackground = parseBackground(attributes);
					if (parsedBackground != null)
						elements.add(parsedBackground);
					else
						elements.add(parsePoint("rect", attributes, "x", "y"));
				}
				// pull out the x and y from <text>
				else if (attributes.getIndex("x") >= 0 && attributes.getIndex("y") >= 0)
					elements.add(parsePoint(tagName, attributes, "x", "y"));
				// pull out the cx and cy from <circle>
				else if (attributes.getIndex("cx") >= 0 && attributes.getIndex("cy") >= 0)
					elements.add(parsePoint(tagName, attributes, "cx", "cy"));
				// for everything else, just treat it as an invariant String
				else
					elements.add(new Content(formatAttributes(tagName, attributes)));
			}

			@Override
			public void endElement(String uri, String localName, String qName) {
				// save the </> tag as an invariant string
				elements.add(new Content(format("</%s>", qName)));
				// and drop any transform that may have been been on that element
				transformStack.pop();
			}

			@Override
			public void characters(char[] ch, int start, int length) {
				// just dump the provided characters into currentFormatString
				StringBuilder characterString = new StringBuilder();
				for (int i = 0; i < length; i++) {
					characterString.append(ch[start + i]);
				}
				elements.add(new Content(characterString.toString()));
			}

			private double[] parseTransform(String transString) {
				double xScale = 1, yScale = 1, xTrans = 0, yTrans = 0;
				int i;
				while ((i = transString.indexOf('(')) >= 0) {
					int j = transString.indexOf(')');
					String type = transString.substring(0, i).trim();
					String argString = transString.substring(i + 1, j);
					String[] args = argString.split("[,\\s]+");
					switch (type) {
						case "matrix":
							xScale = parseDouble(args[0]);
							yScale = parseDouble(args[3]);
							xTrans = parseDouble(args[4]);
							yTrans = parseDouble(args[5]);
							break;
						case "translate":
							if (args.length == 1) {
								xTrans = yTrans = parseDouble(args[0]);
							}
							else {
								xTrans = parseDouble(args[0]);
								yTrans = parseDouble(args[1]);
							}
							break;
						case "scale":
							if (args.length == 1) {
								xScale = yScale = parseDouble(args[0]);
							}
							else {
								xScale = parseDouble(args[0]);
								yScale = parseDouble(args[1]);
							}
							break;
					} // I'm not worrying about shear and rotation because I don't want to
					transString = transString.substring(j + 1);
				}
				return new double[]{xScale, yScale, xTrans, yTrans};
			}

			private SVGHeader parseSVGHeader(AttributesImpl attributes) {
				if (attributes.getValue("width") == null)
					attributes.addAttribute("", "", "width", "", "360");
				if (attributes.getValue("height") == null)
					attributes.addAttribute("", "", "height", "", "180");
				double displayWidth = parseDouble(attributes.getValue("width"));
				double displayHeight = parseDouble(attributes.getValue("height"));
				if (attributes.getValue("viewBox") == null)
					attributes.addAttribute("", "", "viewBox", "",
					                        format("%f %f %f %f", 0., 0., displayWidth, displayHeight));
				String[] values = attributes.getValue("viewBox").split("\\s", 4);
				double vbMinX = parseDouble(values[0]);
				double vbMinY = parseDouble(values[1]);
				double vbWidth = parseDouble(values[2]);
				double vbHeight = parseDouble(values[3]);

				attributes.setValue(attributes.getIndex("width"), "%1$.6g");
				attributes.setValue(attributes.getIndex("height"), "%2$.6g");
				attributes.setValue(attributes.getIndex("viewBox"), "%3$.6g %4$.6g %5$.6g %6$.6g");
				String formatSpecifier = formatAttributes("svg", attributes);

				header = new SVGHeader(formatSpecifier, displayWidth, displayHeight, vbMinX, vbMinY, vbWidth, vbHeight);
				return header;
			}

			private Background parseBackground(AttributesImpl attributes) {
				if (attributes.getValue("x") == null || attributes.getValue("y") == null ||
				    attributes.getValue("width") == null || attributes.getValue("height") == null)
					return null;
				double x = parseDouble(attributes.getValue("x"));
				double y = parseDouble(attributes.getValue("y"));
				double width = parseDouble(attributes.getValue("width"));
				double height = parseDouble(attributes.getValue("height"));
				// apply the coordinate transformation
				double[] transform = transformStack.peek();
				if (transform == null)
					throw new IllegalArgumentException("there will always be a nonnull transform to peek at");
				x = x*transform[0] + transform[2];
				y = y*transform[1] + transform[3];
				width = width*transform[0];
				height = height*transform[1];
				// check if it exactly covers the viewbox (it's okay if it's vertically inverted)
				if ((x == header.vbMinX && y == header.vbMinY && width == header.vbWidth && height == header.vbHeight) ||
				    (x == header.vbMinX && y == header.vbMinY + header.vbHeight && width == header.vbWidth && height == -header.vbHeight))
					return new Background(attributes, new BoundingBox(x, x + width, y, y + height));
				else
					return null;
			}

			private Point parsePoint(String tagName, AttributesImpl attributes, String xName, String yName) {
				if (attributes.getValue(xName) == null || attributes.getValue(yName) == null)
					throw new IllegalArgumentException(format("this <%s> seems to be missing some important parameters; where's %s and %s?", tagName, xName, yName));
				// parse the coordinates
				double x = parseDouble(attributes.getValue(xName));
				attributes.setValue(attributes.getIndex(xName), "%1$.6g");
				double y = parseDouble(attributes.getValue(yName));
				attributes.setValue(attributes.getIndex(yName), "%2$.6g");
				String formatSpecifier = formatAttributes(tagName, attributes);
				// apply the coordinate transformations
				double[] transform = transformStack.peek();
				if (transform == null)
					throw new IllegalArgumentException("there will always be a nonnull transform to peek at");
				x = x*transform[0] + transform[2]; //apply the transformation
				y = y*transform[1] + transform[3];
				x = linInterp(x, header.vbMinX, header.vbMinX + header.vbWidth, -PI, PI); //scale to radians
				y = linInterp(y, header.vbMinY + header.vbHeight, header.vbMinY, -PI/2, PI/2);
				// put it all together and return it
				return new Point(formatSpecifier, x, y);
			}

			private Path parsePath(AttributesImpl attributes) {
				// start by extracting the "d" attribute
				String d;
				if (attributes.getValue("d") != null) {
					d = attributes.getValue("d");
					attributes.setValue(attributes.getIndex("d"), "%1$s");
				}
				else
					d = "";

				// form the format specifier from the other attributes
				String formatSpecifier = formatAttributes("path", attributes);

				// then parse the d
				List<Command> commands = new ArrayList<>();
				int i = 0;
				double[] lastMove = {0, 0}; //for closepaths
				double[] last = {0, 0}; //for relative coordinates
				while (i < d.length()) {
					char type = d.charAt(i);
					i ++;
					int start = i;
					while (i < d.length() && !isCommandTypeLetter(d.charAt(i)))
						i ++;
					String argString = d.substring(start, i);
					argString = argString.replaceAll("([0-9.])-", "$1,-"); //this is necessary because some Adobe products leave out delimiters between negative numbers

					String[] argStrings;
					if (argString.trim().isEmpty()) 	argStrings = new String[0];
					else 								argStrings = argString.trim().split("[\\s,]+");
					final double[] args;

					// replaceInText arcs with lines because I don't want to deal with those
					if (type == 'a' || type == 'A') {
						type += (type == 'a') ? 'l' : 'L';
						argStrings = new String[] {argStrings[3], argStrings[4]};
					}
					// convert horizontal and vertical lines to general lines to keep things simple
					if (type == 'h' || type == 'H' || type == 'v' || type == 'V') {
						final int direcIdx = (type%32 == 8) ? 0 : 1;
						args = new double[] {last[0], last[1]};
						if (type <= 'Z') 	args[direcIdx] = parseDouble(argStrings[0]); //uppercase (absolute)
						else 				args[direcIdx] += parseDouble(argStrings[0]); //lowercase (relative)
						last[direcIdx] = args[direcIdx];
						type = 'L';
					}
					// also replaceInText closepath with an explicit line to the last moveto
					else if (type == 'z' || type == 'Z') {
						args = new double[] {lastMove[0], lastMove[1]};
						type = 'L';
					}
					// now that everything is a curve or a line, extract the coordinates
					else {
						args = new double[argStrings.length];
						for (int j = 0; j < args.length; j ++) {
							args[j] = parseDouble(argStrings[j]); //parse the coordinate

							if (type >= 'a')
								args[j] += last[j%2]; //account for relative commands
							last[j%2] = args[j];
						}
						if (type >= 'a') //make all letters uppercase (because it's all absolute from now on)
							type -= 32;
					}
					if (type == 'M') { //make note of the last moveto, so we can interpret closepaths properly
						lastMove[0] = args[args.length-2];
						lastMove[1] = args[args.length-1];
					}

					double[] transform = transformStack.peek();
					if (transform == null)
						throw new IllegalArgumentException("there will always be a nonnull transform to peek at");
					for (int j = 0; j < args.length; j ++) {
						if (!isFinite(args[j]))
							throw new IllegalArgumentException("uhh... "+type+argString);
						if (j%2 == 0) {
							args[j] = args[j]*transform[0] + transform[2]; //apply the transformation
							args[j] = linInterp(args[j], header.vbMinX, header.vbMinX+header.vbWidth,
							                    -PI, PI); //scale to radians
						}
						else {
							args[j] = args[j]*transform[1] + transform[3];
							args[j] = linInterp(args[j], header.vbMinY+header.vbHeight, header.vbMinY, //keep in mind that these are paired longitude-latitude
							                    -PI/2, PI/2); //not latitude-longitude, as they are elsewhere
						}
					}
					commands.add(new Command(type, args));
				}

				return new Path(formatSpecifier, commands);
			}
		};

		parser.parse(new BufferedInputStream(Files.newInputStream(file.toPath())), handler);

		if (header == null)
			throw new IllegalArgumentException("the SVG didn't start with a header");

		int numVertices = 0;
		for (SVGElement element: elements)
			if (element instanceof Path)
				numVertices += ((Path)element).commands.size();
		this.numVertices = numVertices;
	}


	/**
	 * construct an SVG map by passing a preparsed list of SVGElements
	 */
	public SVGMap(List<SVGElement> elements, int numVertices) {
		this.elements = elements;
		this.numVertices = numVertices;
		for (SVGElement element: elements)
			if (element instanceof SVGHeader)
				this.header = (SVGHeader) element;
	}
	
	
	
	public int getNumVertices() {
		return numVertices;
	}


	public int getNumElements() {
		return elements.size();
	}
	
	
	@Override
	public Iterator<SVGElement> iterator() {
		return elements.iterator();
	}
	
	
	/**
	 * write this SVG image to a file
	 * @throws IOException if it has trouble writing to the file
	 */
	public void save(File file) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(file));

		for (SVGElement element: this)
			out.write(element.toString());
		out.close();
	}


	/**
	 * @return the left boundary of the viewBox
	 */
	public double getVbMinX() {
		return this.header.vbMinX;
	}

	/**
	 * @return the upper boundary of the viewBox
	 */
	public double getVbMinY() {
		return this.header.vbMinY;
	}

	/**
	 * @return the right boundary of the viewBox
	 */
	public double getVbMaxX() {
		return this.header.vbMinX + this.header.vbWidth;
	}

	/**
	 * @return the lower boundary of the viewBox
	 */
	public double getVbMaxY() {
		return this.header.vbMinY + this.header.vbHeight;
	}

	/**
	 * @return the viewBox width
	 */
	public double getVbWidth() {
		return this.header.vbWidth;
	}

	/**
	 * @return the viewBox height
	 */
	public double getVbHeight() {
		return this.header.vbHeight;
	}

	/**
	 * @return the SVG's display width
	 */
	public double getDisplayWidth() {
		return this.header.width;
	}

	/**
	 * @return the SVG's display height
	 */
	public double getDisplayHeight() {
		return this.header.height;
	}


	public static Path breakWraps(Path continuous, double inSize) { //break excessively long commands, as they are likely wrapping over a discontinuity
		if (continuous.commands.size() <= 2) 	return continuous;
		List<Command> broken = new LinkedList<>();
		final double lengthThreshold = inSize*MAX_EDGE_LENGTH;
		double[] lens = {NaN, NaN, NaN}; //the revolving array of command lengths
		for (int i = 0; i < continuous.commands.size(); i ++) {
			if (i < continuous.commands.size()-1 && continuous.commands.get(i+1).type != 'M')
				lens[2] = hypot( //compute this next length
						continuous.commands.get(i+1).args[0] - continuous.commands.get(i).args[0],
						continuous.commands.get(i+1).args[1] - continuous.commands.get(i).args[1]);
			else
				lens[2] = NaN;
			
			char type = continuous.commands.get(i).type;
			if (lens[1] >= lengthThreshold && // check it against an absolute threshold
					(isNaN(lens[0]) || lens[1] > 20*lens[0]) && //and compare it to the last two lengths
					(isNaN(lens[2]) || lens[1] > 20*lens[2])) //if both sides are far longer or nonexistent
				type = 'M'; //break this line
			
			broken.add(new Command(type, continuous.commands.get(i).args.clone()));
			lens[0] = lens[1];
			lens[1] = lens[2];
		}
		
		return new Path(continuous.formatSpecifier, broken);
	}
	
	
	public static Path closePaths(Path open) { //replaceInText plain loops with 'Z's and combine connected parts
		if (open.commands.size() <= 1) 	return open;
		List<List<Command>> parts = new ArrayList<>();
		List<Command> currentPart = null;
		for (Command cmd: open.commands) { //start by breaking the Path into parts,
			if (cmd.type == 'M') { //separated by movetos
				if (currentPart != null)
					parts.add(currentPart);
				currentPart = new LinkedList<Command>();
			}
			else if (currentPart == null)
				throw new RuntimeException(format(
						"this path should start with 'M' (and we should have caught that by now): '%s'", open));
			currentPart.add(cmd);
		}
		parts.add(currentPart);
		
		List<Command> closed = new LinkedList<Command>();
		for (int i = 0; i < parts.size(); i ++) { //now look through those parts
			List<Command> partI = parts.get(i);
			if (partI.size() > 1
					&& Arrays.equals(partI.get(0).args, partI.get(partI.size()-1).args)) { //if it is self-enclosing
				partI.set(partI.size()-1, new Command('Z', new double[0])); //give it a closepath and send it on its way
			}
			else { //if it is open
				for (int j = i+1; j < parts.size(); j ++) { //look to see if there is anything that completes it
					List<Command> partJ = parts.get(j);
					if (Arrays.equals(partI.get(0).args, partJ.get(partJ.size()-1).args)) { //if so,
						partI.remove(0); //remove the useless moveto
						partJ.addAll(partI); //combine them
						partI = partJ;
						parts.remove(j); //don't look at J anymone; it has been absorbed.
						break;
					}
				}
			}
			closed.addAll(partI); //now turn in whatever you've got
		}
		return new Path(open.formatSpecifier, closed);
	}
	
	
	private static boolean isCommandTypeLetter(char c) {
		return (c >= 'A' && c <= 'Z' && c != 'E') || (c >= 'a' && c <= 'z' && c != 'e');
	}


	private static String formatDouble(double d) { //format a number with the minimum number of digits, but at most 3
		String str = format("%.3f", d);
		if (str.contains("."))
			str = str.replaceFirst("0+$", "");
		if (str.endsWith("."))
			str = str.substring(0, str.length() - 1);
		return str;
	}


	private static String formatAttributes(String tagName, Attributes attributes) {
		StringBuilder string = new StringBuilder("<").append(tagName);
		for (int i = 0; i < attributes.getLength(); i++)
			string.append(" ").append(attributes.getQName(i)).append("=\"").append(attributes.getValue(i)).append("\"");
		string.append(">");
		return string.toString();
	}


	/**
	 * a representation of some bytes from an SVG file
	 */
	public interface SVGElement {
		String toString();
	}

	/**
	 * some generic bytes from an SVG file represented as a simple String
	 */
	public static class Content implements SVGElement {
		public final String content;

		public Content(String content) {
			this.content = content;
		}

		public Content projected() {
			return this;
		}

		public String toString() {
			return this.content;
		}
	}

	/**
	 * the &lt;svg ...&gt; tag at the top of the SVG file, with dimension information extracted
	 */
	public static class SVGHeader implements SVGElement {
		public final String formatSpecifier;
		public final double width;
		public final double height;
		public final double vbMinX;
		public final double vbMinY;
		public final double vbWidth;
		public final double vbHeight;

		public SVGHeader(String formatSpecifier, double width, double height, double vbMinX, double vbMinY, double vbWidth, double vbHeight) {
			this.formatSpecifier = formatSpecifier;
			this.width = width;
			this.height = height;
			this.vbMinX = vbMinX;
			this.vbMinY = vbMinY;
			this.vbWidth = vbWidth;
			this.vbHeight = vbHeight;
		}

		public String toString() {
			return format(formatSpecifier, width, height, vbMinX, vbMinY, vbWidth, vbHeight);
		}
	}

	/**
	 * a shape that represents not a particular geographic feature but the outline of the full world map
	 */
	public static class Background implements SVGElement {
		public final BoundingBox bounds;
		public final AttributesImpl attributes;

		public Background(AttributesImpl attributes, BoundingBox bounds) {
			this.bounds = bounds;
			this.attributes = attributes;
		}

		public String toString() {
			AttributesImpl fullAttributes = new AttributesImpl(attributes);
			fullAttributes.setValue(fullAttributes.getIndex("x"), Double.toString(bounds.xMin));
			fullAttributes.setValue(fullAttributes.getIndex("y"), Double.toString(bounds.yMin));
			fullAttributes.setValue(fullAttributes.getIndex("width"), Double.toString(bounds.width));
			fullAttributes.setValue(fullAttributes.getIndex("height"), Double.toString(bounds.height));
			return formatAttributes("rect", fullAttributes);
		}
	}

	/**
	 * an SVG tag containing a single pair of coordinates that can be projected
	 */
	public static class Point implements SVGElement {
		public final String formatSpecifier;
		public final double x;
		public final double y;

		public Point(String formatSpecifier, double x, double y) {
			this.formatSpecifier = formatSpecifier;
			this.x = x;
			this.y = y;
		}

		public Point projected() {
			return this;
		}

		public String toString() {
			return format(formatSpecifier, x, y);
		}
	}


	/**
	 * a &lt;path d=...&gt; tag, with the d String extracted and parsed into a modifiable form
	 */
	public static class Path implements SVGElement {
		public final String formatSpecifier;
		public final List<Command> commands;

		public Path(String formatSpecifier, List<Command> commands) {
			this.formatSpecifier = formatSpecifier;
			this.commands = commands;
		}
		
		public String toString() {
			StringBuilder s = new StringBuilder();
			for (Command c: commands)
				s.append(c.toString()).append(" ");
			return format(formatSpecifier, s.toString().trim());
		}
	}
	
	
	/**
	 * An SVG path command, like line or bezier curve or whatever
	 * @author jkunimune
	 */
	public static class Command {
		final public char type; //M, L, C, etc.
		final public double[] args; //the coordinates that go with it
		
		public Command(char type, double[] args) {
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
