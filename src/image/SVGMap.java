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
import utils.Quantity;

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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

import static java.lang.Double.NaN;
import static java.lang.Double.parseDouble;
import static java.lang.Math.PI;
import static java.lang.Math.hypot;
import static java.lang.String.format;
import static java.util.Locale.US;
import static utils.Math2.linInterp;
import static utils.Math2.max;
import static utils.Quantity.parseQuantity;

/**
 * An input equirectangular map based on an SVG file
 * 
 * @author jkunimune
 */
public class SVGMap implements Iterable<SVGMap.SVGElement>, SavableImage {
	
	public static final double[] NULL_TRANSFORM = {1, 1, 0, 0};

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
		elements = new ArrayList<>(); // the list elements
		elements.add(new Content( // start by guessing the xml declaration, since SaxParser always skips it
				"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"));

		final SAXParser parser = SAXParserFactory.newInstance().newSAXParser();
		
		final DefaultHandler handler = new DefaultHandler() {
			private final Stack<String> tagStack = new Stack<String>();
			private final Stack<double[]> transformStack = new Stack<double[]>();

			@Override
			public void startDocument() {
				transformStack.push(NULL_TRANSFORM);
			}

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

				SVGElement element;
				// pull out width, height, and viewBox from <svg>
				if (tagName.equals("svg"))
					element = parseSVGHeader(attributes);
				// pull out the d from <path>
				else if (tagName.equals("path"))
					element = parsePath(attributes);
				// convert certain rectangles to Boundary objects
				else if (tagName.equals("rect")) {
					element = parseBackground(attributes);
					if (element != null)
						tagName = "path";
					else
						element = parsePoint("rect", attributes, "x", "y");
				}
				// pull out the x and y from <text>
				else if (attributes.getIndex("x") >= 0 && attributes.getIndex("y") >= 0)
					element = parsePoint(tagName, attributes, "x", "y");
				// pull out the cx and cy from <circle>
				else if (attributes.getIndex("cx") >= 0 && attributes.getIndex("cy") >= 0)
					element = parsePoint(tagName, attributes, "cx", "cy");
				// for everything else, just treat it as an invariant String
				else
					element = new Content(formatAttributes(tagName, attributes));

				elements.add(element);
				tagStack.add(tagName);
			}

			@Override
			public void endElement(String uri, String localName, String qName) {
				// save the </> tag as an invariant string
				elements.add(new Content(format("</%s>", tagStack.pop())));
				// and drop any transform that may have been been on that element
				transformStack.pop();
			}

			@Override
			public void characters(char[] characters, int start, int length) {
				// just dump the provided characters into currentFormatString
				elements.add(new Content(escapeHTML(
						new String(characters).substring(start, start + length))));
			}
			
//			@Override  /* only available in Java 14+ */
//			public void declaration(String version, String encoding, String standalone) {
//				elements.add(new Content(format(
//						"<?xml version=\"%s\" encoding=\"%s\" standalone=\"%s\"?>\n",
//						version, encoding, standalone)));
//			}
			
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
				Quantity displayWidth = parseQuantity(attributes.getValue("width"));
				Quantity displayHeight = parseQuantity(attributes.getValue("height"));
				if (attributes.getValue("viewBox") == null)
					attributes.addAttribute(
							"", "", "viewBox", "",
							format(US, "%f %f %f %f",
							       0., 0., displayWidth.value, displayHeight.value));
				String[] values = attributes.getValue("viewBox").split("\\s", 4);
				double vbMinX = parseDouble(values[0]);
				double vbMinY = parseDouble(values[1]);
				double vbWidth = parseDouble(values[2]);
				double vbHeight = parseDouble(values[3]);

				attributes.setValue(attributes.getIndex("width"), "%1$.6g%2$s");
				attributes.setValue(attributes.getIndex("height"), "%3$.6g%4$s");
				attributes.setValue(attributes.getIndex("viewBox"), "%5$.6g %6$.6g %7$.6g %8$.6g");
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
				// check if it exactly covers the viewBox (it's okay if it's vertically inverted)
				if ((x == header.vbMinX && y == header.vbMinY && width == header.vbWidth && height == header.vbHeight) ||
				    (x == header.vbMinX && y == header.vbMinY + header.vbHeight && width == header.vbWidth && height == -header.vbHeight)) {
					// if so, remove the dimensions and stick the remainder in a Background Object
					attributes = new AttributesImpl(attributes);
					attributes.removeAttribute(attributes.getIndex("x"));
					attributes.removeAttribute(attributes.getIndex("y"));
					attributes.removeAttribute(attributes.getIndex("width"));
					attributes.removeAttribute(attributes.getIndex("height"));
					return new Background(attributes, null);
				}
				else
					return null;
			}

			private GeographicPoint parsePoint(String tagName, AttributesImpl attributes, String xName, String yName) {
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
				return new GeographicPoint(formatSpecifier, x, y);
			}

			private GeographicPath parsePath(AttributesImpl attributes) {
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

				// then parse the d and simplify the result
				List<Path.Command> path = Path.parse(d);
				
				// start by breaking up any commands with multiple sets of arguments
				for (int i = path.size() - 1; i >= 0; i --) {
					char type = path.get(i).type;
					int expectedArgsLength;
					if (type == 'Z' || type == 'z')
						expectedArgsLength = 0;
					else if (type == 'A' || type == 'a')
						expectedArgsLength = 7;
					else if (type == 'C' || type == 'c')
						expectedArgsLength = 6;
					else if (type == 'Q' || type == 'q' || type == 'S' || type == 's')
						expectedArgsLength = 4;
					else if (type == 'H' || type == 'h' || type == 'V' || type == 'v')
						expectedArgsLength = 1;
					else
						expectedArgsLength = 2;
					double[] args = path.get(i).args;
					if (args.length < expectedArgsLength)
						throw new IllegalArgumentException(format(
								"this '%s' command only has %d elements, but I think '%s's should have %d.",
								type, args.length, type, expectedArgsLength));
					if (args.length > expectedArgsLength) {
						if (expectedArgsLength == 0 || args.length%expectedArgsLength != 0)
							throw new IllegalArgumentException(format(
									"this '%s' command has %d elements, which is not a whole multiple of %d like it should be.",
									type, args.length, expectedArgsLength));
						path.remove(i);
						for (int j = 0; j < args.length/expectedArgsLength; j ++) {
							if (j > 0 && type == 'M')
								type = 'L'; //when multiple pairs of args come after 'M', the subsequent commands are 'L'
							double[] subArgs = Arrays.copyOfRange(
									args, j*expectedArgsLength, (j + 1)*expectedArgsLength);
							path.add(i + j, new Path.Command(type, subArgs));
						}
					}
				}
				
				// then convert fancy commands into simpler commands
				double[] lastMove = {0, 0}; //for closepaths
				double[] lastPoint = {0, 0}; //for H, V, and relative coordinates
				double[] lastControlPoint = {0, 0}; // for T and S commands
				for (int i = 0; i < path.size(); i ++) {
					char type = path.get(i).type; //get the type
					double[] args = path.get(i).args;
					
					// replace relative commands with absolute commands
					if (type >= 'a') {
						if (type == 'h') //for h, it's relative to the last x coordinate
							args[0] += lastPoint[0];
						else if (type == 'v') //for v, it's relative to the last y coordinate
							args[0] += lastPoint[1];
						else if (type == 'a') { //for arcs, only the last two args are relative
							args[5] += lastPoint[0];
							args[6] += lastPoint[1];
						}
						else { //for everything else, shift all the args in pairs
							for (int j = 0; j < args.length; j += 2) {
								args[j] += lastPoint[0];
								args[j + 1] += lastPoint[1];
							}
						}
						type -= 32; //and then make it a capital letter to indicate that it is now absolute
					}
					// replace arcs with lines because I don't want to deal with those
					if (type == 'A') {
						type = 'L';
						args = new double[] {args[5], args[6]};
					}
					// convert horizontal and vertical lines to general lines to keep things simple
					if (type == 'H') {
						args = new double[] {args[0], lastPoint[1]};
						type = 'L';
					}
					if (type == 'V') {
						args = new double[] {lastPoint[0], args[0]};
						type = 'L';
					}
					// also replace closepath with an explicit line to the last moveto
					else if (type == 'Z') {
						args = new double[] {lastMove[0], lastMove[1]};
						type = 'L';
					}
					// also convert 'T' Bezier curves to explicit 'Q' Bezier curves
					else if (type == 'T') {
						args = new double[] {
							2*lastPoint[0] - lastControlPoint[0],
							2*lastPoint[1] - lastControlPoint[1],
							args[0], args[1],
						};
						type = 'Q';
					}
					// also convert 'S' Bezier curves to explicit 'C' Bezier curves
					else if (type == 'S') {
						args = new double[] {
							2*lastPoint[0] - lastControlPoint[0],
							2*lastPoint[1] - lastControlPoint[1],
							args[0], args[1],
							args[2], args[3],
						};
						type = 'C';
					}
					lastPoint = Arrays.copyOfRange(args, args.length - 2, args.length); //make note of the endpoint of this command
					if (type == 'M') //make note of the last moveto, so we can interpret closepaths properly
						lastMove = lastPoint;
					if (type == 'Q' || type == 'C') // make note of the last control point, if there were any
						lastControlPoint = Arrays.copyOfRange(args, args.length - 4, args.length - 2);
					else
						lastControlPoint = lastPoint;
					path.set(i, new Path.Command(type, args));
				}

				// apply the transform
				double[] transform = transformStack.peek();
				if (transform == null)
					throw new IllegalArgumentException("there will always be a nonnull transform to peek at");
				path = Path.transformed(transform[0], transform[1], transform[2], transform[3], path);
				path = Path.translated(-header.vbMinX, -header.vbMinY, path);
				path = Path.transformed(2*PI/header.vbWidth, -PI/header.vbHeight, -PI, PI/2, path);

				return new GeographicPath(formatSpecifier, path);
			}
		};

		parser.parse(new BufferedInputStream(Files.newInputStream(file.toPath())), handler);

		if (header == null)
			throw new IllegalArgumentException("the SVG didn't start with a header");

		int numVertices = 0;
		for (SVGElement element: elements)
			if (element instanceof GeographicPath)
				numVertices += ((GeographicPath)element).commands.size();
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
	public Quantity getDisplayWidth() {
		return this.header.width;
	}

	/**
	 * @return the SVG's display height
	 */
	public Quantity getDisplayHeight() {
		return this.header.height;
	}


	/**
	 * modify a path to heuristicly remove line erroneus line segments caused by the path crossing an interruption.
	 * it looks for lines that are too long compared to the map size and compared to the adjacent segments.
	 * @param continuous the path containing erroneus segments to be removed.
	 * @param inSize the total length scale of the map
	 * @param strict if strict is true, ony lines with a sufficient absolute length and relative length will be cut.
	 *               if it is not, a really long relative or absolute length alone is sufficient for cutting.
	 * @return the path containing breaks
	 */
	public static GeographicPath breakWraps(GeographicPath continuous, double inSize, boolean strict) { //break excessively long commands, as they are likely wrapping over a discontinuity
		if (continuous.commands.size() <= 2) 	return continuous;
		List<Path.Command> broken = new ArrayList<>();
		double[] lens = {0, 0, 0}; //the revolving array of command lengths
		for (int i = 0; i < continuous.commands.size(); i ++) {
			if (i < continuous.commands.size()-1 && continuous.commands.get(i+1).type != 'M')
				lens[2] = hypot( //compute this next length
						continuous.commands.get(i+1).args[0] - continuous.commands.get(i).args[0],
						continuous.commands.get(i+1).args[1] - continuous.commands.get(i).args[1]);
			else
				lens[2] = NaN;
			
			char type = continuous.commands.get(i).type;
			// check it against both an absolute threshold and its two neighbor lengths
			boolean tooLong = lens[1] >= inSize*0.05 && lens[1] >= max(lens[0], lens[2])*20;
			if (!strict)
				tooLong |= lens[1] >= inSize*0.25;
			if (tooLong) //if it's suspiciusly long
				type = 'M'; //break this line
			
			broken.add(new Path.Command(type, continuous.commands.get(i).args.clone()));
			lens[0] = lens[1];
			lens[1] = lens[2];
		}

		return new GeographicPath(continuous.formatSpecifier, broken);
	}
	
	
	public static GeographicPath closePaths(GeographicPath open) { //replace plain loops with 'Z's and combine connected parts
		if (open.commands.size() <= 1) 	return open;
		List<List<Path.Command>> parts = new ArrayList<>();
		List<Path.Command> currentPart = null;
		for (Path.Command cmd: open.commands) { //start by breaking the Path into parts,
			if (cmd.type == 'M') { //separated by movetos
				if (currentPart != null)
					parts.add(currentPart);
				currentPart = new ArrayList<Path.Command>();
			}
			else if (currentPart == null)
				throw new RuntimeException(format(
						"this path should start with 'M' (and we should have caught that by now): '%s'", open));
			currentPart.add(cmd);
		}
		parts.add(currentPart);
		
		List<Path.Command> closed = new ArrayList<Path.Command>();
		for (int i = 0; i < parts.size(); i ++) { //now look through those parts
			List<Path.Command> partI = parts.get(i);
			if (partI.size() > 1
					&& Arrays.equals(partI.get(0).args, partI.get(partI.size()-1).args)) { //if it is self-enclosing
				partI.set(partI.size()-1, new Path.Command('Z')); //give it a closepath and send it on its way
			}
			else { //if it is open
				for (int j = i+1; j < parts.size(); j ++) { //look to see if there is anything that completes it
					List<Path.Command> partJ = parts.get(j);
					if (Arrays.equals(partI.get(0).args, partJ.get(partJ.size()-1).args)) { //if so,
						partI.remove(0); //remove the useless moveto
						partJ.addAll(partI); //combine them
						partI = partJ;
						parts.remove(j); //don't look at J anymore; it has been absorbed.
						break;
					}
				}
			}
			closed.addAll(partI); //now turn in whatever you've got
		}
		return new GeographicPath(open.formatSpecifier, closed);
	}
	
	
	private static String formatAttributes(String tagName, Attributes attributes) {
		StringBuilder string = new StringBuilder("<").append(tagName);
		for (int i = 0; i < attributes.getLength(); i++)
			string.append(" ").append(attributes.getQName(i)).append("=\"").append(attributes.getValue(i)).append("\"");
		string.append(">");
		return string.toString();
	}
	
	
	public static String escapeHTML(String raw) {
		StringBuilder escaped = new StringBuilder();
		for (int i = 0; i < raw.length(); i ++) {
			char c = raw.charAt(i);  // don't forget to escape special characters
			if (c >= 128 || c == '\'' || c == '"' || c == '<' || c == '>' || c == '&') // some characters must be escaped here
				escaped.append(format("&#%d;", (int) c));
			else
				escaped.append(c);
		}
		return escaped.toString();
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

		public String toString() {
			return this.content;
		}
	}

	/**
	 * the &lt;svg ...&gt; tag at the top of the SVG file, with dimension information extracted
	 */
	public static class SVGHeader implements SVGElement {
		public final String formatSpecifier;
		public final Quantity width;
		public final Quantity height;
		public final double vbMinX;
		public final double vbMinY;
		public final double vbWidth;
		public final double vbHeight;

		public SVGHeader(String formatSpecifier, Quantity width, Quantity height,
		                 double vbMinX, double vbMinY, double vbWidth, double vbHeight) {
			this.formatSpecifier = formatSpecifier;
			this.width = width;
			this.height = height;
			this.vbMinX = vbMinX;
			this.vbMinY = vbMinY;
			this.vbWidth = vbWidth;
			this.vbHeight = vbHeight;
		}

		public String toString() {
			return format(US, formatSpecifier, width.value, width.units, height.value, height.units,
			              vbMinX, vbMinY, vbWidth, vbHeight);
		}
	}

	/**
	 * a shape that represents not a particular geographic feature but the outline of the full world map
	 */
	public static class Background implements SVGElement {
		/** the shape under which the background should be added (null means this is an element in an unprojected SVGMap) */
		public final List<Path.Command> shape;
		public final AttributesImpl attributes;

		public Background(AttributesImpl attributes, List<Path.Command> shape) {
			this.shape = shape;
			this.attributes = attributes;
		}

		public String toString() {
			AttributesImpl fullAttributes = new AttributesImpl(attributes);
			fullAttributes.addAttribute("", "", "d", "", Path.toString(shape));
			return formatAttributes("path", fullAttributes);
		}
	}

	/**
	 * an SVG tag containing a single pair of coordinates that can be projected
	 */
	public static class GeographicPoint implements SVGElement {
		public final String formatSpecifier;
		public final double x;
		public final double y;

		public GeographicPoint(String formatSpecifier, double x, double y) {
			this.formatSpecifier = formatSpecifier;
			this.x = x;
			this.y = y;
		}

		public String toString() {
			return format(US, formatSpecifier, x, y);
		}
	}


	/**
	 * a &lt;path d=...&gt; tag, with the d String extracted and parsed into a modifiable form
	 */
	public static class GeographicPath implements SVGElement {
		public final String formatSpecifier;
		public final List<Path.Command> commands;

		public GeographicPath(String formatSpecifier, List<Path.Command> commands) {
			this.formatSpecifier = formatSpecifier;
			this.commands = commands;
		}
		
		public String toString() {
			return format(US, formatSpecifier, Path.toString(commands));
		}
	}
}
