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
import static utils.Math2.linInterp;
import static utils.Math2.max;

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
			public void characters(char[] ch, int start, int length) {
				// just dump the provided characters into currentFormatString
				StringBuilder characterString = new StringBuilder();
				for (int i = 0; i < length; i++) {
					characterString.append(ch[start + i]);
				}
				elements.add(new Content(characterString.toString()));
			}
			
			@Override
			public void declaration(String version, String encoding, String standalone) {
				elements.add(new Content(String.format(
						"<?xml version=\"%s\" encoding=\"%s\" standalone=\"%s\"?>\n",
						version, encoding, standalone)));
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
				double[] lastMove = {0, 0}; //for closepaths
				double[] last = {0, 0}; //for relative coordinates
				for (int i = 0; i < path.size(); i ++) {
					char type = path.get(i).type;
					double[] args = path.get(i).args;
					// replace arcs with lines because I don't want to deal with those
					if (type == 'a' || type == 'A') {
						type = (type == 'a') ? 'l' : 'L';
						args = new double[] {args[5], args[6]};
					}
					// convert horizontal and vertical lines to general lines to keep things simple
					if (type == 'h' || type == 'H' || type == 'v' || type == 'V') {
						final int direcIdx = (type%32 == 8) ? 0 : 1;
						double[] newArgs = new double[] {last[0], last[1]};
						if (type <= 'Z')
							newArgs[direcIdx] = args[0]; //uppercase (absolute)
						else
							newArgs[direcIdx] += args[0]; //lowercase (relative)
						args = newArgs;
						last[direcIdx] = args[direcIdx];
						type = 'L';
					}
					// also replace closepath with an explicit line to the last moveto
					else if (type == 'z' || type == 'Z') {
						args = new double[] {lastMove[0], lastMove[1]};
						type = 'L';
					}
					else {
						for (int j = 0; j < args.length; j ++) {
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
	public double getDisplayWidth() {
		return this.header.width;
	}

	/**
	 * @return the SVG's display height
	 */
	public double getDisplayHeight() {
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
						parts.remove(j); //don't look at J anymone; it has been absorbed.
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

		public GeographicPoint projected() {
			return this;
		}

		public String toString() {
			return format(formatSpecifier, x, y);
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
			return format(formatSpecifier, Path.toString(commands));
		}
	}
}
