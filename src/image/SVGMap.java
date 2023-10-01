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
import utils.SAXUtils;

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
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

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
public class SVGMap implements Iterable<SVGMap.Path> {
	
	public static final double[] NULL_TRANSFORM = {1, 1, 0, 0};
	
	private static final double MAX_EDGE_LENGTH = 1/20.; // cut lines that are more than this far across the map
	
	private final List<Path> paths; //the set of closed curves in this image
	private final List<String> formatStrings; //the stuff that goes between the curve descriptions, like '"></path><path d="'
	private double vbMinX, vbMinY, vbWidth, vbHeight; //the SVG viewBox
	private double svgWidth, svgHeight; //the actual SVG dimensions
	private final int length; //the total number of path commands, for optimization purposes
	
	private static final Map<String, String> ATTRIBUTE_PLACEHOLDERS;
	static {
		Map<String, String> map = new HashMap<String, String>();
		map.put("width", "hmxMLwhWHeqMA8Ba");
		map.put("height", "VlMBunXsmQUtmCw4");
		map.put("viewBox", "UrFo1q9niPDkKSNC"); //attributes of the SVG object to change
		ATTRIBUTE_PLACEHOLDERS = Collections.unmodifiableMap(map);
	}
	
	
	
	public SVGMap(File file) throws IOException, SAXException, ParserConfigurationException {
		paths = new LinkedList<Path>(); // the list of geometric elements
		formatStrings = new LinkedList<String>(); // the list of strings that go between those Paths
		
		final SAXParser parser = SAXParserFactory.newInstance().newSAXParser();
		
		final DefaultHandler handler = new DefaultHandler() {
			private final Deque<double[]> transformStack =
					new ArrayDeque<double[]>(Collections.singleton(NULL_TRANSFORM));
			private StringBuilder currentFormatString = new StringBuilder();

			@Override
			public InputSource resolveEntity(String publicId, String systemId) {
				return new InputSource(new StringReader("")); //ignore all external references - we don't need to validate
			}

			@Override
			public void startElement(
					String uri, String localName, String tagName, Attributes immutableAttributes
			) throws SAXException {
				AttributesImpl attributes = new AttributesImpl(immutableAttributes);

				// first parse any transformation embedded in this element
				if (attributes.getIndex("transform") >= 0) {
					transformStack.push(parseTransform(attributes.getValue("transform")));
					attributes.removeAttribute(attributes.getIndex("transform"));
				}
				else {
					transformStack.push(transformStack.peek());
				}

				// then go thru the <> tag and pull out any important attributes
				currentFormatString.append("<").append(tagName);

				// pull out width, height, and viewBox from <svg>
				if (tagName.equals("svg")) {
					parseSVGHeader(attributes);
					for (String attr_name : ATTRIBUTE_PLACEHOLDERS.keySet()) { //now insert the placeholders for the format string
						if (attributes.getIndex(attr_name) >= 0)
							attributes.setValue(attributes.getIndex(attr_name),
							                    ATTRIBUTE_PLACEHOLDERS.get(attr_name));
						else
							attributes.addAttribute("", "", attr_name, "",
							                        ATTRIBUTE_PLACEHOLDERS.get(attr_name));
					}
					currentFormatString.append(parseAttributes(attributes));
				}
				// pull out the d from <path>
				else if (tagName.equals("path")) {
					currentFormatString.append(" d=\"");
					formatStrings.add(currentFormatString.toString());
					// represent it as a Path object
					try {
						paths.add(parsePath(attributes.getValue("d")));
					} catch (Exception e) {
						throw new SAXException(e.getLocalizedMessage(), null);
					}
					currentFormatString = new StringBuilder("\"");
					currentFormatString.append(parseAttributes(attributes, "d"));
				}
				// pull out the x and y from <rect> or <text>
				else if (attributes.getIndex("x") >= 0 && attributes.getIndex("y") >= 0) {
					formatStrings.add(currentFormatString.toString());
					// represent it as a Path object with one vertex
					paths.add(parsePoint(
							Double.parseDouble(attributes.getValue("x")),
							Double.parseDouble(attributes.getValue("y")),
							'P'));
					currentFormatString = new StringBuilder(
							parseAttributes(attributes, "x", "y"));
				}
				// pull out the cx and cy from <circle>
				else if (attributes.getIndex("cx") >= 0 && attributes.getIndex("cy") >= 0) {
					formatStrings.add(currentFormatString.toString()); //points are represented as single-point paths
					// represent it as a Path object with one vertex
					paths.add(parsePoint(
							Double.parseDouble(attributes.getValue("cx")),
							Double.parseDouble(attributes.getValue("cy")),
							'O'));
					currentFormatString = new StringBuilder(
							parseAttributes(attributes, "cx", "cy"));
				}
				// for everything else, you can just leave all of the attributes in the formatString
				else {
					currentFormatString.append(parseAttributes(attributes));
				}

				currentFormatString.append(">");
			}

			@Override
			public void endElement(String uri, String localName, String qName) {
				// save the </> tag
				currentFormatString.append(format("</%s>", qName));
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
				currentFormatString.append(characterString);
			}

			@Override
			public void endDocument() {
				formatStrings.add(currentFormatString.toString());
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
							xScale = Double.parseDouble(args[0]);
							yScale = Double.parseDouble(args[3]);
							xTrans = Double.parseDouble(args[4]);
							yTrans = Double.parseDouble(args[5]);
							break;
						case "translate":
							if (args.length == 1) {
								xTrans = yTrans = Double.parseDouble(args[0]);
							}
							else {
								xTrans = Double.parseDouble(args[0]);
								yTrans = Double.parseDouble(args[1]);
							}
							break;
						case "scale":
							if (args.length == 1) {
								xScale = yScale = Double.parseDouble(args[0]);
							}
							else {
								xScale = Double.parseDouble(args[0]);
								yScale = Double.parseDouble(args[1]);
							}
							break;
					} // I'm not worrying about shear and rotation because I don't want to
					transString = transString.substring(j + 1);
				}
				return new double[]{xScale, yScale, xTrans, yTrans};
			}

			private void parseSVGHeader(Attributes attributes) {
				if (attributes.getValue("width") != null)
					svgWidth = Double.parseDouble(attributes.getValue("width"));
				else
					svgWidth = 360.;
				if (attributes.getValue("height") != null)
					svgHeight = Double.parseDouble(attributes.getValue("height"));
				else
					svgHeight = 180.;

				if (attributes.getValue("viewBox") != null) {
					String[] values = attributes.getValue("viewBox").split("\\s", 4);
					vbMinX = Double.parseDouble(values[0]);
					vbMinY = Double.parseDouble(values[1]);
					vbWidth = Double.parseDouble(values[2]);
					vbHeight = Double.parseDouble(values[3]);
				}
				else {
					vbWidth = svgWidth;
					vbHeight = svgHeight;
					vbMinX = 0;
					vbMinY = 0;
				}
			}

			private Path parsePath(String pathD) {
				return new Path(pathD, transformStack.peek(),
				                vbMinX, vbMinY, vbWidth, vbHeight);
			}

			private Path parsePoint(double x, double y, char commandLetter) {
				double[] transform = transformStack.peek();
				if (transform == null)
					throw new RuntimeException("the transform stack ran dry");
				x = transform[0]*x + transform[2]; //apply any necessary transformation
				y = transform[1]*y + transform[3];
				x = linInterp(x, vbMinX, vbMinX + vbWidth, -PI, PI);
				y = linInterp(y, vbMinY + vbHeight, vbMinY, -PI/2, PI/2);
				double[] coords = {x, y};
				return new Path(new Command(commandLetter, coords));
			}

			private String parseAttributes(Attributes attributes, String... except) {
				StringBuilder attributeString = new StringBuilder();
				attributeIteration:
				for (int i = 0; i < attributes.getLength(); i++) {
					for (String excludedAttributeName: except)
						if (attributes.getQName(i).equals(excludedAttributeName))
							continue attributeIteration;
					attributeString.append(format(
							" %s=\"%s\"", attributes.getQName(i), attributes.getValue(i)));
				}
				return attributeString.toString();
			}
		};

		parser.parse(new BufferedInputStream(Files.newInputStream(file.toPath())), handler);

		int totalLength = 0;
		for (Path path: paths)
			totalLength += path.size();
		this.length = totalLength;
	}
	
	
	private SVGMap(List<Path> paths, List<String> formatStrings, double vbMinX, double vbMinY,
	               double vbWidth, double vbHeight, double svgWidth, double svgHeight, int length) {
		this.paths = paths;
		this.formatStrings = formatStrings;
		this.vbMinX = vbMinX;
		this.vbMinY = vbMinY;
		this.vbWidth = vbWidth;
		this.vbHeight = vbHeight;
		this.svgWidth = svgWidth;
		this.svgHeight = svgHeight;
		this.length = length;
	}
	
	
	
	public int length() {
		return length;
	}
	
	
	public int numCurves() {
		return paths.size();
	}
	
	
	@Override
	public Iterator<Path> iterator() {
		return paths.iterator();
	}
	
	
	/**
	 * @return the result of replacing all instances of target in the format strings
	 */
	public SVGMap replace(CharSequence target, CharSequence replacement) {
		List<String> newFormat = new LinkedList<String>();
		for (String f: this.formatStrings)
			newFormat.add(f.replace(target, replacement));
		return new SVGMap(
				paths, newFormat, vbMinX, vbMinY, vbWidth, vbHeight, svgWidth, svgHeight, length);
	}
	
	
	public void save(Iterable<Path> paths, File file, double inMinX, double inMaxY, double inWidth,
			double inHeight) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(file));
		
		final Iterator<String> formatIterator = formatStrings.iterator();
		final Iterator<Path> curveIterator = paths.iterator();
		
		out.write(SAXUtils.encode(replacePlaceholders(formatIterator.next(), inWidth/inHeight)));
		while (curveIterator.hasNext()) {
			out.write(closePaths(
					breakWraps(curveIterator.next(), max(inWidth, inHeight))
				).toString(
					inMinX, inMaxY, vbMinX, vbMinY,
					max(vbWidth, vbHeight)/max(inWidth, inHeight)));
			out.write(SAXUtils.encode(formatIterator.next()));
		}
		out.close();
	}
	
	
	private String replacePlaceholders(String str, double desiredAspectRatio) { //change the width, height, and viewBox attributes to match the new aspect ratio
		double outVBWidth, outVBHeight, outSVGWidth, outSVGHeight;
		if (desiredAspectRatio > 1) {
			outVBWidth = max(this.vbWidth, this.vbHeight);
			outVBHeight = outVBWidth/desiredAspectRatio;
			outSVGWidth = max(this.svgWidth, this.svgHeight);
			outSVGHeight = outSVGWidth/desiredAspectRatio;
		}
		else {
			outVBHeight = max(this.vbWidth, this.vbHeight);
			outVBWidth = outVBHeight*desiredAspectRatio;
			outSVGHeight = max(this.svgWidth, this.svgHeight);
			outSVGWidth = outSVGHeight*desiredAspectRatio;
		}
		String viewBoxString = format("%s %s %s %s",
		                              formatDouble(vbMinX), formatDouble(vbMinY),
		                              formatDouble(outVBWidth), formatDouble(outVBHeight));
		str = str.replace(ATTRIBUTE_PLACEHOLDERS.get("viewBox"), viewBoxString);
		str = str.replace(ATTRIBUTE_PLACEHOLDERS.get("width"), formatDouble(outSVGWidth));
		str = str.replace(ATTRIBUTE_PLACEHOLDERS.get("height"), formatDouble(outSVGHeight));
		return str;
	}
	
	
	private Path breakWraps(Path continuous, double inSize) { //break excessively long commands, as they are likely wrapping over a discontinuity
		if (continuous.size() <= 2) 	return continuous;
		Path broken = new Path();
		final double lengthThreshold = inSize*MAX_EDGE_LENGTH;
		double[] lens = {Double.NaN, Double.NaN, Double.NaN}; //the revolving array of command lengths
		for (int i = 0; i < continuous.size(); i ++) {
			if (i < continuous.size()-1 && continuous.get(i+1).type != 'M')
				lens[2] = hypot( //compute this next length
						continuous.get(i+1).args[0] - continuous.get(i).args[0],
						continuous.get(i+1).args[1] - continuous.get(i).args[1]);
			else
				lens[2] = Double.NaN;
			
			char type = continuous.get(i).type;
			if (lens[1] >= lengthThreshold && // check it against an absolute threshold
					(Double.isNaN(lens[0]) || lens[1] > 20*lens[0]) && //and compare it to the last two lengths
					(Double.isNaN(lens[2]) || lens[1] > 20*lens[2])) //if both sides are far longer or nonexistent
				type = 'M'; //break this line
			
			broken.add(new Command(type, continuous.get(i).args.clone()));
			lens[0] = lens[1];
			lens[1] = lens[2];
		}
		
		return broken;
	}
	
	
	private Path closePaths(Path open) { //replace plain loops with 'Z's and combine connected parts
		if (open.size() <= 1) 	return open;
		List<Path> parts = new ArrayList<Path>();
		Path currentPart = null;
		for (Command cmd: open) { //start by breaking the Path into parts,
			if (cmd.type == 'M') { //separated by movetos
				if (currentPart != null)
					parts.add(currentPart);
				currentPart = new Path();
			}
			else if (currentPart == null)
				throw new RuntimeException(format(
						"this path should start with 'M' (and we should have caught that by now): '%s'", open));
			currentPart.add(cmd);
		}
		parts.add(currentPart);
		
		Path closed = new Path();
		for (int i = 0; i < parts.size(); i ++) { //now look through those parts
			Path partI = parts.get(i);
			if (partI.size() > 1
					&& Arrays.equals(partI.get(0).args, partI.get(partI.size()-1).args)) { //if it is self-enclosing
				partI.set(partI.size()-1, new Command('Z', new double[0])); //give it a closepath and send it on its way
			}
			else { //if it is open
				for (int j = i+1; j < parts.size(); j ++) { //look to see if there is anything that completes it
					Path partJ = parts.get(j);
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
		return closed;
	}
	
	
	private static boolean isNonELetter(char c) {
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
	
	
	
	/**
	 * An svg path String, stored in a modifiable form
	 * @author jkunimune
	 */
	public static class Path extends ArrayList<Command> {
		
		private static final long serialVersionUID = 8635895911857570332L;
		
		public Path() {
			super();
		}
		
		public Path(Command... commands) {
			this();
			this.addAll(Arrays.asList(commands));
		}

		public Path(String d, double[] transform,
				double vbMinX, double vbMinY, double vbWidth, double vbHeight) { //I don't know if this is bad coding practice, but I don't really want to find a way to gracefully catch all possible errors into some more specific Exception class
			super();
			
			int i = 0;
			double[] lastMove = {0, 0}; //for closepaths
			double[] last = {0, 0}; //for relative coordinates
			while (i < d.length()) {
				char type = d.charAt(i);
				i ++;
				int start = i;
				while (i < d.length() && !isNonELetter(d.charAt(i)))
					i ++;
				String argString = d.substring(start, i);
				argString = argString.replaceAll("([0-9.])-", "$1,-"); //this is necessary because some Adobe products leave out delimiters between negative numbers

				String[] argStrings;
				if (argString.trim().isEmpty()) 	argStrings = new String[0];
				else 								argStrings = argString.trim().split("[\\s,]+");
				final double[] args;
				
				if (type == 'a' || type == 'A') {
					type += (type == 'a') ? 'l' : 'L'; //change this to a line; I don't want to deal with arcs
					argStrings = new String[] {argStrings[3], argStrings[4]};
				}
				if (type == 'h' || type == 'H' || type == 'v' || type == 'V') { //convert these to 'L'
					final int direcIdx = (type%32 == 8) ? 0 : 1;
					args = new double[] {last[0], last[1]};
					if (type <= 'Z') 	args[direcIdx] = Double.parseDouble(argStrings[0]); //uppercase (absolute)
					else 				args[direcIdx] += Double.parseDouble(argStrings[0]); //lowercase (relative)
					last[direcIdx] = args[direcIdx];
					type = 'L';
				}
				else if (type == 'z' || type == 'Z') { //change this to 'L', too
					args = new double[] {lastMove[0], lastMove[1]};
					type = 'L';
				}
				else {
					args = new double[argStrings.length];
					for (int j = 0; j < args.length; j ++) {
						args[j] = Double.parseDouble(argStrings[j]); //parse the coordinate
						
						if (type >= 'a')
							args[j] += last[j%2]; //account for relative commands
						last[j%2] = args[j];
					}
					if (type >= 'a') //make all letters uppercase
						type -= 32;
				}
				if (type == 'M') { //make note, so we can interpret closepaths properly
					lastMove[0] = args[args.length-2];
					lastMove[1] = args[args.length-1];
				}
				
				for (int j = 0; j < args.length; j ++) {
					if (!Double.isFinite(args[j]))
						throw new IllegalArgumentException("uhh... "+type+argString);
					if (j%2 == 0) {
						args[j] = args[j]*transform[0] + transform[2]; //apply the transformation
						args[j] = linInterp(args[j], vbMinX, vbMinX+vbWidth,
						                    -PI, PI); //scale to radians
					}
					else {
						args[j] = args[j]*transform[1] + transform[3];
						args[j] = linInterp(args[j], vbMinY+vbHeight, vbMinY, //keep in mind that these are paired longitude-latitude
						                     -PI/2, PI/2); //not latitude-longitude, as they are elsewhere
					}
				}
				
				this.add(new Command(type, args));
			}
		}
		
		public String toString(
				double inMinX, double inMaxY, double outMinX, double outMinY, double outScale) {
			StringBuilder s = new StringBuilder();
			for (Command c: this)
				s.append(c.toString(inMinX, inMaxY, outMinX, outMinY, outScale)).append(" ");
			return s.toString();
		}
	}
	
	
	/**
	 * An SVG command, like line or bezier curve or whatever
	 * @author jkunimune
	 */
	public static class Command {
		final public char type; //M, L, C, etc. This will never be lowercase
		final public double[] args; //the absolute coordinates that go with it
		
		public Command(char type, double[] args) {
			this.type = type;
			this.args = args;
		}

		public String toString() {
			return this.toString(-1, -1, 0, 0, 1);
		}
		
		public String toString(
				double inMinX, double inMaxY, double outMinX, double outMinY, double outScale) {
			if (type == 'O') //'O' and 'P' are special; specific points
				return format(" cx=\"%s\" cy=\"%s\"",
				              formatDouble(outMinX + (args[0] - inMinX)*outScale),
				              formatDouble(outMinY + (inMaxY - args[1])*outScale));
			else if (type == 'P')
				return format(" x=\"%s\" cy=\"%s\"",
				              formatDouble(outMinX + (args[0] - inMinX)*outScale),
				              formatDouble(outMinY + (inMaxY - args[1])*outScale));
			else {
				StringBuilder s = new StringBuilder(Character.toString(type));
				for (int i = 0; i < args.length; i ++) {
					if (i%2 == 0)
						s.append(formatDouble(outMinX + (args[i] - inMinX)*outScale)).append(",");
					else
						s.append(formatDouble(outMinY + (inMaxY - args[i])*outScale)).append(",");
				}
				return s.substring(0, max(1, s.length()-1));
			}
		}
	}
	
}
