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

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
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

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;
import org.xml.sax.helpers.DefaultHandler;

import utils.Math2;
import utils.SAXUtils;

/**
 * An input equirectangular map based on an SVG file
 * 
 * @author jkunimune
 */
public class SVGMap implements Iterable<SVGMap.Path> {
	
	public static final double[] NULL_TRANSFORM = {1, 1, 0, 0};
	
	private static final double MAX_EDGE_LENGTH = 1/20.; // cut lines that are more than this far across the map
	
	private List<Path> paths; //the set of closed curves in this image
	private List<String> format; //the stuff that goes between the curve descriptions, probably important for something.
	private double vbMinX, vbMinY, vbWidth, vbHeight; //the SVG viewBox
	private double svgWidth, svgHeight; //the actual SVG dimensions
	private int length; //the total number of path commands, for optimization purposes
	
	private static final Map<String, String> ATTRIBUTE_PLACEHOLDERS;
	static {
		Map<String, String> map = new HashMap<String, String>();
		map.put("width", "hmxMLwhWHeqMA8Ba");
		map.put("height", "VlMBunXsmQUtmCw4");
		map.put("viewBox", "UrFo1q9niPDkKSNC"); //attributes of the SVG object to change
		ATTRIBUTE_PLACEHOLDERS = Collections.unmodifiableMap(map);
	}
	
	
	
	public SVGMap(File file) throws IOException, SAXException, ParserConfigurationException {
		paths = new LinkedList<Path>();
		format = new LinkedList<String>();
		
		final SAXParser parser = SAXParserFactory.newInstance().newSAXParser();
		
		final DefaultHandler handler = new DefaultHandler() {
			private Deque<double[]> transformStack =
					new ArrayDeque<double[]>(Collections.singleton(NULL_TRANSFORM));
			private String currentFormatString = "";
			
			@Override
			public InputSource resolveEntity(String publicId, String systemId) {
				return new InputSource(new StringReader("")); //ignore all external references - we don't need to validate
			}
			
			@Override
			public void startElement(
					String uri, String localName, String qName, Attributes attributes) throws SAXException {
				currentFormatString += "<"+qName;
				
				if (attributes.getIndex("transform") >= 0)
					attributes = parseTransform(attributes);
				else
					parseTransform();
				
				if (qName.equals("svg")) {
					attributes = parseViewBox(attributes);
				}
				else if (qName.equals("path")) {
					try {
						attributes = parsePath(attributes);
					} catch (Exception e) {
						throw new SAXException(e.getLocalizedMessage(), null);
					}
				}
				else if (attributes.getIndex("x") >= 0 && attributes.getIndex("y") >= 0) {
					attributes = parsePoint(attributes, false);
				}
				else if (attributes.getIndex("cx") >= 0 && attributes.getIndex("cy") >= 0) {
					attributes = parsePoint(attributes, true);
				}
				
				for (int i = 0; i < attributes.getLength(); i ++)
					currentFormatString +=
							" "+attributes.getQName(i)+"=\""+attributes.getValue(i)+"\"";
				currentFormatString += ">";
			}
			
			@Override
			public void endElement(String uri, String localName, String qName) {
				currentFormatString += "</"+qName+">";
				transformStack.pop();
			}
			
			@Override
			public void characters(char[] ch, int start, int length) {
				for (int i = 0; i < length; i ++) {
					char c = ch[start+i];
//					if (c >= 128 || c == '\'' || c == '"' || c == '<' || c == '>' || c == '&') // some characters must be escaped here
//						currentFormatString += "&#" + (int)c + ";";
//					else
						currentFormatString += c;
				}
			}
			
			@Override
			public void endDocument() {
				format.add(currentFormatString);
			}
			
			private Attributes parseViewBox(Attributes attrs) {
				if (attrs.getValue("width") != null)
					svgWidth = Double.parseDouble(attrs.getValue("width"));
				else
					svgWidth = 360.;
				if (attrs.getValue("height") != null)
					svgHeight = Double.parseDouble(attrs.getValue("height"));
				else
					svgHeight =  180.;
				
				if (attrs.getValue("viewBox") != null) {
					String[] values = attrs.getValue("viewBox").split("\\s", 4);
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
				
				AttributesImpl modAttrs = new AttributesImpl(attrs); //now insert the placeholders for the format string
				for (String qName: ATTRIBUTE_PLACEHOLDERS.keySet()) {
					if (modAttrs.getIndex(qName) >= 0)
						modAttrs.setValue(modAttrs.getIndex(qName), ATTRIBUTE_PLACEHOLDERS.get(qName));
					else
						modAttrs.addAttribute("", "", qName, "", ATTRIBUTE_PLACEHOLDERS.get(qName));
				}
				return modAttrs;
			}
			
			private void parseTransform() {
				transformStack.push(transformStack.peek());
			}
			
			private Attributes parseTransform(Attributes attrs) {
				String transString = attrs.getValue("transform");
				double xScale = 1, yScale = 1, xTrans = 0, yTrans = 0;
				int i;
				while ((i = transString.indexOf('(')) >= 0) {
					int j = transString.indexOf(')');
					String type = transString.substring(0, i).trim();
					String argString = transString.substring(i+1, j);
					String[] args = argString.split("[,\\s]+");
					if (type.equals("matrix")) {
						xScale = Double.parseDouble(args[0]);
						yScale = Double.parseDouble(args[3]);
						xTrans = Double.parseDouble(args[4]);
						yTrans = Double.parseDouble(args[5]);
					}
					else if (type.equals("translate")) {
						if (args.length == 1) {
							xTrans = yTrans = Double.parseDouble(args[0]);
						}
						else {
							xTrans = Double.parseDouble(args[0]);
							yTrans = Double.parseDouble(args[1]);
						}
					}
					else if (type.equals("scale")) {
						if (args.length == 1) {
							xScale = yScale = Double.parseDouble(args[0]);
						}
						else {
							xScale = Double.parseDouble(args[0]);
							yScale = Double.parseDouble(args[1]);
						}
					} //I'm not worrying about shear and rotation because I don't want to
					transString = transString.substring(j+1);
				}
				transformStack.push(new double[] {xScale, yScale, xTrans, yTrans});
				return SAXUtils.removeAttribute(attrs, "transform");
			}
			
			private Attributes parsePath(Attributes attrs) throws Exception {
				currentFormatString += " d=\"";
				format.add(currentFormatString);
				paths.add(new Path(attrs.getValue("d"), transformStack.peek(),
						vbMinX, vbMinY, vbWidth, vbHeight));
				currentFormatString = "\"";
				length += paths.get(paths.size()-1).size();
				return SAXUtils.removeAttribute(attrs, "d");
			}
			
			private Attributes parsePoint(Attributes attrs, boolean center) {
				format.add(currentFormatString); //points are represented as single-point paths
				double[] transform = transformStack.peek();
				double[] coords = new double[2];
				coords[0] = Double.parseDouble(attrs.getValue(center ? "cx" : "x")); //get the coordinates from the attributes
				coords[1] = Double.parseDouble(attrs.getValue(center ? "cy" : "y"));
				coords[0] = transform[0]*coords[0] + transform[2]; //apply any necessary transformation
				coords[1] = transform[1]*coords[1] + transform[3];
				coords[0] = Math2.linInterp(coords[0], vbMinX, vbMinX+vbWidth, -Math.PI, Math.PI);
				coords[1] = Math2.linInterp(coords[1], vbMinY+vbHeight, vbMinY, -Math.PI/2, Math.PI/2);
				paths.add(new Path(new Command(center ? 'O' : 'P', coords)));
				currentFormatString = "";
				length += 1;
				if (center)
					return SAXUtils.removeAttribute(attrs, "cx", "cy");
				else
					return SAXUtils.removeAttribute(attrs, "x", "y");
			}
		};
		
		parser.parse(new BufferedInputStream(new FileInputStream(file)), handler);
	}
	
	
	private SVGMap(List<Path> paths, List<String> format, double vbMinX, double vbMinY,
			double vbWidth, double vbHeight, double svgWidth, double svgHeight, int size) {
		this.paths = paths;
		this.format = format;
		this.vbMinX = vbMinX;
		this.vbMinY = vbMinY;
		this.vbWidth = vbWidth;
		this.vbHeight = vbHeight;
		this.svgWidth = svgWidth;
		this.svgHeight = svgHeight;
		this.length = size;
	}
	
	
	
	public int length() {
		return this.length;
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
		for (String f: this.format)
			newFormat.add(f.replace(target, replacement));
		return new SVGMap(
				paths, newFormat, vbMinX, vbMinY, vbWidth, vbHeight, svgWidth, svgHeight, length);
	}
	
	
	public void save(List<Path> paths, File file, double inMinX, double inMaxY, double inWidth,
			double inHeight) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(file));
		
		final Iterator<String> formatIterator = format.iterator();
		final Iterator<Path> curveIterator = paths.iterator();
		
		out.write(SAXUtils.encode(replacePlaceholders(formatIterator.next(), inWidth/inHeight)));
		while (curveIterator.hasNext()) {
			out.write(closePaths(
					breakWraps(curveIterator.next(), Math.max(inWidth, inHeight))
				).toString(
					inMinX, inMaxY, vbMinX, vbMinY,
					Math.max(vbWidth, vbHeight)/Math.max(inWidth, inHeight)));
			out.write(SAXUtils.encode(formatIterator.next()));
		}
		out.close();
	}
	
	
	private String replacePlaceholders(String str, double desiredAspectRatio) { //change the width, height, and viewBox attributes to match the new aspect ratio
		double outVBWidth, outVBHeight, outSVGWidth, outSVGHeight;
		if (desiredAspectRatio > 1) {
			outVBWidth = Math.max(this.vbWidth, this.vbHeight);
			outVBHeight = outVBWidth/desiredAspectRatio;
			outSVGWidth = Math.max(this.svgWidth, this.svgHeight);
			outSVGHeight = outSVGWidth/desiredAspectRatio;
		}
		else {
			outVBHeight = Math.max(this.vbWidth, this.vbHeight);
			outVBWidth = outVBHeight*desiredAspectRatio;
			outSVGHeight = Math.max(this.svgWidth, this.svgHeight);
			outSVGWidth = outSVGHeight*desiredAspectRatio;
		}
		String viewBoxString = formatDouble(vbMinX) + " " + formatDouble(vbMinY) + " "
				+ formatDouble(outVBWidth) + " " + formatDouble(outVBHeight);
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
				lens[2] = Math.hypot( //compute this next length
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
	
	
	private static String formatDouble(double d) { //format numbers just the way I want them
		String str = String.format("%04d", (int)Math.round(d*1000));
		str = str.substring(0, str.length()-3) + "." + str.substring(str.length()-3);
		return str.replaceFirst("\\.?0*$", "");
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
		
		public Path(String d, double vbWidth, double vbHeight) throws Exception {
			this(d, new double[] {1,1,0,0}, 0, 0, vbWidth, vbHeight);
		}
		
		public Path(String d, double[] transform,
				double vbMinX, double vbMinY, double vbWidth, double vbHeight) throws Exception { //I don't know if this is bad coding practice, but I don't really want to find a way to gracefully catch all possible errors into some more specific Exception class
			super();
			
			int i = 0;
			double[] lastMove = {0, 0}; //for closepaths
			double[] last = {0, 0}; //for relative coordinates
			while (i < d.length()) {
				char type = d.charAt(i);
				String argString = "";
				while (i+1 < d.length() && !isNonELetter(d.charAt(i+1))) {
					i ++;
					argString += d.charAt(i);
				}
				argString = argString.replaceAll("([0-9\\.])-", "$1,-"); //this is necessary because some Adobe products leave out delimiters between negative numbers
				i ++;
				
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
						args[j] = Math2.linInterp(args[j], vbMinX, vbMinX+vbWidth,
								-Math.PI, Math.PI); //scale to radians
					}
					else {
						args[j] = args[j]*transform[1] + transform[3];
						args[j] = Math2.linInterp(args[j], vbMinY+vbHeight, vbMinY, //keep in mind that these are paired longitude-latitude
								-Math.PI/2, Math.PI/2); //not latitude-longitude, as they are elsewhere
					}
				}
				
				this.add(new Command(type, args));
			}
		}
		
		public String toString(
				double inMinX, double inMaxY, double outMinX, double outMinY, double outScale) {
			String s = "";
			for (Command c: this)
				s += c.toString(inMinX, inMaxY, outMinX, outMinY, outScale)+" ";
			return s;
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
		
		public Command(Command command) {
			this(command.type, command.args.clone());
		}
		
		public String toString() {
			return this.toString(-1, -1, 0, 0, 1);
		}
		
		public String toString(
				double inMinX, double inMaxY, double outMinX, double outMinY, double outScale) {
			if (type == 'O') //'O' and 'P' are special; specific points
				return " cx=\""+formatDouble(outMinX + (args[0]-inMinX)*outScale)+
						"\" cy=\""+formatDouble(outMinY + (inMaxY-args[1])*outScale)+"\"";
			else if (type == 'P')
				return " x=\""+formatDouble(outMinX + (args[0]-inMinX)*outScale)+
						"\" y=\""+formatDouble(outMinY + (inMaxY-args[1])*outScale)+"\"";
			
			String s = Character.toString(type);
			for (int i = 0; i < args.length; i ++) {
				if (i%2 == 0)
					s += formatDouble(outMinX + (args[i]-inMinX)*outScale)+",";
				else
					s += formatDouble(outMinY + (inMaxY-args[i])*outScale)+",";
			}
			return s.substring(0, Math.max(1, s.length()-1));
		}
	}
	
}
