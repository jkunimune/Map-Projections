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
package utils;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
import java.util.function.DoubleConsumer;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.Locator;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * An input equirectangular map based on an SVG file
 * 
 * @author jkunimune
 */
public class SVGMap implements Iterable<SVGMap.Path> {
	
	public double[] NULL_TRANSFORM = {1, 1, 0, 0};
	
	private List<Path> paths; //the set of closed curves in this image
	private List<String> format; //the stuff that goes between the curve descriptions, probably important for something.
	private double minX, minY, width, height; //the SVG viewBox
	private int size; //the total number of path commands, for optimization purposes
	
	
	
	public SVGMap(File file) throws IOException, SAXException, ParserConfigurationException {
		paths = new LinkedList<Path>();
		format = new LinkedList<String>();
		
		final SAXParser parser = SAXParserFactory.newInstance().newSAXParser();
		
		final DefaultHandler handler = new DefaultHandler() {
			
			private Stack<double[]> transformStack = new Stack<double[]>();
			private String currentFormatString = "";
			
			@Override
			public InputSource resolveEntity(String publicId, String systemId) {
				return new InputSource(new StringReader("")); //ignore all external references - we don't need to validate
			}
			
			@Override
			public void notationDecl(String name, String publicId, String systemId) {
				System.out.println("notatien: "+name+", "+publicId+", "+systemId);
			}
			
			@Override
			public void setDocumentLocator(Locator locator) {
				System.out.println("The locator is now "+locator);
			}
			
			@Override
			public void startDocument() {
				System.out.println("begin!");
			}
			
			@Override
			public void unparsedEntityDecl(
					String name, String publicId, String systemId, String notationName) {
				System.out.println("unparsed entity: "+name+", "+publicId+", "+systemId+", "+notationName);
			}
			
			@Override
			public void startElement(
					String uri, String localName, String qName, Attributes attributes) {
				currentFormatString += "<"+qName;
				
				if (attributes.getIndex("transform") >= 0)
					parseTransform(attributes.getValue("transform"));
				else
					parseTransform();
				
				if (qName.equals("path"))
					parsePath(attributes.getValue("d"));
				if (qName.equals("svg"))
					parseViewBox(attributes.getValue("viewBox"));
				
				for (int i = 0; i < attributes.getLength(); i ++)
					if (!attributes.getQName(i).equals("d") && //d is already taken care of
							!attributes.getQName(i).equals("transform")) //there shall be no transforms in the final output
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
				for (int i = 0; i < length; i ++)
					currentFormatString += ch[start+i];
			}
			
			@Override
			public void endDocument() {
				format.add(currentFormatString);
			}
			
			private void parseViewBox(String viewBox) {
				if (viewBox != null) {
					String[] values = viewBox.split("\\s", 4);
					minX = Double.parseDouble(values[0]);
					minY = Double.parseDouble(values[1]);
					width = Double.parseDouble(values[2]);
					height = Double.parseDouble(values[3]);
				}
			}
			
			private void parseTransform() {
				if (transformStack.isEmpty())
					transformStack.push(NULL_TRANSFORM);
				else
					transformStack.push(transformStack.peek());
			}
			
			private void parseTransform(String transString) {
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
			}
			
			private void parsePath(String d) {
				currentFormatString += " d=\"";
				format.add(currentFormatString);
				paths.add(new Path(d, transformStack.peek(), minX, minY, width, height));
				currentFormatString = "\"";
				size += paths.get(paths.size()-1).length();
			}
		};
		
		parser.parse(new BufferedInputStream(new FileInputStream(file)), handler);
	}
	
	
	
	public int size() {
		return this.size;
	}
	
	
	public int numCurves() {
		return paths.size();
	}
	
	
	@Override
	public Iterator<Path> iterator() {
		return paths.iterator();
	}
	
	
	public void save(List<Path> paths, File file, DoubleConsumer tracker) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(file));
		
		int i = 0;
		final Iterator<String> formatIterator = format.iterator();
		final Iterator<Path> curveIterator = paths.iterator();
		while (curveIterator.hasNext()) {
			out.write(formatIterator.next());
			out.write(curveIterator.next().toString(
					minX, minY, Math.min(width, height), Math.min(width, height)));
			tracker.accept((double)i/paths.size());
			i ++;
		}
		out.write(formatIterator.next());
		out.close();
	}
	
	
	public static boolean isLetter(char c) {
		return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
	}
	
	
	
	/**
	 * An svg path String, stored in a modifiable form
	 * @author jkunimune
	 */
	public static class Path implements Iterable<Command> {
		
		final private List<Command> commands;
		
		public Path() {
			commands = new ArrayList<Command>();
		}
		
		public Path(String d, double vbWidth, double vbHeight) {
			this(d, new double[] {1,1,0,0}, 0, 0, vbWidth, vbHeight);
		}
		
		public Path(String d, double[] transform,
				double vbMinX, double vbMinY, double vbWidth, double vbHeight) {
			commands = new ArrayList<Command>();
			int i = 0;
			double[] last = new double[] {0,0}; //for relative coordinates
			while (i < d.length()) {
				char type = d.charAt(i);
				String argString = "";
				while (i+1 < d.length() && !isLetter(d.charAt(i+1))) {
					i ++;
					argString += d.charAt(i);
				}
				i ++;
				
				String[] argStrings;
				if (argString.trim().isEmpty()) 	argStrings = new String[0];
				else 								argStrings = argString.trim().split("[\\s,]+");
				final double[] args;
				
				if (type == 'a' || type == 'A') {
					type += 11;
					argStrings = new String[] {argStrings[3], argStrings[4]};
				}
				if (type == 'h' || type == 'H' || type == 'v' || type == 'V') { //convert these to 'L'
					final int chgIdx = (type%32 == 8) ? 0 : 1;
					args = new double[] {last[0], last[1]};
					if (type <= 'Z') 	args[chgIdx] = Double.parseDouble(argStrings[0]); //uppercase (absolute)
					else 				args[chgIdx] += Double.parseDouble(argStrings[0]); //lowercase (relative)
					last[chgIdx] = args[chgIdx];
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
				
				for (int j = 0; j < args.length; j ++) {
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
				
				commands.add(new Command(type, args));
			}
		}
		
		public Iterator<Command> iterator() {
			return commands.iterator();
		}
		
		public int length() {
			return commands.size();
		}
		
		public boolean add(Command c) {
			return commands.add(c);
		}
		
		public String toString() {
			return this.toString(-1, -1, 2, 2);
		}
		
		public String toString(double minX, double minY, double width, double height) {
			String s = "";
			for (Command c: commands)
				s += c.toString(minX, minY, width, height)+" ";
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
		
		public String toString() {
			return this.toString(-1, -1, 2, 2);
		}
		
		public String toString(double minX, double minY, double width, double height) {
			String s = type+" ";
			for (int i = 0; i < args.length; i ++) {
				if (i%2 == 0)
					s += Math2.linInterp(args[0], -1, 1, minX, minX+width)+",";
				else
					s += Math2.linInterp(args[1], -1, 1, minY+height, minY)+",";
			}
			return s.substring(0, s.length()-1);
		}
	}
}
