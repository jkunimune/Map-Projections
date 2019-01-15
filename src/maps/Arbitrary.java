/**
 * MIT License
 * 
 * Copyright (c) 2018 Justin Kunimune
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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import maps.Projection.Property;
import maps.Projection.Type;

/**
 * A class for completely arbitrary projections, where every square degree can be specified anywhere on the plane.
 * 
 * @author Justin Kunimune
 */
public class Arbitrary {
	
	private static final int X = 0, Y = 1;
	
	
	public static final ArbitraryProjection DANSEIJI_O = new ArbitraryProjection(
			"Danseiji O", "The optimal conventional lenticular map.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiO.csv");
	
	
	public static final ArbitraryProjection DANSEIJI_I = new ArbitraryProjection(
			"Danseiji I", "The optimal conventional equal-area map.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiI.csv");
	
	
	public static final ArbitraryProjection DANSEIJI_II = new ArbitraryProjection(
			"Danseiji II", "Like Danseiji O, but with more weight given to shapes.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiII.csv");
	
	
	public static final ArbitraryProjection DANSEIJI_III = new ArbitraryProjection(
			"Danseiji III", "A map optimised to move distortion from the continents into the oceans.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiIII.csv");
	
	
	public static final ArbitraryProjection DANSEIJI_IV = new ArbitraryProjection(
			"Danseiji IV", "A map optimised to show off the continents by compressing the oceans.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiIV.csv");
	
	
	public static final ArbitraryProjection DANSEIJI_V = new ArbitraryProjection(
			"Danseiji V", "A map optimised to show off the continents by compressing the oceans.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiV.csv");
	
	
	public static final ArbitraryProjection DANSEIJI_VI = new ArbitraryProjection(
			"Danseiji VI", "A map where area is approximately proportional to population.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiVI.csv");
	
	
	public static final ArbitraryProjection DANSEIJI_VII = new ArbitraryProjection(
			"Danseiji VII", "A compromise conventional map, where both physical area and population affect size.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiVII.csv");
	
	
	public static final ArbitraryProjection DANSEIJI_VIII = new ArbitraryProjection(
			"Danseiji VIII", "A compromise unconventional map, where both physical area and population affect size.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiVIII.csv");
	
	
	public static final ArbitraryProjection DANSEIJI_IX = new ArbitraryProjection(
			"Danseiji IX", "A map preserving the oceans over landmasses.",
			true, Type.OTHER, Property.COMPROMISE, "danseijiIV.csv");
	
	
	
	private static class ArbitraryProjection extends Projection {
		
		private String filename; // the data filename
		private double[][][][] cells; // the values of the corner of each cell
		private int[][] cellShapes; // the slope of each cell
		private double[][][] pixels; // the pixel values, for inverse mapping
		private double[][] edge; // the indices of the edge vertices
		
		public ArbitraryProjection(String title, String description, boolean interrupted, Type type, Property property, String filename) {
			super(title, description, 0, 0, interrupted ? 0b1010 : 0b1011, type, property, 4);
			this.filename = filename;
		}
		
		
		public void setParameters(double... params) throws IllegalArgumentException { // these maps don't actually have parameters, but this is the best place to load files
			if (cells != null)
				return; // but don't do it if you've already done it
			
			BufferedReader in = null;
			try {
				in = new BufferedReader(new FileReader(String.format("src/data/%s", filename))); // parsing the input mesh is pretty simple
				String[] row = in.readLine().split(","); // get the header
				double[][] vertices = new double[Integer.parseInt(row[0])][2];
				cells = new double[Integer.parseInt(row[1])][Integer.parseInt(row[2])][][];
				cellShapes = new int[cells.length][cells[0].length];
				edge = new double[Integer.parseInt(row[3])][];
				pixels = new double[Integer.parseInt(row[4])][Integer.parseInt(row[5])][2];
				width = Double.parseDouble(row[6]);
				height = Double.parseDouble(row[7]);
				
				for (int i = 0; i < vertices.length; i ++) { // do the vertex coordinates
					row = in.readLine().split(",");
					for (int j = 0; j < vertices[i].length; j ++)
						vertices[i][j] = Double.parseDouble(row[j]);
				}
				
				for (int i = 0; i < cells.length; i ++) { // get the cell vertices
					for (int j = 0; j < cells[i].length; j ++) {
						row = in.readLine().split(",");
						cellShapes[i][j] = Integer.parseInt(row[0]);
						cells[i][j] = new double[row.length-1][];
						for (int k = 1; k < row.length; k ++)
							cells[i][j][k-1] = vertices[Integer.parseInt(row[k])];
					}
				}
				
				for (int i = 0; i < edge.length; i ++) { // the edge
					row = in.readLine().split(",");
					edge[i] = vertices[Integer.parseInt(row[0])];
				}
				
				for (int i = 0; i < pixels.length; i ++) { // the pixels
					for (int j = 0; j < pixels[i].length; j ++) {
						row = in.readLine().split(",");
						for (int k = 0; k < pixels[i][j].length; k ++)
							pixels[i][j][k] = Double.parseDouble(row[k]);
					}
				}
			} catch (IOException | NullPointerException | ArrayIndexOutOfBoundsException e) {
				cells = new double[][][][] {{{{0,0},{0,0},{0,0},{0,0}}}};
				cellShapes = new int[][] {{0}};
				edge = new double[][] {{0,0}};
				pixels = new double[][][] {{{0,0}}};
				width = 0;
				height = 0;
				e.printStackTrace();
				throw new IllegalArgumentException("Missing or corrupt data file for "+this.getName());
			} finally {
				try {
					if (in != null)	in.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
		
		public double[] project(double lat, double lon) {
			int i = (int)((Math.PI/2-lat)/Math.PI*cells.length); // map it to the array
			i = Math.max(0, Math.min(cells.length-1, i)); // coerce it into bounds
			int j = (int)((lon+Math.PI)/(2*Math.PI)*cells[i].length);
			j = Math.max(0, Math.min(cells[i].length-1, j));
			double cs = (Math.PI/2-lat)/Math.PI*cells.length - i; // do linear interpolation
			double ce = (lon+Math.PI)/(2*Math.PI)*cells[i].length - j;
			double cn = 1 - cs;
			double cw = 1 - ce;
			double[][] v; // which vertices we use depends on some things
			double[] c;
			
			if (cellShapes[i][j] < 0) { // for negative sloped cells,
				if (ce >= cs) { // is it in the ne,
					v = new double[][] {cells[i][j][0], cells[i][j][1], cells[i][j][5]};
					c = new double[] {         1-cw-cs,             cw,             cs};
				}
				else { // or the sw?
					v = new double[][] {cells[i][j][3], cells[i][j][4], cells[i][j][2]};
					c = new double[] {         1-ce-cn,             ce,             cn};
				}
			}
			else if (cellShapes[i][j] > 0) { // for positive sloped cells,
				if (ce >= cn) { // is it in the se,
					v = new double[][] {cells[i][j][5], cells[i][j][0], cells[i][j][4]};
					c = new double[] {         1-cn-cw,             cn,             cw};
				}
				else { // or the nw?
					v = new double[][] {cells[i][j][2], cells[i][j][3], cells[i][j][1]};
					c = new double[] {         1-cs-ce,             cs,             ce};
				}
			}
			else { // for single-element cells
				v = new double[][] {cells[i][j][0], cells[i][j][1], cells[i][j][2], cells[i][j][3]}; // we can just use all of them
				c = new double[] {           cn*ce,          cn*cw,          cs*cw,          cs*ce};
			}
			
			double[] coords = new double[] {0, 0};
			for (int k = 0; k < v.length; k ++) {
				coords[X] += c[k]*v[k][X];
				coords[Y] += c[k]*v[k][Y];
			}
			return coords;
		}
		
		
		public double[] inverse(double x, double y) { // this linear interpolation is much simpler
			boolean inside = false;
			for (int i = 0; i < edge.length; i ++) {
				double x0 = edge[i][0], y0 = edge[i][1]; // for each segment of the edge
				double x1 = edge[(i+1)%edge.length][0], y1 = edge[(i+1)%edge.length][1];
				if ((y0 > y) != (y1 > y)) // if the two points fall on either side of a rightward ray from (X,Y)
					if ((y-y0)/(y1-y0)*(x1-x0)+x0 > x) // and the line between them intersects our ray right of (X,Y)
						inside = !inside; // toggle the boolean
			}
			
			double i = (height/2 - y)/height*(pixels.length-1);
			int i0 = Math.min((int)i, pixels.length-2);
			double cy = i - i0;
			double j = (x + width/2)/width*(pixels[i0].length-1);
			int j0 = Math.min((int)j, pixels[i0].length-2);
			double cx = j - j0;
			
			double X = 0, Y = 0, Z = 0;
			for (int di = 0; di <= 1; di ++) {
				for (int dj = 0; dj <= 1; dj ++) {
					double weight = ((di == 0) ? 1-cy : cy)*((dj == 0) ? 1-cx : cx);
					double phiV = pixels[i0+di][j0+dj][0], lamV = pixels[i0+di][j0+dj][1];
					X += weight*Math.cos(phiV)*Math.cos(lamV);
					Y += weight*Math.cos(phiV)*Math.sin(lamV);
					Z += weight*Math.sin(phiV);
				}
			}
			double phi = Math.atan2(Z, Math.hypot(X, Y)), lam = Math.atan2(Y, X);
			
			if (!inside)	lam += 2*Math.PI; // signal that this point is outside the normal map, if necessary
			return new double[] {phi, lam};
		}
	}
}
