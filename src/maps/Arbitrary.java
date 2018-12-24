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
import utils.Math2;

/**
 * A class for completely arbitrary projections, where every square degree can be specified anywhere on the plane.
 * 
 * @author Justin Kunimune
 */
public class Arbitrary {
	
	private static final int NE = 0, NW = 1, SW = 2, SE = 3;
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
	
	
	
	private static class ArbitraryProjection extends Projection {
		
		private double[][] vertices; // the vertex x-y coordinates
		private int[][][] cells; // the indices of the corner of each cell
		
		public ArbitraryProjection(String title, String description, boolean interrupted, Type type, Property property, String filename) {
			super(title, description, 0, 0, interrupted ? 0b1010 : 0b1011, type, property, 4);
			
			BufferedReader in = null;
			try {
				in = new BufferedReader(new FileReader(String.format("src/data/%s", filename))); // parsing the input mesh is pretty simple
				String[] row = in.readLine().split(","); // get the header
				vertices = new double[Integer.parseInt(row[0])][2];
				cells = new int[Integer.parseInt(row[1])][Integer.parseInt(row[2])][4];
				width = Double.parseDouble(row[3]);
				height = Double.parseDouble(row[4]);
				for (int i = 0; i < vertices.length; i ++) { // do the vertex coordinates
					row = in.readLine().split(",");
					for (int j = 0; j < vertices[i].length; j ++)
						vertices[i][j] = Double.parseDouble(row[j]);
				}
				for (int i = 0; i < cells.length; i ++) { // get the cell vertices
					for (int j = 0; j < cells[i].length; j ++) {
						row = in.readLine().split(",");
						for (int k = 0; k < cells[i][j].length; k ++)
							cells[i][j][k] = Integer.parseInt(row[k]);
					}
				} // and skip the edge; it's not relevant here
			} catch (IOException | NullPointerException | ArrayIndexOutOfBoundsException e) {
				System.err.println("Could not load mesh: "+e);
				width = 2;
				height = 2;
				vertices = new double[][] {{1,1},{-1,1},{-1,-1},{1,-1}};
				cells = new int[][][] {{{0,1,2,3}}};
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
			double[] ne = vertices[cells[i][j][NE]], nw = vertices[cells[i][j][NW]],
					sw = vertices[cells[i][j][SW]], se = vertices[cells[i][j][SE]];
			return new double[] {
					cs*cw*sw[X] + cs*ce*se[X] + cn*cw*nw[X] + cn*ce*ne[X],
					cs*cw*sw[Y] + cs*ce*se[Y] + cn*cw*nw[Y] + cn*ce*ne[Y] };
		}
		
		
		public double[] inverse(double x, double y) {
			return WinkelTripel.WINKEL_TRIPEL.inverse(
					x*WinkelTripel.WINKEL_TRIPEL.getWidth()/this.getWidth(),
					y*WinkelTripel.WINKEL_TRIPEL.getHeight()/this.getHeight(), Projection.NORTH_POLE, 40);
		}
	}
}
