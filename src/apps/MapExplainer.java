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
package apps;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import javax.imageio.ImageIO;

import maps.Azimuthal;
import maps.Conic;
import maps.Cylindrical;
import maps.Misc;
import maps.MyProjections;
import maps.Projection;
import maps.Pseudocylindrical;
import maps.Robinson;
import maps.Tetrahedral;
import maps.Tobler;
import maps.WinkelTripel;
import utils.ImageUtils;

/**
 * A program to generate and output an HTML snippet listing and explaining all of my maps
 * 
 * @author jkunimune
 */
public class MapExplainer {
	
	public static final int IMG_WIDTH = 800;
	public static final int SMOOTHING = 3;
	
	public static final Projection[][] ALL_PROJECTIONS = {
			{ Cylindrical.MERCATOR, Cylindrical.PLATE_CARREE, Cylindrical.GALL_PETERS,
					Cylindrical.HOBO_DYER, Cylindrical.BEHRMANN, Cylindrical.LAMBERT,
					Cylindrical.GALL, Azimuthal.STEREOGRAPHIC.transverse(),
					Azimuthal.POLAR.transverse(), Azimuthal.EQUAL_AREA.transverse(),
					Azimuthal.GNOMONIC.transverse(), Azimuthal.ORTHOGRAPHIC.transverse(),
					Conic.LAMBERT, Conic.EQUIDISTANT, Conic.ALBERS, Tetrahedral.LEE,
					Tetrahedral.TETRAGRAPH, Tetrahedral.AUTHAGRAPH, Pseudocylindrical.SINUSOIDAL,
					Pseudocylindrical.MOLLWEIDE, Misc.HAMMER, Tobler.TOBLER, Misc.AITOFF,
					Misc.VAN_DER_GRINTEN, Robinson.ROBINSON, WinkelTripel.WINKEL_TRIPEL,
					Misc.PEIRCE_QUINCUNCIAL.transverse(), Misc.GUYOU, Misc.TWO_POINT_EQUIDISTANT,
					Misc.HAMMER_RETROAZIMUTHAL },
			
			{ MyProjections.PSEUDOSTEREOGRAPHIC, MyProjections.HYPERELLIPOWER,
					Tetrahedral.TETRAPOWER, Tetrahedral.TETRAFILLET,
					MyProjections.TWO_POINT_EQUALIZED.transverse() } };
	
	
	public static void main(String[] args) throws IOException {
		
		new File("images").mkdirs();
		final BufferedImage inputSkew = ImageIO.read(new File("input/Tissot-alt2.jpg"));
		final BufferedImage inputPole = ImageIO.read(new File("input/Tissot-alt1.jpg"));
		final PrintStream out = System.out;
		
		for (Projection[] projs: ALL_PROJECTIONS) {
			out.println("<h2>Map Projections</h2>");
			
			for (Projection proj: projs) {
				out.println("<div class=row>");
				out.println("  <div class=\"col-xs-12 col-sm-10 col-md-8 col-lg-6\">");
				out.println("    <h3>"+proj.getName()+"</h3>");
				out.println("  </div>");
				out.println("</div>");
				
				ImageIO.write(
						makeImage(proj, proj.hasAspect() ? inputSkew : inputPole),
						"gif", new File("images/"+proj+".gif"));
				
				out.println("<div class=row>");
				out.println("  <div class=\"col-xs-12 col-sm-10 col-md-4 col-lg-3\">");
				out.println("    <p>"+proj.getDescription()+"</p>");
				
				out.println("    <dl>");
				out.println("      <dt>Geometry:&nbsp;</dt>");
				out.println("      <dd>"+proj.getType().getName()+"</dd>");
				out.println("      <dt>Property:&nbsp;</dt>");
				out.println("      <dd>"+proj.getProperty().getName()+"</dd>");
				out.println("      <dt>Continuous:&nbsp;</dt>");
				out.println("      <dd>"+ (proj.isContinuous() ? "Yes":"No") +"</dd>");
				out.println("      <dt>Shows entire world:&nbsp;</dt>");
				out.println("      <dd>"+ (proj.isFinite() ? "Yes":"No") +"</dd>");
				out.println("      <dt>Closed-form solution:&nbsp;</dt>");
				out.println("      <dd>"+ (proj.isSolveable() ? "Yes":"No") +"</dd>");
				out.println("      <dt>Closed-form inverse:&nbsp;</dt>");
				out.println("      <dd>"+ (proj.isInvertable() ? "Yes":"No") +"</dd>");
				out.println("    </dl>");
				out.println("  </div>");
				
				out.println("  <div class=\"col-xs-12 col-sm-10 col-md-4 col-lg-3\">");
				out.println("    <img src=\"images/"+proj+".gif\" class=\"map\" " +
						"alt=\"Sorry. This map is not available at the moment. " +
						"Please leave a message after the beep.\" " +
						"title=\"I'll think of something to put here later.\">");
				out.println("  </div>");
				out.println("</div>");
			}
		}
		
	}
	
	
	private static BufferedImage makeImage(Projection proj, BufferedImage input) {
		proj.setParameters(proj.getDefaultParameters());
		final BufferedImage out = new BufferedImage(
				IMG_WIDTH, (int) (IMG_WIDTH/proj.getAspectRatio()), BufferedImage.TYPE_INT_ARGB);
		for (int y = 0; y < out.getHeight(); y ++) {
			for (int x = 0; x < out.getWidth(); x ++) {
				int[] colors = new int[SMOOTHING*SMOOTHING];
				for (int dy = 0; dy < SMOOTHING; dy ++) {
					for (int dx = 0; dx < SMOOTHING; dx ++) {
						final double X = 2*(x+(dx+.5)/SMOOTHING)/out.getWidth()-1;
						final double Y = 1-2*(y+(dy+.5)/SMOOTHING)/out.getHeight();
						colors[SMOOTHING*dy+dx] = ImageUtils.getArgb(
								proj.inverse(X, Y), input);
					}
				}
				out.setRGB(x, y, ImageUtils.blend(colors));
			}
		}
		return out;
	}
}
