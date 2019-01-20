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

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import image.PixelMap;
import image.SavableImage;
import javafx.application.Application;
import javafx.concurrent.Task;
import javafx.stage.Stage;
import maps.Arbitrary;
import maps.ArbitraryPseudocylindrical;
import maps.Azimuthal;
import maps.Conic;
import maps.Cylindrical;
import maps.Lenticular;
import maps.Misc;
import maps.Octohedral;
import maps.Polyhedral;
import maps.Projection;
import maps.Pseudocylindrical;
import maps.Snyder;
import maps.Tobler;
import maps.WinkelTripel;

/**
 * A program to generate and output an HTML snippet listing and explaining all of my maps
 * 
 * @author Justin Kunimune
 */
public class MapExplainer extends Application {
	
	public static final int IMG_WIDTH = 800;
	public static final int SMOOTHING = 3;
	
	public static final Projection[][] ALL_PROJECTIONS = {
			{ Cylindrical.MERCATOR, Cylindrical.PLATE_CARREE, Cylindrical.GALL_ORTHOGRAPHIC,
					Cylindrical.HOBO_DYER, Cylindrical.BEHRMANN, Cylindrical.LAMBERT,
					Cylindrical.GALL_STEREOGRAPHIC, Azimuthal.STEREOGRAPHIC.transverse(),
					Azimuthal.POLAR.transverse(), Azimuthal.EQUAL_AREA.transverse(),
					Azimuthal.GNOMONIC.transverse(), Azimuthal.ORTHOGRAPHIC.transverse(),
					Conic.LAMBERT, Conic.EQUIDISTANT, Conic.ALBERS, Pseudocylindrical.SINUSOIDAL,
					Pseudocylindrical.MOLLWEIDE, Tobler.TOBLER, Lenticular.AITOFF,
					Lenticular.VAN_DER_GRINTEN, ArbitraryPseudocylindrical.ROBINSON,
					WinkelTripel.WINKEL_TRIPEL, Polyhedral.AUTHAGRAPH,
					Polyhedral.LEE_TETRAHEDRAL_RECTANGULAR, Octohedral.CAHILL_KEYES,
					Octohedral.CAHILL_CONCIALDI, Misc.PEIRCE_QUINCUNCIAL.transverse(), Snyder.GS50,
					Misc.TWO_POINT_EQUIDISTANT, Misc.HAMMER_RETROAZIMUTHAL, Misc.FLAT_EARTH },
			{ Arbitrary.DANSEIJI_O, Arbitrary.DANSEIJI_IV, Arbitrary.DANSEIJI_V, Polyhedral.TETRAGRAPH,
					Polyhedral.AUTHAPOWER, Polyhedral.ACTUAUTHAGRAPH } };
	
	private PixelMap inputSkew, inputPole, inputNone, inputEast, inputWest;
	
	
	public static void main(String[] args) {
		launch(args);
	}
	
	
	public void start(Stage stage) throws Exception {
		new File("images").mkdirs();
		try {
			inputSkew = new PixelMap(new File("input/Advanced/Tissot Oblique.jpg"));
			inputPole = new PixelMap(new File("input/Advanced/Tissot Standard.jpg"));
			inputEast = new PixelMap(new File("input/Advanced/Tissot Shift East.jpg"));
			inputWest = new PixelMap(new File("input/Advanced/Tissot Shift West.jpg"));
			inputNone = new PixelMap(new File("input/Advanced/Tissot Not.jpg"));
		} catch(IOException e) {
			e.printStackTrace();
			stop();
			return;
		}
		
		final PrintStream out = System.out;
		for (Projection[] projs: ALL_PROJECTIONS) {
			out.println("<h1>Map Projections</h1>");
			
			for (Projection proj: projs) {
				proj.setParameters(proj.getDefaultParameters());
				out.println("<h2>"+proj.getName()+"</h2>");
				
				PixelMap input = inputSkew;
				if (!proj.hasAspect() || proj == Polyhedral.AUTHAGRAPH)
					input = inputPole;
				if (proj == Octohedral.CAHILL_KEYES) //some projections only look good in standard aspect
					input = inputWest;
				if (proj == Octohedral.CAHILL_CONCIALDI)
					input = inputEast;
				if (proj == Misc.FLAT_EARTH) //or can't be seen with indicatrices
					input = inputNone;
				
				Task<SavableImage> task = MapDesignerRaster.calculateTask(
						IMG_WIDTH, (int)(IMG_WIDTH/proj.getAspectRatio()), 2, input, proj, null,
						false, 0, null);
				task.setOnSucceeded((event) -> {
					try {
						task.getValue().save(new File("images/"+proj+".gif"));
					} catch (IOException e) {
						e.printStackTrace();
					}
				});
				task.run();
				
				out.println("<div class=\"row\">");
				out.println("  <div class=\"col-left\">");
				out.println("    <p>"+proj.getDescription()+"</p>");
				
				out.println("    <dl>");
				out.println("      <dt>Rating:&nbsp;</dt>");
				out.print("      <dd>");
				for (int i = 0; i < proj.getRating(); i ++)
					out.print("&#9733;");
				for (int i = 0; i < 4-proj.getRating(); i ++)
					out.print("&#9734;");
				out.println("</dt>");
				out.println("      <dt>Geometry:&nbsp;</dt>");
				out.println("      <dd>"+proj.getType().getName()+"</dd>");
				out.println("      <dt>Property:&nbsp;</dt>");
				out.println("      <dd>"+proj.getProperty().getName()+"</dd>");
				out.println("      <dt>Uninterrupted:&nbsp;</dt>");
				out.println("      <dd>"+ (proj.isContinuous() ? "Yes":"No") +"</dd>");
				out.println("      <dt>Shows entire world:&nbsp;</dt>");
				out.println("      <dd>"+ (proj.isFinite() ? "Yes":"No") +"</dd>");
				out.println("      <dt>Closed-form solution:&nbsp;</dt>");
				out.println("      <dd>"+ (proj.isSolveable() ? "Yes":"No") +"</dd>");
				out.println("      <dt>Closed-form inverse:&nbsp;</dt>");
				out.println("      <dd>"+ (proj.isInvertable() ? "Yes":"No") +"</dd>");
				out.println("    </dl>");
				out.println("    <br>");
				out.println("  </div>");
				
				out.println("  <div class=\"col-right\">");
				out.println("    <img src=\"images/"+proj+".gif\" class=\"map\" " +
						"alt=\"Sorry. This map is not available at the moment. " +
						"Please leave a message after the beep.\" " +
						"title=\"TODO\">");
				out.println("  </div>");
				out.println("</div>");
				out.println("<br>");
				
				out.println();
			}
		}
		
		stop();
	}
}
