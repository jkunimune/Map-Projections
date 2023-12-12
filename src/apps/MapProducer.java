/**
 * MIT License
 * 
 * Copyright (c) 2017 Justin KuPermission is hereby granted, free of charge, to any person obtaining a copy
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

import image.SVGMap;
import image.SavableImage;
import javafx.application.Application;
import javafx.concurrent.Task;
import javafx.stage.Stage;
import maps.ArbitraryPseudocylindrical;
import maps.Azimuthal;
import maps.Conic;
import maps.Cylindrical;
import maps.EqualEarth;
import maps.Gyorffy;
import maps.Lenticular;
import maps.Misc;
import maps.Octohedral;
import maps.Polyhedral;
import maps.Projection;
import maps.Pseudocylindrical;
import maps.Snyder;
import maps.Tobler;
import maps.WinkelTripel;

import static java.lang.Math.PI;
import static java.lang.Math.toRadians;

/**
 * A script that automatically converts one input into a bunch of different
 * projections, with multi-threading.
 * 
 * @author Justin Kunimune
 */
public class MapProducer extends Application {
	
	public static final Projection[][] ALL_PROJECTIONS = {
			{ 
					Cylindrical.PLATE_CARREE, Misc.CASSINI,
					Cylindrical.MERCATOR,
					Cylindrical.MERCATOR.withAspect("Transverse Mercator", 0, PI/2, -PI/2),
					Cylindrical.GALL_STEREOGRAPHIC, Cylindrical.MILLER, Cylindrical.LAMBERT,
					Cylindrical.BEHRMANN, Cylindrical.HOBO_DYER, Cylindrical.GALL_ORTHOGRAPHIC,
					Pseudocylindrical.SINUSOIDAL, Pseudocylindrical.MOLLWEIDE,
					Pseudocylindrical.ECKERT_IV, Pseudocylindrical.KAVRAYSKIY_VII,
					ArbitraryPseudocylindrical.ROBINSON, ArbitraryPseudocylindrical.NATURAL_EARTH,
					Tobler.TOBLER, Lenticular.AITOFF, Lenticular.HAMMER, WinkelTripel.WINKEL_TRIPEL,
					Lenticular.VAN_DER_GRINTEN, EqualEarth.EQUAL_EARTH, Conic.EQUIDISTANT,
					Conic.LAMBERT, Conic.ALBERS, Azimuthal.POLAR, Azimuthal.GNOMONIC,
					Azimuthal.EQUAL_AREA, Azimuthal.STEREOGRAPHIC, Azimuthal.ORTHOGRAPHIC,
					Azimuthal.PERSPECTIVE, Misc.TWO_POINT_EQUIDISTANT,
					Misc.PEIRCE_QUINCUNCIAL.transverse("Adams Doubly-Periodic"),
					Polyhedral.DYMAXION, Misc.HAMMER_RETROAZIMUTHAL, Snyder.GS50,
					Azimuthal.STEREOGRAPHIC.withAspect("Oblique Stereographic", .5, 2.5, -2.5),
					Cylindrical.MERCATOR.withAspect("Oblique Mercator", .5, 2.5, 2.5), Misc.BONNE,
					Misc.BRAUN_CONIC, Octohedral.CAHILL_CONCIALDI, Octohedral.KEYES_STANDARD,
					Lenticular.EISENLOHR, Gyorffy.E, Lenticular.LAGRANGE, Lenticular.POLYCONIC,
					Polyhedral.VAN_LEEUWEN, Pseudocylindrical.WAGNER_II, Pseudocylindrical.WAGNER_V,
					Lenticular.WAGNER_VIII, Pseudocylindrical.HOMOLOSINE_INTERRUPTED },
			{
					Misc.PEIRCE_QUINCUNCIAL, Misc.GUYOU, Polyhedral.LEE_TETRAHEDRAL_TRIANGULAR,
					Octohedral.CONFORMAL_CAHILL, Octohedral.WATERMAN } };
	public static final double[] ctrMerids = {0, toRadians(-20)};
	
	
	public static void main(String[] args) {
		launch(args);
	}


	public void start(Stage stage) throws Exception {
		new File("images").mkdirs();
		
		SVGMap[] inputs = { new SVGMap(new File("input/Advanced/Tissot Wikipedia +0.svg")),
				new SVGMap(new File("input/Advanced/Tissot Wikipedia -20.svg")) };
		for (int i = 0; i < 2; i ++) {
			double[] pole = {PI/2, 0, ctrMerids[i]};
			for (Projection proj: ALL_PROJECTIONS[i]) {
				System.out.println(proj);
				
				proj.initialize(proj.getDefaultParameters());
				
				Task<SavableImage> task =
						MapDesignerVector.calculateTask(1, inputs[i], proj, pole, null);
				task.setOnSucceeded((event) -> {
					try {
						task.getValue().save(new File("images/"+proj+".svg"));
					} catch (IOException e) {
						e.printStackTrace();
					}
				});
				task.run();
			}
		}
		
		stop(); // NOTE: I haven't actually set it to terminate when all the tasks finish, so you might have to do that yourself.
	}
}
