package apps;

import java.io.File;
import java.io.IOException;

import image.SVGMap;
import image.SavableImage;
import javafx.application.Application;
import javafx.concurrent.Task;
import javafx.stage.Stage;
import maps.Arbitrary;
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

/**
 * A script that automatically manufactures a slew of projections in standard aspect.
 * 
 * @author jkunimune
 */
public class MapProducer extends Application {
	
	public static final Projection[][] ALL_PROJECTIONS = {
			{ Cylindrical.PLATE_CARREE,
					Cylindrical.PLATE_CARREE.withAspect("Cassini", 0,Math.PI/2, -Math.PI/2),
					Cylindrical.MERCATOR,
					Cylindrical.MERCATOR.withAspect("Transverse Mercator", 0,Math.PI/2,-Math.PI/2),
					Cylindrical.GALL_STEREOGRAPHIC, Cylindrical.MILLER, Cylindrical.LAMBERT,
					Cylindrical.BEHRMANN, Cylindrical.HOBO_DYER, Cylindrical.GALL_ORTHOGRAPHIC,
					Pseudocylindrical.SINUSOIDAL, Pseudocylindrical.MOLLWEIDE,
					Pseudocylindrical.ECKERT_IV, Pseudocylindrical.KAVRAYSKIY_VII,
					Arbitrary.ROBINSON, Arbitrary.NATURAL_EARTH, Tobler.TOBLER, Lenticular.AITOFF,
					Lenticular.HAMMER, WinkelTripel.WINKEL_TRIPEL, Lenticular.VAN_DER_GRINTEN,
					Conic.EQUIDISTANT, Conic.LAMBERT, Conic.ALBERS, Azimuthal.POLAR,
					Azimuthal.GNOMONIC, Azimuthal.EQUAL_AREA, Azimuthal.STEREOGRAPHIC,
					Azimuthal.ORTHOGRAPHIC, Azimuthal.PERSPECTIVE, Misc.TWO_POINT_EQUIDISTANT,
					Misc.PEIRCE_QUINCUNCIAL.transverse("Adams Doubly-Periodic"),
					Polyhedral.DYMAXION, Misc.HAMMER_RETROAZIMUTHAL, Snyder.GS50,
					Azimuthal.STEREOGRAPHIC.withAspect("Oblique Stereographic", .5,2.5,-2.5),
					Cylindrical.MERCATOR.withAspect("Oblique Mercator", .5,2.5,2.5) },
			{ Misc.PEIRCE_QUINCUNCIAL, Misc.GUYOU, Polyhedral.LEE_TETRAHEDRAL_TRIANGULAR,
					Octohedral.CAHILL_KEYES, Octohedral.WATERMAN } };
	
	
	public static void main(String[] args) {
		launch(args);
	}


	public void start(Stage stage) throws Exception {
		new File("images").mkdirs();
		
		SVGMap[] inputs = { new SVGMap(new File("input/Advanced/Tissot Wikipedia +0.svg")),
				new SVGMap(new File("input/Advanced/Tissot Wikipedia -20.svg")) };
		double[] ctrMerids = {0, Math.toRadians(-20)};
		for (int i = 0; i < 2; i ++) {
			double[] pole = {Math.PI/2, 0, ctrMerids[i]};
			for (Projection proj: ALL_PROJECTIONS[i]) {
				System.out.println(proj);
				
				if (proj == Tobler.TOBLER)
					proj.setParameters(37.07, 0., 3.);
				else if (proj == WinkelTripel.WINKEL_TRIPEL)
					proj.setParameters(50.46);
				else
					proj.setParameters(proj.getDefaultParameters());
				
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
		
		stage.close();
	}
}
