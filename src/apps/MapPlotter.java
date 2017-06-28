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

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.chart.XYChart.Series;
import javafx.stage.Stage;
import maps.Projection;
import util.Math2;

/**
 * A simple script that creates an annotated ScatterPlot of map projections
 * 
 * @author jkunimune
 */
public class MapPlotter extends Application {

	
	private static final Projection[] COMMON = { Projection.STEREOGRAPHIC, Projection.LAMBERT_CONIC,
			Projection.MERCATOR, Projection.POLAR, Projection.E_D_CONIC, Projection.EQUIRECTANGULAR,
			Projection.E_A_AZIMUTH, Projection.ALBERS, Projection.HOBODYER, Projection.VAN_DER_GRINTEN,
			Projection.ROBINSON, Projection.WINKEL_TRIPEL, Projection.GNOMONIC, Projection.ORTHOGRAPHIC };
	private static final Projection[] UNCOMMON = { Projection.PEIRCE_QUINCUNCIAL, Projection.LEE, Projection.TOBLER,
			Projection.MOLLWEIDE, Projection.AITOFF };
	private static final Projection[] INVENTED = { Projection.HYPERELLIPOWER, Projection.TETRAPOWER,
			Projection.TETRAFILLET };
	
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage root) {
		final ScatterChart<Number, Number> plot =
				new ScatterChart<Number, Number>(
				new NumberAxis("Size distortion", 0, 3, 0.2),
				new NumberAxis("Shape distortion", 0, 3, 0.2));
		final double[][][] points = Projection.globe(.02);
		
		plotProjections(plot, COMMON, "Common", points);
		plotProjections(plot, UNCOMMON, "Hipster", points);
		plotProjections(plot, INVENTED, "Mine", points);
		
		root.setTitle("Map Projections");
		root.setScene(new Scene(plot));
		root.show();
	}
	
	
	private static void plotProjections(ScatterChart<Number, Number> chart,
			Projection[] projections, String name, double[][][] points) {
		final Series<Number, Number> series = new Series<Number, Number>();
		series.setName(name);
		
		for (Projection projection: projections) {
			final double distortion[][][] = projection.calculateDistortion(points);
			final double sizeDist = Math2.stdDev(distortion[0]);
			final double anglDist = Math2.mean(distortion[1]);
			series.getData().add(new Data<Number, Number>(sizeDist, anglDist));
		}
		
		chart.getData().add(series);
	}

}
