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
import java.util.LinkedList;
import java.util.List;

import javax.imageio.ImageIO;

import javafx.application.Application;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.control.Label;
import javafx.scene.layout.AnchorPane;
import javafx.scene.layout.StackPane;
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
			Projection.MERCATOR, Projection.POLAR, Projection.E_D_CONIC, Projection.PLATE_CARREE,
			Projection.E_A_AZIMUTH, Projection.ALBERS, Projection.HOBO_DYER, Projection.VAN_DER_GRINTEN,
			Projection.ROBINSON, Projection.WINKEL_TRIPEL };
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
		final AnchorPane overlay = new AnchorPane();
		final StackPane stack = new StackPane(plot, overlay);
		final List<Label> labels = new LinkedList<Label>();
		final List<double[]> coords = new LinkedList<double[]>();
		final double[][][] points = Projection.globe(.02);
		
		plotProjections(plot, overlay, labels, coords, COMMON, "Common", points);
		plotProjections(plot, overlay, labels, coords, UNCOMMON, "Hipster", points);
		plotProjections(plot, overlay, labels, coords, INVENTED, "Mine", points);
		
		root.setTitle("Map Projections");
		root.setScene(new Scene(stack));
		root.show();
		
		addLabels(plot, labels, coords);
		
		try {
			ImageIO.write(
					SwingFXUtils.fromFXImage(stack.snapshot(null, null), null),
					"png", new File("output/graph - plotter.png"));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	private static void plotProjections(ScatterChart<Number, Number> chart,
			AnchorPane overlay, List<Label> labels, List<double[]> coords,
			Projection[] projections, String name, double[][][] points) {
		final Series<Number, Number> series = new Series<Number, Number>();
		series.setName(name);
		
		for (Projection projection: projections) {
			final double[] params = projection.getDefaultParameters();
			final double distortion[][][] = projection.calculateDistortion(points, params);
			final double sizeDist = Math2.stdDev(distortion[0]);
			final double anglDist = Math2.mean(distortion[1]);
			series.getData().add(new Data<Number, Number>(sizeDist, anglDist));
			final Label lbl = new Label(projection.getName());
			overlay.getChildren().add(lbl);
			labels.add(lbl);
			coords.add(new double[] {sizeDist, anglDist});
		}
		
		chart.getData().add(series);
	}
	
	
	private static void addLabels(ScatterChart<Number,Number> chart, List<Label> labels, List<double[]> coords) {
		for (int i = 0; i < labels.size(); i ++) {
			AnchorPane.setLeftAnchor(labels.get(i),
					chart.getXAxis().localToParent(chart.getXAxis().getDisplayPosition(coords.get(i)[0]), 0).getX() + chart.getPadding().getLeft());
			AnchorPane.setBottomAnchor(labels.get(i),
					chart.getHeight()-chart.getYAxis().localToParent(0, chart.getYAxis().getDisplayPosition(coords.get(i)[1])).getY() + chart.getPadding().getBottom());
		}
	}

}
