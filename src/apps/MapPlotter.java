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
import java.util.Timer;
import java.util.TimerTask;

import javax.imageio.ImageIO;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.embed.swing.SwingFXUtils;
import javafx.geometry.Side;
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.control.Label;
import javafx.scene.layout.AnchorPane;
import javafx.scene.layout.StackPane;
import javafx.stage.Stage;
import maps.Meshed;
import maps.ArbitraryPseudocylindrical;
import maps.Azimuthal;
import maps.Cylindrical;
import maps.EqualEarth;
import maps.Lenticular;
import maps.Misc;
import maps.Octohedral;
import maps.Polyhedral;
import maps.Projection;
import maps.Pseudocylindrical;
import maps.WinkelTripel;

/**
 * A simple script that creates an annotated ScatterPlot of map projections
 * 
 * @author jkunimune
 */
public class MapPlotter extends Application {

	private static final double DECIBEL = Math.log(10)/10;
	
	private static final double GLOBE_RES = .005;
	
	private static final Projection[] CYLINDRICAL = { Cylindrical.MERCATOR,
			Cylindrical.PLATE_CARREE, Cylindrical.GALL_ORTHOGRAPHIC,
			Cylindrical.GALL_STEREOGRAPHIC };
	private static final Projection[] AZIMUTHAL = { Azimuthal.POLAR };
	private static final Projection[] PSEUDOCYL = { Pseudocylindrical.MOLLWEIDE,
			ArbitraryPseudocylindrical.ROBINSON, ArbitraryPseudocylindrical.NATURAL_EARTH,
			Pseudocylindrical.KAVRAYSKIY_VII, EqualEarth.EQUAL_EARTH };
	private static final Projection[] LENTICULAR = { Lenticular.AITOFF, Lenticular.VAN_DER_GRINTEN,
			WinkelTripel.WINKEL_TRIPEL, Meshed.DANSEIJI_N, Meshed.DANSEIJI_I,
			Meshed.DANSEIJI_II };
	private static final Projection[] TETRAHEDRAL = { Polyhedral.LEE_TETRAHEDRAL_RECTANGULAR,
			Polyhedral.AUTHAGRAPH, Polyhedral.VAN_LEEUWEN };
	private static final Projection[] CHEATY = { Pseudocylindrical.LEMONS, Octohedral.CAHILL_KEYES,
			Polyhedral.DYMAXION, Octohedral.CAHILL_CONCIALDI, Meshed.DANSEIJI_V };
	private static final Projection[] OTHER = { Misc.PEIRCE_QUINCUNCIAL };
	
	
	private StackPane stack;
	
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage root) {
		final ScatterChart<Number, Number> plot =
				new ScatterChart<Number, Number>(
				new NumberAxis("RMS scale distortion", 0, 4, .5),
				new NumberAxis("RMS shape distortion", 0, 4, .5));
		plot.setLegendSide(Side.RIGHT);
		final AnchorPane overlay = new AnchorPane();
		stack = new StackPane(plot, overlay);
		
		final List<Label> labels = new LinkedList<Label>();
		final List<Data<Number,Number>> data = new LinkedList<Data<Number,Number>>();
		final double[][][] points = Projection.globe(GLOBE_RES);
		
		plotProjections(plot, overlay, labels, data, AZIMUTHAL, "Azimuthal ", points);
		plotProjections(plot, overlay, labels, data, CYLINDRICAL, "Cylindrical ", points);
		plotProjections(plot, overlay, labels, data, PSEUDOCYL, "Pseudocylindrical ", points);
		plotProjections(plot, overlay, labels, data, LENTICULAR, "Lenticular ", points);
		plotProjections(plot, overlay, labels, data, TETRAHEDRAL, "Tetrahedral ", points);
		plotProjections(plot, overlay, labels, data, CHEATY, "Interrupted ", points);
		plotProjections(plot, overlay, labels, data, OTHER, "Other ", points);
		
		final ChangeListener<Number> listener = new ChangeListener<Number>() {
			final Timer timer = new Timer();
			TimerTask task = null;
			@Override
			public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
				if (task != null) 	task.cancel();
				task = new TimerTask() {
					public void run() {
						drawLabelsAndSave(plot, labels, data);
					}
				};
				timer.schedule(task, 1000);
			}
		};
		plot.widthProperty().addListener(listener);
		plot.heightProperty().addListener(listener);
		
		root.setTitle("Map Projections");
		root.setScene(new Scene(stack));
		root.show();
	}
	
	
	private static void plotProjections(ScatterChart<Number, Number> chart,
			AnchorPane overlay, List<Label> labels, List<Data<Number,Number>> data,
			Projection[] projections, String name, double[][][] points) {
		final Series<Number, Number> series = new Series<Number, Number>();
		series.setName(name);
		
		for (Projection projection: projections) {
			System.out.print(projection+": ");
			final double[] params = projection.getDefaultParameters();
			final double distortion[] = projection.avgDistortion(points, params);
			final Data<Number, Number> datum = new Data<Number, Number>(
					distortion[0]/DECIBEL, distortion[1]/DECIBEL);
			series.getData().add(datum);
			final Label lbl = new Label(projection.getName());
			overlay.getChildren().add(lbl);
			labels.add(lbl);
			data.add(datum);
			System.out.println(distortion[0]+", "+distortion[1]);
		}
		
		chart.getData().add(series);
	}
	
	
	private void drawLabelsAndSave(ScatterChart<Number,Number> chart,
			List<Label> labels, List<Data<Number,Number>> data) {
		for (int i = 0; i < labels.size(); i ++) {
			AnchorPane.setLeftAnchor(labels.get(i),
					chart.getXAxis().localToParent(chart.getXAxis().getDisplayPosition(data.get(i).getXValue()), 0).getX() + chart.getPadding().getLeft() + 3);
			AnchorPane.setTopAnchor(labels.get(i),
					chart.getYAxis().localToParent(0, chart.getYAxis().getDisplayPosition(data.get(i).getYValue())).getY() + chart.getPadding().getTop() - labels.get(i).getHeight()
					);
		}
		Platform.runLater(() -> {
			try {
				ImageIO.write(
						SwingFXUtils.fromFXImage(stack.snapshot(null, null), null),
						"png", new File("output/graph - plotter.png"));
			} catch (IOException e) {
				e.printStackTrace();
			}
		});
	}

}
