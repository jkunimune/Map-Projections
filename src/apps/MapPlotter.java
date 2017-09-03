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
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.control.Label;
import javafx.scene.layout.AnchorPane;
import javafx.scene.layout.StackPane;
import javafx.stage.Stage;
import maps.Azimuthal;
import maps.Cylindrical;
import maps.Misc;
import maps.MyProjections;
import maps.Projection;
import maps.Robinson;
import maps.Tetrahedral;
import maps.Tobler;
import maps.WinkelTripel;

/**
 * A simple script that creates an annotated ScatterPlot of map projections
 * 
 * @author jkunimune
 */
public class MapPlotter extends Application {

	private static final double GLOBE_RES = .002;
	
	private static final Projection[] COMMON = { Azimuthal.STEREOGRAPHIC, Cylindrical.MERCATOR,
			Azimuthal.POLAR, Cylindrical.PLATE_CARREE, Azimuthal.EQUAL_AREA,
			Cylindrical.GALL_PETERS, Cylindrical.GALL, Misc.VAN_DER_GRINTEN, Robinson.ROBINSON,
			WinkelTripel.WINKEL_TRIPEL, Azimuthal.ORTHOGRAPHIC };
	private static final Projection[] UNCOMMON = { Misc.PEIRCE_QUINCUNCIAL, Tetrahedral.LEE,
			Tobler.TOBLER, Cylindrical.BEHRMANN, Misc.AITOFF };
	private static final Projection[] INVENTED = { MyProjections.PSEUDOSTEREOGRAPHIC,
			MyProjections.HYPERELLIPOWER, Tetrahedral.TETRAPOWER, Tetrahedral.TETRAFILLET,
			MyProjections.TWO_POINT_EQUALIZED };
	
	
	private StackPane stack;
	
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage root) {
		final ScatterChart<Number, Number> plot =
				new ScatterChart<Number, Number>(
				new NumberAxis("Size distortion", 0, 2.4, 0.2),
				new NumberAxis("Shape distortion", 0, 1.2, 0.2));
		final AnchorPane overlay = new AnchorPane();
		stack = new StackPane(plot, overlay);
		
		final List<Label> labels = new LinkedList<Label>();
		final List<Data<Number,Number>> data = new LinkedList<Data<Number,Number>>();
		final double[][][] points = Projection.globe(GLOBE_RES);
		
		plotProjections(plot, overlay, labels, data, COMMON, "Common", points);
		plotProjections(plot, overlay, labels, data, UNCOMMON, "Hipster", points);
		plotProjections(plot, overlay, labels, data, INVENTED, "Mine", points);
		
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
			final Data<Number, Number> datum = new Data<Number, Number>(distortion[0], distortion[1]);
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
					chart.getXAxis().localToParent(chart.getXAxis().getDisplayPosition(data.get(i).getXValue()), 0).getX() + chart.getPadding().getLeft()
					+ 3);
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
