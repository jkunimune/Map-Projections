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
import java.util.ArrayList;
import java.util.List;
import java.util.function.UnaryOperator;

import javax.imageio.ImageIO;

import dialogs.ProgressBarDialog;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.concurrent.Task;
import javafx.embed.swing.SwingFXUtils;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.image.PixelWriter;
import javafx.scene.image.WritableImage;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
import javafx.scene.input.KeyEvent;
import javafx.scene.layout.HBox;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import maps.Projection;
import util.Math2;

/**
 * An application to analyze the characteristics of map projections
 * 
 * @author Justin Kunimune
 */
public class MapAnalyzer extends Application {

	private static final int CONT_WIDTH = 300;
	private static final int IMG_WIDTH = 500;
	private static final int CHART_WIDTH = 400;
	
	
	private static final KeyCombination ctrlS = new KeyCodeCombination(KeyCode.S, KeyCodeCombination.CONTROL_DOWN);
	
	
	private static final Projection[] PROJ_ARR = { Projection.MERCATOR, Projection.EQUIRECTANGULAR,
			Projection.GALL_PETERS, Projection.HOBODYER, Projection.LAMBERT_CYLIND, Projection.GALL,
			Projection.STEREOGRAPHIC, Projection.POLAR, Projection.E_A_AZIMUTH, Projection.ORTHOGRAPHIC,
			Projection.GNOMONIC, Projection.LAMBERT_CONIC, Projection.E_D_CONIC, Projection.ALBERS, Projection.LEE,
			Projection.TETRAGRAPH, Projection.SINUSOIDAL, Projection.MOLLWEIDE, Projection.HAMMER, Projection.TOBLER,
			Projection.VAN_DER_GRINTEN, Projection.ROBINSON, Projection.WINKEL_TRIPEL, Projection.PEIRCE_QUINCUNCIAL,
			Projection.GUYOU, Projection.MAGNIFIER, Projection.EXPERIMENT };
	
	
	private Stage stage;
	private FileChooser saver;
	private ComboBox<Projection> projectionChooser;
	private Text projectionDesc;
	private Button calculate, saveMap, saveCharts;
	private Label avgSizeDistort, avgShapeDistort;
	private ImageView output;
	private VBox charts;
	private BarChart<String, Number> sizeChart, shapeChart;
	
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage primaryStage) {
		stage = primaryStage;
		stage.setTitle("Map Analyzer");
		
		final VBox layout = new VBox();
		layout.setSpacing(5);
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(CONT_WIDTH);
		
		Label lbl = new Label("Projection:");
		ObservableList<Projection> items = FXCollections.observableArrayList(PROJ_ARR);
		projectionChooser = new ComboBox<Projection>(items);
		projectionChooser.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				projectionDesc.setText(projectionChooser.getValue().getDescription());
			}
		});
		projectionChooser.setValue(Projection.MERCATOR);
		layout.getChildren().add(new HBox(3, lbl, projectionChooser));
		
		projectionDesc = new Text(projectionChooser.getValue().getDescription());
		projectionDesc.setWrappingWidth(CONT_WIDTH);
		layout.getChildren().add(projectionDesc);
		
		calculate = new Button("Calculate");
		calculate.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				calculateMap();
			}
		});
		calculate.setTooltip(new Tooltip(
				"Calculate the distortion for this map."));
		calculate.setDefaultButton(true);
		
		saver = new FileChooser();
		saver.setInitialDirectory(new File("output"));
		saver.setInitialFileName("myMap.jpg");
		saver.setTitle("Save Image");
		saver.getExtensionFilters().addAll(
				new FileChooser.ExtensionFilter("JPG", "*.jpg"),
				new FileChooser.ExtensionFilter("PNG", "*.png"));
		
		saveMap = new Button("Save Image...");
		saveMap.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				startFinalizingMap();
			}
		});
		saveMap.setTooltip(new Tooltip("Save the distortion graphic."));
		stage.addEventHandler(KeyEvent.KEY_RELEASED, new EventHandler<KeyEvent>() {	//ctrl-S saves
			public void handle(KeyEvent event) {
				if (ctrlS.match(event))	saveMap.fire();
			}
		});
		
		saveCharts = new Button("Save Chart...");
		saveCharts.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				saveImage(charts.snapshot(null, null), null);
			}
		});
		
		HBox box = new HBox(5, calculate, saveMap, saveCharts);
		box.setAlignment(Pos.CENTER);
		layout.getChildren().add(box);
		
		layout.getChildren().add(new Separator());
		
		avgSizeDistort = new Label("...");
		avgShapeDistort = new Label("..."); //TODO: reorder?
		lbl = new Label("Blue areas are dilated, red areas are compressed, and black areas are stretched.");
		lbl.setWrapText(true);
		
		VBox bxo = new VBox(3,
				new HBox(new Label("Average size distortion: "),avgSizeDistort),
				new HBox(new Label("Average shape distortion: "),avgShapeDistort),
				lbl);
		bxo.setAlignment(Pos.CENTER_LEFT);
		layout.getChildren().add(bxo);
		
		output = new ImageView();
		output.setFitWidth(IMG_WIDTH);
		output.setFitHeight(IMG_WIDTH);
		output.setPreserveRatio(true);
		
		sizeChart = new BarChart<String, Number>(new CategoryAxis(), new NumberAxis());
		sizeChart.setPrefWidth(CHART_WIDTH);
		sizeChart.setPrefHeight(IMG_WIDTH/2);
		sizeChart.getXAxis().setLabel("Scale factor");
		sizeChart.setBarGap(0);
		sizeChart.setCategoryGap(0);
		sizeChart.setAnimated(false);
		sizeChart.setLegendVisible(false);
		
		shapeChart = new BarChart<String, Number>(new CategoryAxis(), new NumberAxis());
		shapeChart.setPrefWidth(CHART_WIDTH);
		shapeChart.setPrefHeight(IMG_WIDTH/2);
		shapeChart.getXAxis().setLabel("Stretch factor");
		shapeChart.setBarGap(0);
		shapeChart.setCategoryGap(0);
		shapeChart.setAnimated(false);
		shapeChart.setLegendVisible(false);
		
		charts = new VBox(sizeChart, shapeChart);
		
		final HBox gui = new HBox(layout, output, charts);
		
		new Thread(() -> {
			calculate.fire();
		}).start();
		
		gui.setAlignment(Pos.CENTER);
		gui.setSpacing(10);
		StackPane.setMargin(gui, new Insets(10));
		stage.setScene(new Scene(new StackPane(gui)));
		stage.show();
	}
	
	
	private void calculateMap() {
		calculate.setDisable(true);
		new Thread(new Task<Void>() {
			protected Void call() {
				try {
					Platform.runLater(() -> {
						sizeChart.getData().clear();
						shapeChart.getData().clear();
						
						avgSizeDistort.setText("...");
						avgShapeDistort.setText("...");
					});
					
					final Projection p = projectionChooser.getValue();
					final double[][][] distortionM =
							calculateDistortion(map(250, p), p::project);
					
					output.setImage(makeGraphic(distortionM));
					
					final double[][][] distortionG =
							calculateDistortion(globe(0.02), p::project);
					
					Platform.runLater(() -> {
						sizeChart.getData().add(histogram(distortionG[0],
								-2,2,14, true));
						shapeChart.getData().add(histogram(distortionG[1],
								0,1.6,14, false));
						
						avgSizeDistort.setText(format(Math2.stdDev(distortionG[0])));
						avgShapeDistort.setText(format(Math2.mean(distortionG[1])));
					});
					calculate.setDisable(false);
					return null;
				} catch (Exception e) {
					e.printStackTrace();
					return null;
				}
			}
		}).start();
	}
	
	
	private void startFinalizingMap() {
		final Projection p = projectionChooser.getValue();
		
		ProgressBarDialog pBar = new ProgressBarDialog();
		pBar.show();
		new Thread(() -> {
			final double[][][] distortion = calculateDistortion(
					map(1000, p), p::project, pBar);
			Image graphic = makeGraphic(distortion);
			Platform.runLater(() -> saveImage(graphic, pBar));
		}).start();
	}
	
	
	public static double[][][] calculateDistortion(double[][][] points, UnaryOperator<double[]> p) {
		return calculateDistortion(points, p, null);
	}
	
	
	public static double[][][] calculateDistortion(double[][][] points, UnaryOperator<double[]> p,
			ProgressBarDialog pBar) { //calculate both kinds of distortion over the given region
		double[][][] output = new double[2][points.length][points[0].length]; //the distortion matrix
		
		for (int y = 0; y < points.length; y ++) {
			for (int x = 0; x < points[y].length; x ++) {
				if (points[y][x] != null) {
					final double[] distortions = getDistortionAt(points[y][x], p);
					output[0][y][x] = distortions[0]; //the output matrix has two layers:
					output[1][y][x] = distortions[1]; //area and angular distortion
				}
				else {
					output[0][y][x] = Double.NaN;
					output[1][y][x] = Double.NaN; //NaN means no map here
				}
			}
			if (pBar != null)
				pBar.setProgress((double)(y+1)/points.length);
		}
		
		final double avgArea = Math2.mean(output[0]); //don't forget to normalize output[0] so the average is zero
		for (int y = 0; y < output[0].length; y ++)
			for (int x = 0; x < output[0][y].length; x ++)
				output[0][y][x] -= avgArea;
		
		return output;
	}
	
	
	private static double[] getDistortionAt(double[] s0, UnaryOperator<double[]> p) { //calculate both kinds of distortion at the given point
		final double[] output = new double[2];
		final double dx = 1e-6;
		
		final double[] s1 = { s0[0], s0[1]+dx/Math.cos(s0[0]) }; //consider a point slightly to the east
		final double[] s2 = { s0[0]+dx, s0[1] }; //and slightly to the north
		final double[] p0 = p.apply(s0);
		final double[] p1 = p.apply(s1);
		final double[] p2 = p.apply(s2);
		
		final double dA = 
				(p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]);
		output[0] = Math.log(Math.abs(dA/(dx*dx))); //the zeroth output is the size (area) distortion
		if (Math.abs(output[0]) > 25)
			output[0] = Double.NaN; //discard outliers
		
		final double s1ps2 = Math.hypot((p1[0]-p0[0])+(p2[1]-p0[1]), (p1[1]-p0[1])-(p2[0]-p0[0]));
		final double s1ms2 = Math.hypot((p1[0]-p0[0])-(p2[1]-p0[1]), (p1[1]-p0[1])+(p2[0]-p0[0]));
		final double factor = Math.abs((s1ps2-s1ms2)/(s1ps2+s1ms2)); //there's some linear algebra behind this formula. Don't worry about it.
		
		//output[1] = Math.asin(1-factor);
		output[1] = (1-factor)/Math.sqrt(factor);
		//if (output[1] >= Math.PI/2)
		//	output[1] = Double.NaN;
		
		return output;
	}
	
	
	private static Image makeGraphic(double[][][] distortion) {
		WritableImage output = new WritableImage(distortion[0][0].length, distortion[0].length);
		PixelWriter writer = output.getPixelWriter();
		for (int y = 0; y < distortion[0].length; y ++) {
			for (int x = 0; x < distortion[0][y].length; x ++) {
				final double sizeDistort = distortion[0][y][x];
				final double shapeDistort = distortion[1][y][x];
				if (Double.isNaN(sizeDistort) || Double.isNaN(shapeDistort)) {
					writer.setArgb(x, y, 0);
					continue;
				}
				
				final int r, g, b;
				if (sizeDistort < 0) { //if compressing
					r = (int)(255.9*Math.exp(-shapeDistort/1.5));
					g = (int)(255.9*Math.exp(-shapeDistort/1.5)*Math.exp(sizeDistort/1.5));
					b = g;
				}
				else { //if dilating
					r = (int)(255.9*Math.exp(-shapeDistort/1.5)*Math.exp(-sizeDistort/1.5));
					g = r;
					b = (int)(255.9*Math.exp(-shapeDistort/1.5));
				}
				
				final int argb = ((((((0xFF)<<8)+r)<<8)+g)<<8)+b;
				writer.setArgb(x, y, argb);
			}
		}
		return output;
	}
	
	
	private void saveImage(Image img, ProgressBarDialog pBar) { // call from the main thread!
		if (pBar != null)
			pBar.close();
		
		final File f = saver.showSaveDialog(stage);
		if (f != null) {
			new Thread(() -> {
				try {
					saveMap.setDisable(true);
					ImageIO.write(SwingFXUtils.fromFXImage(img,null), "png", f);
					saveMap.setDisable(false);
				} catch (IOException e) {}
			}).start();
		}
	}


	public static double[][][] map(int size, Projection proj) { //generate a matrix of coordinates based on a map projection
		final int w = size, h = (int)(size/proj.getAspectRatio());
		double[][][] output = new double[h][w][2]; //the coordinate matrix
		
		for (int y = 0; y < output.length; y ++)
			for (int x = 0; x < output[y].length; x ++)
				output[y][x] = proj.inverse(2.*(x+.5)/w - 1, 1 - 2.*(y+.5)/h); //s0 is this point on the sphere
		
		return output;
	}
	
	
	public static double[][][] globe(double dt) { //generate a matrix of coordinates based on the sphere
		List<double[]> points = new ArrayList<double[]>();
		for (double phi = -Math.PI/2+dt/2; phi < Math.PI/2; phi += dt) { // make sure phi is never exactly +-tau/4
			for (double lam = -Math.PI; lam < Math.PI; lam += dt/Math.cos(phi)) {
				points.add(new double[] {phi, lam});
			}
		}
		return new double[][][] {points.toArray(new double[0][])};
	}
	
	
	private static final Series<String, Number> histogram(double[][] values,
			double min, double max, int num, boolean logarithmic) {
		int[] hist = new int[num+1]; //this array is the histogram values for min, min+dx, ..., max-dx, max
		int tot = 0;
		for (double[] row: values) {
			for (double x: row) {
				if (Double.isFinite(x)) {
					final int i = (int)Math.round((x-min)/(max-min)*num);
					if (i >= 0 && i <= num)
						hist[i] ++;
					tot ++;
				}
			}
		}
		Series<String, Number> output = new Series<String, Number>();
		for (int i = 0; i <= num; i ++) {
			double x = i*(max-min)/num+min;
			if (logarithmic)	x = Math.exp(x);
			else				x = 1/(1-Math.sin(x)); //this is a bit nonsensical and sketch. Don't worry about it.
			output.getData().add(new Data<String, Number>(
					Double.toString(Math.round(100*x)/100.),
					(double)hist[i]/tot*100));
		}
		return output;
	}
	
	
	
	private static final String format(double d) {
		if (d < 1000)
			return Double.toString(Math.round(d*1000.)/1000.);
		else
			return "1000+";
	}

}