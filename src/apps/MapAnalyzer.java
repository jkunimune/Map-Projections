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
			Projection.GALL_PETERS, Projection.HOBODYER, Projection.BEHRMANN, Projection.LAMBERT_CYLIND,
			Projection.GALL, Projection.STEREOGRAPHIC, Projection.POLAR, Projection.E_A_AZIMUTH,
			Projection.ORTHOGRAPHIC, Projection.GNOMONIC, Projection.LAMBERT_CONIC, Projection.E_D_CONIC,
			Projection.ALBERS, Projection.LEE, Projection.TETRAGRAPH, Projection.SINUSOIDAL, Projection.MOLLWEIDE,
			Projection.HAMMER, Projection.TOBLER, Projection.AITOFF, Projection.VAN_DER_GRINTEN, Projection.ROBINSON,
			Projection.WINKEL_TRIPEL, Projection.PEIRCE_QUINCUNCIAL, Projection.GUYOU, Projection.MAGNIFIER,
			Projection.EXPERIMENT };
	
	
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
					
					final Projection proj = projectionChooser.getValue();
					final double[][][] distortionM =
							proj.calculateDistortion(proj.map(250));
					
					output.setImage(makeGraphic(distortionM));
					
					final double[][][] distortionG =
							proj.calculateDistortion(Projection.globe(0.02));
					
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
			final double[][][] distortion = p.calculateDistortion(
					p.map(1000), pBar);
			Image graphic = makeGraphic(distortion);
			Platform.runLater(() -> saveImage(graphic, pBar));
		}).start();
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