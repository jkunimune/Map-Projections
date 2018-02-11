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
import java.util.function.DoubleUnaryOperator;

import javafx.concurrent.Task;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.SnapshotResult;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.image.ImageView;
import javafx.scene.image.PixelWriter;
import javafx.scene.image.WritableImage;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import maps.Projection;
import utils.Flag;
import utils.Math2;
import utils.MutableDouble;
import utils.Procedure;
import utils.SavableImage;

/**
 * An application to analyse the characteristics of map projections
 * 
 * @author Justin Kunimune
 */
public class MapAnalyzer extends MapApplication {

	public static final void main(String[] args) {
		launch(args);
	}
	
	
	
	private static final double LN_10 = Math.log(10);
	
	private static final int CHART_WIDTH = 420;
	private static final int FINE_SAMP_NUM = 2048;
	private static final double GLOBE_RES = .01;
	
	private static final FileChooser.ExtensionFilter[] RASTER_TYPES = {
			new FileChooser.ExtensionFilter("PNG", "*.png"),
			new FileChooser.ExtensionFilter("JPG", "*.jpg","*.jpeg","*.jpe","*.jfif"),
			new FileChooser.ExtensionFilter("GIF", "*.gif") };
	
	private Flag cropAtIDL;
	private MutableDouble graticuleSpacing;
	private Text avgSizeDistort, avgShapeDistort;
	private ImageView mapDisplay;
	private Region charts;
	private BarChart<String, Number> sizeChart;
	private BarChart<String, Number> shapeChart;
	
	
	
	public MapAnalyzer() {
		super("Map Analyzer");
	}
	
	
	@Override
	public void start(Stage root) {
		super.start(root);
		new Thread(calculateGraphicForUpdate()).start();
	}
	
	
	@Override
	protected Region makeWidgets() {
		this.cropAtIDL = new Flag(true);
		this.graticuleSpacing = new MutableDouble();
		final Region projectionSelector = buildProjectionSelector(Procedure.NONE);
		final Region parameterSelector = buildParameterSelector(Procedure.NONE);
		final Region optionPane = buildOptionPane(cropAtIDL, graticuleSpacing);
		final Region textDisplay = buildTextDisplay();
		final Region updateBtn = buildUpdateButton("Calculate", this::calculateGraphicForUpdate);
		final Region saveMapBtn = buildSaveButton(true, "map", RASTER_TYPES,
				RASTER_TYPES[0], ()->true, this::calculateGraphicForSaving);
		final Region savePltBtn = buildSaveButton(true, "plots", RASTER_TYPES,
				RASTER_TYPES[0], ()->true, this::calculatePlot);
		final HBox buttons = new HBox(H_SPACE, updateBtn, saveMapBtn, savePltBtn);
		buttons.setAlignment(Pos.CENTER);
		
		final VBox layout = new VBox(V_SPACE,
				projectionSelector, parameterSelector, new Separator(),
				optionPane, new Separator(), buttons, new Separator(), textDisplay);
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(GUI_WIDTH);
		
		this.mapDisplay = new ImageView();
		this.mapDisplay.setFitWidth(IMG_SIZE);
		this.mapDisplay.setFitHeight(IMG_SIZE);
		this.mapDisplay.setPreserveRatio(true);
		final StackPane pane = new StackPane(mapDisplay);
		pane.setMinWidth(IMG_SIZE);
		
		this.charts = buildDistortionHistograms();
		
		final HBox gui = new HBox(MARGIN, layout, pane, charts);
		gui.setAlignment(Pos.CENTER);
		StackPane.setMargin(gui, new Insets(MARGIN));
		
		return gui;
	}
	
	
	private Region buildTextDisplay() {
		this.avgSizeDistort = new Text("\u2026");
		this.avgShapeDistort = new Text("\u2026");
		final Text txt = new Text("Blue areas are dilated, red areas are compressed, and black areas are stretched.");
		txt.setWrappingWidth(GUI_WIDTH);
		
		VBox box = new VBox(3,
				new HBox(new Label("Average size distortion: "),avgSizeDistort),
				new HBox(new Label("Average shape distortion: "),avgShapeDistort),
				txt);
		box.setAlignment(Pos.CENTER_LEFT);
		return box;
	}
	
	
	private Region buildDistortionHistograms() {
		this.sizeChart = new BarChart<String, Number>(new CategoryAxis(), new NumberAxis());
		this.sizeChart.setPrefWidth(CHART_WIDTH);
		this.sizeChart.setPrefHeight(IMG_SIZE/2);
		this.sizeChart.getXAxis().setLabel("Scale factor");
		this.sizeChart.setBarGap(0);
		this.sizeChart.setCategoryGap(0);
		this.sizeChart.setAnimated(false);
		this.sizeChart.setLegendVisible(false);
		
		this.shapeChart = new BarChart<String, Number>(new CategoryAxis(), new NumberAxis());
		this.shapeChart.setPrefWidth(CHART_WIDTH);
		this.shapeChart.setPrefHeight(IMG_SIZE/2);
		this.shapeChart.getXAxis().setLabel("Stretch factor");
		this.shapeChart.setBarGap(0);
		this.shapeChart.setCategoryGap(0);
		this.shapeChart.setAnimated(false);
		this.shapeChart.setLegendVisible(false);
		
		return new VBox(5, sizeChart, shapeChart);
	}
	
	
	private Task<SavableImage> calculatePlot() {
		Task<SavableImage> task = new Task<SavableImage>() {
			private SnapshotResult snapRes = null;
			
			protected void running() {
				charts.snapshot((snapRes) -> {
					this.snapRes = snapRes;
					return null;
				}, null, null);
			}
			
			protected SavableImage call() { //yeah, this is weird and dumb, but only because the
				updateProgress(-1, 1); //snapshot method is dumb. Why does it have to be called on
				updateMessage("Converting graph..."); //the GUI thread? It's clearly doing Thread
				while (!isCancelled() && snapRes == null) {} //things, judging by the fact that it
				return SavableImage.savable(snapRes.getImage()); //_has_ a Callback version
			}
		};
		disableWhile(task.runningProperty(), ButtonType.SAVE_GRAPH);
		return task;
	}
	
	
	private Task<SavableImage> calculateGraphicForUpdate() {
		Task<SavableImage> task = calculateGraphic(IMG_SIZE, true);
		disableWhile(task.runningProperty(), ButtonType.UPDATE_MAP, ButtonType.SAVE_MAP);
		return task;
	}
	
	private Task<SavableImage> calculateGraphicForSaving() {
		Task<SavableImage> task = calculateGraphic(FINE_SAMP_NUM, false);
		disableWhile(task.runningProperty(), ButtonType.SAVE_MAP);
		return task;
	}
	
	private Task<SavableImage> calculateGraphic(int imgSize, boolean detailedAnalysis) {
		loadParameters();
		final Projection proj = getProjection();
		final boolean crop = cropAtIDL.isSet();
		final double gratSpace = graticuleSpacing.get();
		
		return new Task<SavableImage>() {
			double[][][] distortionG; //some variables that might get used later
			double sizeDistort, shapeDistort;
			WritableImage graphic;
			
			protected SavableImage call() {
				updateProgress(-1, 1);
				updateMessage("Calculating distortion...");
				
				double[][][] distortionM = proj.calculateDistortion(proj.map(imgSize, crop),
						this::isCancelled, (p) -> updateProgress(p, 2)); //calculate
				if (detailedAnalysis) {
					distortionG = proj.calculateDistortion(Projection.globe(GLOBE_RES));
					sizeDistort = Math2.stdDev(distortionG[0]);
					shapeDistort = Math2.mean(distortionG[1]);
				}
				
				if (isCancelled()) 	return null;
				updateProgress(1, 2);
				updateMessage("Making graphic...");
				
				graphic = new WritableImage(
						distortionM[0][0].length, distortionM[0].length);
				PixelWriter writer = graphic.getPixelWriter();
				for (int y = 0; y < distortionM[0].length; y ++) {
					if (isCancelled()) 	return null;
					updateProgress(1 + (double)y/distortionM[0].length, 2);
					for (int x = 0; x < distortionM[0][y].length; x ++) {
						final double sizeDistort = distortionM[0][y][x], shapeDistort = distortionM[1][y][x];
						final double sizeContour = Math.round(sizeDistort/(LN_10/10))*LN_10/10; //contour the size by decibels
						final double shapeContour = Math.round(shapeDistort/(LN_10/20))*LN_10/20; //contour the size by semidecibels
						if (Double.isNaN(sizeDistort) || Double.isNaN(shapeDistort)) {
							writer.setArgb(x, y, 0);
							continue;
						}
						
						final int r, g, b;
						if (sizeDistort < 0) { //if compressing
							r = (int)(255.9*Math.exp(-shapeContour*.6));
							g = (int)(255.9*Math.exp(-shapeContour*.6)*Math.exp(sizeContour*.6));
							b = g;
						}
						else { //if dilating
							r = (int)(255.9*Math.exp(-shapeContour*.6)*Math.exp(-sizeContour*.6));
							g = r; //I find .6 to be a rather visually pleasing sensitivity
							b = (int)(255.9*Math.exp(-shapeContour*.6));
						}
						
						final int argb = ((((((0xFF)<<8)+r)<<8)+g)<<8)+b;
						writer.setArgb(x, y, argb);
					}
				}
				
				return SavableImage.savable(graphic);
			}
			
			protected void succeeded() {
				mapDisplay.setImage(graphic);
				
				if (detailedAnalysis) {
					sizeChart.getData().clear();
					sizeChart.getData().add(histogram(distortionG[0],
							-LN_10, LN_10, 20, Math::exp));
					shapeChart.getData().clear();
					shapeChart.getData().add(histogram(distortionG[1],
							   0.0, LN_10, 20, Math::exp));
					
					avgSizeDistort.setText(format(sizeDistort/LN_10*10)+"dB");
					avgShapeDistort.setText(format(shapeDistort/LN_10*10)+"dB");
				}
			}
		};
	}
	
	
	private static final Series<String, Number> histogram(double[][] values,
			double min, double max, int num, DoubleUnaryOperator converter) {
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
			double x = converter.applyAsDouble(i*(max-min)/num+min);
			output.getData().add(new Data<String, Number>(
					Double.toString(Math.round(100*x)/100.),
					(double)hist[i]/tot*100));
		}
		return output;
	}
	
	
	
	private static final String format(double d) {
		if (d < 1000)
			return Double.toString(Math.round(d*100.)/100.);
		else
			return "1000+";
	}

}