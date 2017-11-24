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
import java.util.function.DoubleUnaryOperator;

import javax.imageio.ImageIO;

import dialogs.ProgressBarDialog;
import javafx.application.Platform;
import javafx.embed.swing.SwingFXUtils;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.image.Image;
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
import maps.Azimuthal;
import maps.Conic;
import maps.Cylindrical;
import maps.Misc;
import maps.MyProjections;
import maps.Projection;
import maps.Pseudocylindrical;
import maps.Robinson;
import maps.Tetrahedral;
import maps.Tobler;
import maps.WinkelTripel;
import utils.Math2;
import utils.Procedure;

/**
 * An application to analyze the characteristics of map projections
 * 
 * @author Justin Kunimune
 */
public class MapAnalyzer extends MapApplication {

	public static final void main(String[] args) {
		launch(args);
	}
	
	
	
	private static final double LN_10 = Math.log(10);
	
	private static final int CHART_WIDTH = 400;
	private static final int ROUGH_SAMP_NUM = 500;
	private static final int FINE_SAMP_NUM = 1000;
	private static final double GLOBE_RES = .01;
	
	private static final FileChooser.ExtensionFilter[] RASTER_TYPES = {
			new FileChooser.ExtensionFilter("PNG", "*.png"),
			new FileChooser.ExtensionFilter("JPG", "*.jpg","*.jpeg","*.jpe","*.jfif"),
			new FileChooser.ExtensionFilter("GIF", "*.gif") };
	
	private static final Projection[] PROJ_ARR = { Cylindrical.MERCATOR, Cylindrical.PLATE_CARREE,
			Cylindrical.EQUIRECTANGULAR, Cylindrical.GALL_PETERS, Cylindrical.HOBO_DYER,
			Cylindrical.BEHRMANN, Cylindrical.LAMBERT, Cylindrical.EQUAL_AREA, Cylindrical.GALL,
			Azimuthal.STEREOGRAPHIC, Azimuthal.POLAR, Azimuthal.EQUAL_AREA, Azimuthal.GNOMONIC,
			Azimuthal.ORTHOGRAPHIC, Azimuthal.PERSPECTIVE, Conic.LAMBERT, Conic.EQUIDISTANT,
			Conic.ALBERS, Tetrahedral.LEE, Tetrahedral.TETRAGRAPH, Tetrahedral.ACTUAUTHAGRAPH,
			Tetrahedral.AUTHAGRAPH, Pseudocylindrical.SINUSOIDAL, Pseudocylindrical.MOLLWEIDE,
			Tobler.TOBLER, Misc.HAMMER, Misc.AITOFF, Misc.VAN_DER_GRINTEN, Robinson.ROBINSON,
			WinkelTripel.WINKEL_TRIPEL, Misc.PEIRCE_QUINCUNCIAL, Misc.TWO_POINT_EQUIDISTANT,
			Misc.HAMMER_RETROAZIMUTHAL, Pseudocylindrical.LEMONS, MyProjections.EXPERIMENT,
			MyProjections.PSEUDOSTEREOGRAPHIC, MyProjections.HYPERELLIPOWER, Tetrahedral.TETRAPOWER,
			Tetrahedral.TETRAFILLET, MyProjections.TWO_POINT_EQUALIZED };
	
	
	private Button updateBtn;
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
		new Thread(() -> {
			calculateAndUpdate();
		}).start();
	}
	
	
	@Override
	protected Node makeWidgets() {
		final Node projectionSelector = buildProjectionSelector(PROJ_ARR, Procedure.NONE);
		final Node parameterSelector = buildParameterSelector(Procedure.NONE);
		final Node textDisplay = buildTextDisplay();
		this.updateBtn = buildUpdateButton(this::calculateAndUpdate);
		this.updateBtn.setText("Calculate"); //I don't need to follow your darn conventions!
		final Button saveMapBtn = buildSaveButton(true, "map", RASTER_TYPES,
				RASTER_TYPES[0], ()->true, this::calculateAndSaveMap);
		final Button savePltBtn = buildSaveButton(true, "plots", RASTER_TYPES,
				RASTER_TYPES[0], ()->true, this::calculateAndSavePlot);
		final HBox buttons = new HBox(5, updateBtn, saveMapBtn, savePltBtn);
		buttons.setAlignment(Pos.CENTER);
		
		final VBox layout = new VBox(5,
				projectionSelector, parameterSelector, new Separator(),
				buttons, new Separator(), textDisplay);
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(GUI_WIDTH);
		
		this.mapDisplay = new ImageView();
		this.mapDisplay.setFitWidth(IMG_WIDTH);
		this.mapDisplay.setFitHeight(IMG_WIDTH);
		this.mapDisplay.setPreserveRatio(true);
		final StackPane pane = new StackPane(mapDisplay);
		pane.setMinWidth(IMG_WIDTH);
		
		this.charts = buildDistortionHistograms();
		
		final HBox gui = new HBox(10, layout, pane, charts);
		gui.setAlignment(Pos.CENTER);
		StackPane.setMargin(gui, new Insets(10));
		
		return gui;
	}
	
	
	private Node buildTextDisplay() {
		this.avgSizeDistort = new Text("...");
		this.avgShapeDistort = new Text("...");
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
		this.sizeChart.setPrefHeight(IMG_WIDTH/2);
		this.sizeChart.getXAxis().setLabel("Scale factor");
		this.sizeChart.setBarGap(0);
		this.sizeChart.setCategoryGap(0);
		this.sizeChart.setAnimated(false);
		this.sizeChart.setLegendVisible(false);
		
		this.shapeChart = new BarChart<String, Number>(new CategoryAxis(), new NumberAxis());
		this.shapeChart.setPrefWidth(CHART_WIDTH);
		this.shapeChart.setPrefHeight(IMG_WIDTH/2);
		this.shapeChart.getXAxis().setLabel("Stretch factor");
		this.shapeChart.setBarGap(0);
		this.shapeChart.setCategoryGap(0);
		this.shapeChart.setAnimated(false);
		this.shapeChart.setLegendVisible(false);
		
		return new VBox(5, sizeChart, shapeChart);
	}
	
	
	private void calculateAndUpdate() {
		Platform.runLater(() -> {
			sizeChart.getData().clear();
			shapeChart.getData().clear();
			
			avgSizeDistort.setText("...");
			avgShapeDistort.setText("...");
		});
		
		loadParameters();
		final Projection proj = this.getProjection();
		final double[][][] distortionM = proj.calculateDistortion(proj.map(ROUGH_SAMP_NUM));
		
		mapDisplay.setImage(makeGraphic(distortionM));
		
		final double[][][] distortionG = proj.calculateDistortion(Projection.globe(GLOBE_RES));
		
		Platform.runLater(() -> {
				sizeChart.getData().add(histogram(distortionG[0],
						-LN_10, LN_10, 14, Math::exp));
				shapeChart.getData().add(histogram(distortionG[1],
						   0.0, LN_10, 14, Math::exp));
				
				avgSizeDistort.setText(format(Math2.stdDev(distortionG[0])));
				avgShapeDistort.setText(format(Math2.mean(distortionG[1])));
			});
	}
	
	
	private void calculateAndSavePlot(File file, ProgressBarDialog pBar) {
		pBar.setProgress(-1);
		final String filename = file.getName();
		final String extension = filename.substring(filename.lastIndexOf('.')+1);
		try {
			final WritableImage out = new WritableImage(
					(int) charts.getWidth(), (int) charts.getHeight());
			Platform.runLater(() -> charts.snapshot(null,out));
			while (out.getProgress() < 1) {}
			ImageIO.write(SwingFXUtils.fromFXImage(out, null), extension, file); //save
		} catch (IOException e) {
			showError("Failure!",
					"Could not access "+file.getAbsolutePath()+". It's possible that another program has it open.");
		}
	}
	
	
	private void calculateAndSaveMap(File file, ProgressBarDialog pBar) {
		loadParameters();
		final Projection proj = this.getProjection();
		final double[][][] distortion = proj.calculateDistortion(proj.map(FINE_SAMP_NUM), pBar); //calculate
		Image graphic = makeGraphic(distortion);
		
		pBar.setProgress(-1);
		
		final String filename = file.getName();
		final String extension = filename.substring(filename.lastIndexOf('.')+1);
		try {
			ImageIO.write(SwingFXUtils.fromFXImage(graphic,null), extension, file); //save
		} catch (IOException e) {
			showError("Failure!",
					"Could not access "+file.getAbsolutePath()+". It's possible that another program has it open.");
		}
	}
	
	
	private static Image makeGraphic(double[][][] distortion) {
		WritableImage output = new WritableImage(distortion[0][0].length, distortion[0].length);
		PixelWriter writer = output.getPixelWriter();
		for (int y = 0; y < distortion[0].length; y ++) {
			for (int x = 0; x < distortion[0][y].length; x ++) {
				final double sizeDistort = distortion[0][y][x], shapeDistort = distortion[1][y][x];
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
		return output;
	}
	
	
	private static final Series<String, Number> histogram(double[][] values,
			double min, double max, int num, DoubleUnaryOperator labeler) {
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
			double x = labeler.applyAsDouble(i*(max-min)/num+min);
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