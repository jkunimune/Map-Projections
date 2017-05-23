package mapAnalyzer;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

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
import util.ProgressBarDialog;

/**
 * An application to make raster oblique aspects of map projections
 * 
 * @author Justin Kunimune
 */
public class MapAnalyzer extends Application {

	private static final int CONT_WIDTH = 300;
	private static final int IMG_WIDTH = 500;
	private static final int CHART_WIDTH = 400;
	
	
	private static final KeyCombination ctrlS = new KeyCodeCombination(KeyCode.S, KeyCodeCombination.CONTROL_DOWN);
	
	
	private static final String[] PROJ_ARR = { "Equirectangular", "Mercator", "Gall Stereographic",
			"Hobo-Dyer", "Polar", "Stereographic", "Azimuthal Equal-Area", "Orthographic", "Gnomonic",
			"Equidistant Conic", "Conformal Conic", "Albers", "Van der Grinten", "Robinson", "Winkel Tripel", "Mollweide", "Hammer", "Sinusoidal",
			"Pierce Quincuncial", "Guyou", "TetraGraph", "Magnifier", "Experimental" };
	private static final double[] DEFA = { 2, 1, 4/3.0, 1.977, 1, 1, 1, 1, 1,
			2, 2, 2, 1, 1.9716, Math.PI/2, 2, 2, 2, 1, 2, Math.sqrt(3), 1, 1 };
	private static final String[] DESC = { "An equidistant cylindrical map", "A conformal cylindrical map",
			"A compromising cylindrical map", "An equal-area cylindrical map", "An equidistant azimuthal map",
			"A conformal azimuthal map", "An equal-area azimuthal map",
			"Represents earth viewed from an infinite distance",
			"Every straight line on the map is a straight line on the sphere", "An equidistant conic map",
			"A conformal conic map", "An equal-area conic map", "A circular compromise map",
			"A visually pleasing piecewise compromise map", "The compromise map used by National Geographic",
			"An equal-area map shaped like an ellipse", "An equal-area map shaped like an ellipse",
			"An equal-area map shaped like a sinusoid", "A conformal square map that uses complex math",
			"A reorganized version of Pierce Quincuncial and actually the best map ever",
			"A compromising knockoff of the AuthaGraph projection",
			"A novelty map that swells the center to disproportionate scale",
			"What happens when you apply a complex differentiable function to a stereographic projection?" };
	
	
	private Stage stage;
	private FileChooser saver;
	private ComboBox<String> projectionChooser;
	private Text projectionDesc;
	private Button calculate, saveMap;
	private Label avgSizeDistort, avgShapeDistort;
	private ImageView output;
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
		ObservableList<String> items = FXCollections.observableArrayList(PROJ_ARR);
		projectionChooser = new ComboBox<String>(items);
		projectionChooser.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				for (int i = 0; i < PROJ_ARR.length; i ++) {
					if (PROJ_ARR[i].equals(projectionChooser.getValue())) {
						projectionDesc.setText(DESC[i]);
						break;
					}
				}
			}
		});
		projectionChooser.setValue(PROJ_ARR[1]);
		layout.getChildren().add(new HBox(3, lbl, projectionChooser));
		
		projectionDesc = new Text(DESC[1]);
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
		
		HBox box = new HBox(5, calculate, saveMap);
		box.setAlignment(Pos.CENTER);
		layout.getChildren().add(box);
		
		layout.getChildren().add(new Separator());
		
		avgSizeDistort = new Label("...");
		avgShapeDistort = new Label("..."); //TODO: reorder?
		lbl = new Label("Blue areas are dilated, red areas are compressed, and black areas are skewed.");
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
		sizeChart.getXAxis().setLabel("Size Distortion");
		sizeChart.setBarGap(0);
		sizeChart.setCategoryGap(0);
		sizeChart.setAnimated(false);
		
		shapeChart = new BarChart<String, Number>(new CategoryAxis(), new NumberAxis());
		shapeChart.setPrefWidth(CHART_WIDTH);
		shapeChart.setPrefHeight(IMG_WIDTH/2);
		shapeChart.getXAxis().setLabel("Shape Distortion");
		shapeChart.setBarGap(0);
		shapeChart.setCategoryGap(0);
		shapeChart.setAnimated(false);
		
		final HBox gui = new HBox(layout, output, new VBox(sizeChart, shapeChart));
		
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
					
					final String p = projectionChooser.getValue();
					final double[][][] distortionM =
							calculateDistortion(map(250, p), p);
					
					output.setImage(makeGraphic(distortionM));
					
					final double[][][] distortionG =
							calculateDistortion(globe(0.02, p), p);
					
					Platform.runLater(() -> {
						sizeChart.getData().add(histogram(distortionG[0], -2,2,16));
						shapeChart.getData().add(histogram(distortionG[1], 0,1.6,16));
						
						avgSizeDistort.setText(format(stdDev(distortionG[0])));
						avgShapeDistort.setText(format(average(distortionG[1])));
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
		final String p = projectionChooser.getValue();
		
		ProgressBarDialog pBar = new ProgressBarDialog();
		pBar.show();
		new Thread(() -> {
			final double[][][] distortion = calculateDistortion(map(1000,p), p, pBar);
			Image graphic = makeGraphic(distortion);
			Platform.runLater(() -> saveImage(graphic, pBar));
		}).start();
	}
	
	
	private double[][][] calculateDistortion(double[][][] points, String proj) {
		return calculateDistortion(points, proj, null);
	}
	
	
	private double[][][] calculateDistortion(double[][][] points, String proj,
			ProgressBarDialog pBar) { //calculate both kinds of distortion over the given region
		double[][][] output = new double[2][points.length][points[0].length]; //the distortion matrix
		
		for (int y = 0; y < points.length; y ++) {
			for (int x = 0; x < points[y].length; x ++) {
				if (points[y][x] != null) {
					final double[] distortions = getDistortionAt(points[y][x], proj);
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
		
		final double avgArea = average(output[0]); //don't forget to normalize output[0] so the average is zero
		for (int y = 0; y < output[0].length; y ++)
			for (int x = 0; x < output[0][y].length; x ++)
				output[0][y][x] -= avgArea;
		
		return output;
	}
	
	
	private double[] getDistortionAt(double[] s0, String proj) { //calculate both kinds of distortion at the given point
		final double[] output = new double[2];
		final double dx = 1e-5;
		
		final double[] s1 = { s0[0], s0[1]+dx/Math.cos(s0[0]) }; //consider a point slightly to the east
		final double[] s2 = { s0[0]+dx, s0[1] }; //and slightly to the north
		final double[] p0 = vectormaps.MapProjections.project(s0, proj);
		final double[] p1 = vectormaps.MapProjections.project(s1, proj);
		final double[] p2 = vectormaps.MapProjections.project(s2, proj);
		
		final double dA = 
				(p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]);
		output[0] = Math.log(Math.abs(dA/(dx*dx))); //the zeroth output is the size (area) distortion
		if (Math.abs(output[0]) > 25)
			output[0] = Double.NaN; //discard outliers
		
		final double s1ps2 = Math.hypot((p1[0]-p0[0])+(p2[1]-p0[1]), (p1[1]-p0[1])-(p2[0]-p0[0]));
		final double s1ms2 = Math.hypot((p1[0]-p0[0])-(p2[1]-p0[1]), (p1[1]-p0[1])+(p2[0]-p0[0]));
		final double factor = Math.abs((s1ps2-s1ms2)/(s1ps2+s1ms2)); //there's some linear algebra behind this formula. Don't worry about it.
		
		output[1] = Math.asin(1-factor);
		if (output[1] >= Math.PI/2)
			output[1] = Double.NaN;
		
		return output;
	}
	
	
	private Image makeGraphic(double[][][] distortion) {
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
					r = (int)(256*(1-shapeDistort/(Math.PI/2)));
					g = (int)(256*(1-shapeDistort/(Math.PI/2))*(1-Math.min(1, -sizeDistort/3)));
					b = g;
				}
				else { //if dilating
					r = (int)(256*(1-shapeDistort/(Math.PI/2))*(1-Math.min(1, sizeDistort/3)));
					g = r;
					b = (int)(256*(1-shapeDistort/(Math.PI/2)));
				}
				
				final int argb = ((((((0xFF)<<8)+r)<<8)+g)<<8)+b;
				writer.setArgb(x, y, argb);
			}
		}
		return output;
	}
	
	
	private double[][][] map(int size, String proj) { //generate a matrix of coordinates based on a map projection
		int projIdx = 0;
		for (int i = 0; i < PROJ_ARR.length; i ++) {
			if (PROJ_ARR[i].equals(proj)) {
				projIdx = i;
				break;
			}
		}
		
		final int[] dims = {size, (int)(size/DEFA[projIdx])};
		double[][][] output = new double[dims[1]][dims[0]][2]; //the coordinate matrix
		
		for (int y = 0; y < output.length; y ++)
			for (int x = 0; x < output[y].length; x ++)
				output[y][x] = rastermaps.MapProjections.project(x, y, proj, dims); //s0 is this point on the sphere
		
		return output;
	}
	
	
	private double[][][] globe(double dt, String proj) { //generate a matrix of coordinates based on the sphere
		List<double[]> points = new ArrayList<double[]>();
		for (double phi = -Math.PI/2+dt/2; phi < Math.PI/2; phi += dt) { // make sure phi is never exactly +-tau/4
			for (double lam = -Math.PI; lam < Math.PI; lam += dt/Math.cos(phi)) {
				points.add(new double[] {phi, lam});
			}
		}
		return new double[][][] {points.toArray(new double[0][])};
	}
	
	
	public void saveImage(Image img, ProgressBarDialog pBar) { // call from the main thread!
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
	
	
	private final Series<String, Number> histogram(double[][] values,
			double min, double max, int num) {
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
		for (int i = 0; i <= num; i ++)
			output.getData().add(new Data<String, Number>(
					Double.toString(Math.round(10*(i*(max-min)/num+min))/10.),
					(double)hist[i]/tot*100));
		return output;
	}
	
	
	private final double average(double[][] values) { //get the average
		double s = 0, n = 0;
		for (double[] row: values) {
			for (double x: row) {
				if (Double.isFinite(x)) { //ignore NaN values in the average
					s += x;
					n += 1;
				}
			}
		}
		return s/n;
	}
	
	
	private final double stdDev(double[][] values) {
		double s = 0, ss = 0, n = 0;
		for (double[] row: values) {
			for (double x: row) {
				if (Double.isFinite(x)) {
					s += x;
					ss += x*x;
					n += 1;
				}
			}
		}
		return Math.sqrt(ss/n - s*s/n*n);
	}
	
	
	private final String format(double d) {
		if (d < 1000)
			return Double.toString(Math.round(d*1000.)/1000.);
		else
			return "1000+";
	}

}