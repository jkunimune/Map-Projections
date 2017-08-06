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
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import dialogs.ProgressBarDialog;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Button;
import javafx.scene.control.Separator;
import javafx.scene.layout.HBox;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import maps.Projection;
import util.Math2;
import util.Procedure;

/**
 * An application to make vector oblique aspects of map projections
 * 
 * @author Justin Kunimune
 */
public class MapDesignerVector extends MapApplication {

	public static final void main(String[] args) {
		launch(args);
	}
	
	
	
	private static final FileChooser.ExtensionFilter[] VECTOR_TYPES = {
			new FileChooser.ExtensionFilter("SVG", "*.svg") };
	
	private static final Projection[] PROJ_ARR = { Projection.MERCATOR, Projection.EQUIRECTANGULAR,
			Projection.E_A_CYLIND, Projection.GALL, Projection.STEREOGRAPHIC, Projection.POLAR, Projection.E_A_AZIMUTH,
			Projection.ORTHOGRAPHIC, Projection.GNOMONIC, Projection.LAMBERT_CONIC, Projection.E_D_CONIC,
			Projection.ALBERS, Projection.LEE, Projection.TETRAGRAPH, Projection.SINUSOIDAL, Projection.MOLLWEIDE,
			Projection.TOBLER, Projection.AITOFF, Projection.VAN_DER_GRINTEN, Projection.ROBINSON,
			Projection.WINKEL_TRIPEL, Projection.PEIRCE_QUINCUNCIAL, Projection.GUYOU, Projection.TWO_POINT_EQUIDISTANT, Projection.HAMMER_RETROAZIMUTHAL,
			Projection.MAGNIFIER, Projection.EXPERIMENT, Projection.PSEUDOSTEREOGRAPHIC, Projection.HYPERELLIPOWER,
			Projection.TETRAPOWER, Projection.TETRAFILLET, Projection.TETRACHAMFER };
	
	private static final int DEF_MAX_VTX = 5000;
	private static final int FAST_MAX_VTX = 2000;
	
	
	private Node aspectSelector;
	private Button saveBtn;
	private double[] aspect;
	private List<String> format;
	private List<List<double[]>> input;
	private double minX, maxX, minY, maxY;
	private int numVtx;
	private Canvas viewer;
	
	
	
	public MapDesignerVector() {
		super("Map Designer");
	}
	
	
	
	@Override
	public void start(Stage root) {
		super.start(root);
		new Thread(() -> {
			setInput(new File("input/basic.svg")); //this automatically updates the map//TODO: this should cause the buttons to grey out
		}).start();
	}
	
	
	@Override
	protected Node makeWidgets() {
		this.aspect = new double[3];
		final Node inputSelector = buildInputSelector(VECTOR_TYPES,
				VECTOR_TYPES[0], this::setInput);
		final Node projectionSelector = buildProjectionSelector(PROJ_ARR,
				Procedure.concat(this::updateMap, this::hideAspect));
		this.aspectSelector = buildAspectSelector(this.aspect, this::updateMap);
		final Node parameterSelector = buildParameterSelector(this::updateMap);
		this.saveBtn = buildSaveButton(true, "map", VECTOR_TYPES,
				VECTOR_TYPES[0], ()->true, this::calculateAndSaveMap);
		aspectSelector.managedProperty().bind(aspectSelector.visibleProperty());
		
		final VBox layout = new VBox(5,
				inputSelector, new Separator(), projectionSelector,
				new Separator(), aspectSelector, parameterSelector,
				new Separator(), saveBtn);
		
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(GUI_WIDTH);
		
		viewer = new Canvas(IMG_WIDTH, IMG_WIDTH);
		
		final HBox gui = new HBox(10, layout, viewer);
		gui.setAlignment(Pos.CENTER);
		StackPane.setMargin(gui, new Insets(10));
		
		return gui;
	}
	
	
	private void setInput(File file) {
		saveBtn.setDisable(true);
		
		try {
			input = loadSVG(file.getAbsolutePath());
		} catch (IllegalArgumentException e) {
			showError("Unreadable file!",
					"We couldn't read "+file.getAbsolutePath()+". It may be corrupt or an unreadable format.");
		} catch (IOException e) {
			showError("File not found!",
					"We couldn't find "+file.getAbsolutePath()+".");
		} finally {
			saveBtn.setDisable(false);
		}
		
		updateMap();
	}
	
	
	private List<List<double[]>> loadSVG(String filename) throws IOException { // this method is just awful.
		input = new ArrayList<List<double[]>>();
		format = new ArrayList<String>();
		minX = minY = Integer.MAX_VALUE;
		maxX = maxY = Integer.MIN_VALUE;
		numVtx = 0;
		
		BufferedReader in = new BufferedReader(new FileReader(filename));
		String formatStuff = "";
		int c = in.read();
		
		do {
			formatStuff += (char)c;
			if (formatStuff.length() >= 4 &&
					formatStuff.substring(formatStuff.length()-4).equals(" d=\"")) {
				format.add(formatStuff);
				formatStuff = "";
				List<double[]> currentShape = new ArrayList<double[]>();
				c = in.read();
				do {
					if (c == 'Z') {
						input.add(currentShape);
						currentShape = new ArrayList<double[]>();
						format.add("");
					}
					if (isDigit((char) c)) {
						String num = Character.toString((char)c);
						while (isDigit((char) (c = in.read())))
							num += (char) c;
						double x = Double.parseDouble(num);
						if (x < minX)	minX = x;
						if (x > maxX)	maxX = x;
						c = in.read();
						num = Character.toString((char)c);
						while (isDigit((char) (c = in.read())))
							num += (char) c;
						double y = Double.parseDouble(num);
						if (y < minY)	minY = y;
						if (y > maxY)	maxY = y;
						currentShape.add(new double[] {x,y});
						numVtx ++;
					}
					else {
						c = in.read();
					}
				} while (c != '"');
			}
			else {
				c = in.read();
			}
		} while (c >= 0);
		
		format.add(formatStuff);
		in.close();
		
		for (List<double[]> curve: input) {
			for (double[] coords: curve) {
				final double[] radCoords = convCoordsToMathy(coords);
				coords[0] = radCoords[0];
				coords[1] = radCoords[1];
			}
		}
		return input;
	}
	
	
	private void hideAspect() {
		aspectSelector.setVisible(this.getProjection().hasAspect());
	}
	
	
	private void updateMap() {
		int maxVtx = this.getParamsChanging() ? FAST_MAX_VTX : DEF_MAX_VTX;
		final List<List<double[]>> transformed = this.getProjection().transform(
				input, numVtx/maxVtx+1, this.getParams(), aspect.clone(), null);
		Platform.runLater(() -> drawImage(transformed, viewer));
	}
	
	
	private void calculateAndSaveMap(File file, ProgressBarDialog pBar) {
		final List<List<double[]>> transformed = this.getProjection().transform(
				input, 1, this.getParams(), aspect.clone(), pBar::setProgress); //calculate
		saveToSVG(transformed, file, pBar); //save
	}
	
	
	private void drawImage(List<List<double[]>> img, Canvas c) {
		GraphicsContext g = c.getGraphicsContext2D();
		g.clearRect(0, 0, c.getWidth(), c.getHeight());
		for (List<double[]> closedCurve: img) {
			final double[] start = convCoordsToImg(closedCurve.get(0));
			g.beginPath();
			g.moveTo(start[0], start[1]);
			for (int i = 1; i < closedCurve.size()+1; i ++) {
				double[] p0 = convCoordsToImg(closedCurve.get(i-1));
				double[] p1 = convCoordsToImg(closedCurve.get(i%closedCurve.size()));
				if (Math.hypot(p1[0]-p0[0], p1[1]-p0[1]) >= c.getWidth()/8) {
					g.stroke();
					g.beginPath();
					g.moveTo(p1[0], p1[1]);
				}
				else
					g.lineTo(p1[0], p1[1]);
			}
			g.stroke();
		}
	}
	
	
	private void saveToSVG(List<List<double[]>> curves, File file,
			ProgressBarDialog pBar) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(file));
			
			for (int i = 0; i < curves.size(); i ++) {
				out.write(format.get(i));
				String curveCode = "M";
				for (double[] r: curves.get(i)) {
					final double[] p = convCoordsToSVG(r);
					curveCode += p[0]+" "+p[1]+"L";
				}
				out.write(curveCode.substring(0,curveCode.length()-2)+"Z");
				pBar.setProgress((double)i/curves.size());
			}
			out.write(format.get(curves.size()));
			out.close();
			
		} catch (IOException e) {
			showError("Failure!",
					"Could not access "+file.getAbsolutePath()+". It's possible that another program has it open.");
		}
		saveBtn.setDisable(false);
	}
	
	
	private double[] convCoordsToMathy(double... coords) { // changes svg coordinates to radians
		final double NORTHMOST = 1.459095; //TODO: use image height and width
		final double SOUTHMOST = -1.4868809;
		final double EASTMOST = -Math.PI;
		final double WESTMOST = Math.PI;
		return new double[] {Math2.linInterp(coords[1], minY,maxY, SOUTHMOST,NORTHMOST),
				Math2.linInterp(coords[0], minX,maxX, EASTMOST,WESTMOST)};
	}
	
	
	private double[] convCoordsToImg(double... coords) { // changes [-1,1] coordinates to image coordinates
		return new double[] {Math2.linInterp(coords[0], -Math.PI,Math.PI, 0,IMG_WIDTH),
				Math2.linInterp(coords[1], -Math.PI,Math.PI, IMG_WIDTH,0)};
	}
	
	
	private double[] convCoordsToSVG(double... coords) { // changes image coordinates to svg coordinates
		return new double[] {Math2.linInterp(coords[0], -Math.PI,Math.PI, -180,180), //TODO: use image height and width
				Math2.linInterp(coords[1], -Math.PI,Math.PI, -180,180)};
	}
	
	
	private static final boolean isDigit(char c) {
		return c >= '-' && c <= '9';
	}

}
