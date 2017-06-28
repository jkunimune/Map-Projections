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
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import dialogs.ProgressBarDialog;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.concurrent.Task;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.MenuButton;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Separator;
import javafx.scene.control.Slider;
import javafx.scene.control.Spinner;
import javafx.scene.control.Tooltip;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
import javafx.scene.input.KeyEvent;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import maps.Projection;

/**
 * An application to make vector oblique aspects of map projections
 * 
 * @author Justin Kunimune
 */
public class MapDesignerVector extends Application { //TODO: inheritance

	private static final int CONT_WIDTH = 300;
	private static final int IMG_WIDTH = 500;
	
	
	private static final KeyCombination ctrlO = new KeyCodeCombination(KeyCode.O, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination ctrlS = new KeyCodeCombination(KeyCode.S, KeyCodeCombination.CONTROL_DOWN);
	
	
	private static final Projection[] PROJ_ARR = { Projection.MERCATOR, Projection.EQUIRECTANGULAR, Projection.HOBODYER,
			Projection.GALL, Projection.STEREOGRAPHIC, Projection.POLAR, Projection.E_A_AZIMUTH,
			Projection.ORTHOGRAPHIC, Projection.GNOMONIC, Projection.LAMBERT_CONIC, Projection.E_D_CONIC,
			Projection.ALBERS, Projection.LEE, Projection.TETRAGRAPH, Projection.SINUSOIDAL, Projection.MOLLWEIDE,
			Projection.TOBLER, Projection.AITOFF, Projection.VAN_DER_GRINTEN, Projection.ROBINSON,
			Projection.WINKEL_TRIPEL, Projection.PEIRCE_QUINCUNCIAL, Projection.GUYOU, Projection.MAGNIFIER,
			Projection.EXPERIMENT, Projection.HYPERELLIPOWER, Projection.TETRAPOWER, Projection.TETRAFILLET };
	
	private static final String[] AXES = { "Standard", "Transverse", "Center of Mass", "Jerusalem", "Point Nemo",
			"Longest Line", "Longest Line Transverse", "Cylindrical", "Conical", "Tetrahedral", "Quincuncial", "Antipode", "Random" };
	private static final double[][] DEF_ASPECTS = { { 90, 0, 29.9792, 31.7833, 48.8767, -28.5217, -46.4883, -35, -10, 47, 60 },
													{ 0, 0, 31.1344, 35.216, 56.6067, 141.451, 16.5305, -13.6064, 65, -173, -6 },
													{ 0, 0, -32, -35, -45, 161.5, 137, 145, -150, 138, -10 } };
	private static final int DEF_MAX_VTX = 5000;
	
	
	private Stage stage;
	private FileChooser inputChooser, saver;
	private Text inputLabel;
	private Button changeInput;
	private ComboBox<Projection> projectionChooser;
	private Text projectionDesc;
	private Slider[] aspectSliders;
	private Spinner<Double>[] aspectSpinners;
	private Button saveMap;
	private List<String> format;
	private List<List<double[]>> input;
	private double minX, maxX, minY, maxY;
	private int numVtx;
	private Canvas viewer;
	
	
	
	public static void main(String[] args) {
		launch(args);
	}
	
	
	@SuppressWarnings("unchecked")
	@Override
	public void start(Stage primaryStage) {
		stage = primaryStage;
		stage.setTitle("Map Designer");
		
		final VBox layout = new VBox();
		layout.setSpacing(5);
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(CONT_WIDTH);
		
		Label lbl = new Label("Current input:");
		inputLabel = new Text("None");
		layout.getChildren().add(new HBox(3, lbl, inputLabel));
		
		inputChooser = new FileChooser();
		inputChooser.setInitialDirectory(new File("input"));
		inputChooser.setTitle("Choose an input map");
		inputChooser.getExtensionFilters().addAll(
				new FileChooser.ExtensionFilter("SVG", "*.svg"));
		
		changeInput = new Button("Choose input...");
		changeInput.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				chooseInput();
			}
		});
		changeInput.setTooltip(new Tooltip(
				"Change the input image"));
		stage.addEventHandler(KeyEvent.KEY_RELEASED, new EventHandler<KeyEvent>() {	// ctrl-O opens
			public void handle(KeyEvent event) {
				if (ctrlO.match(event))	changeInput.fire();
			}
		});
		layout.getChildren().add(changeInput);
		
		layout.getChildren().add(new Separator());
		
		lbl = new Label("Projection:");
		ObservableList<Projection> items = FXCollections.observableArrayList(PROJ_ARR);
		projectionChooser = new ComboBox<Projection>(items);
		projectionChooser.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				projectionDesc.setText(projectionChooser.getValue().getDescription());
				updateMap(); //TODO: how cool would it be to animate all of the transitions?
			}
		});
		projectionChooser.setPrefWidth(210);
		projectionChooser.setValue(Projection.MERCATOR);
		layout.getChildren().add(new HBox(3, lbl, projectionChooser));
		
		projectionDesc = new Text(projectionChooser.getValue().getDescription());
		projectionDesc.setWrappingWidth(CONT_WIDTH);
		layout.getChildren().add(projectionDesc);
		
		layout.getChildren().add(new Separator());
		
		final MenuButton defAxes = new MenuButton("Aspect Presets");
		for (String preset: AXES) {
			MenuItem m = new MenuItem(preset);
			m.setOnAction(new EventHandler<ActionEvent>() {
				public void handle(ActionEvent event) {
					setAxisByPreset(((MenuItem) event.getSource()).getText());
					updateMap();
				}
			});
			defAxes.getItems().add(m);
		}
		defAxes.setTooltip(new Tooltip(
				"Set the aspect sliders based on a preset"));
		layout.getChildren().add(defAxes);
		
		aspectSliders = new Slider[] {
				new Slider(-90, 90, 90),
				new Slider(-180,180,0),
				new Slider(-180,180,0)
		};
		final Spinner<Double> latSpinner = new Spinner<Double>(-90, 90, 90.0);
		aspectSpinners = (Spinner<Double>[]) Array.newInstance(latSpinner.getClass(), 3);
		aspectSpinners[0] = latSpinner;
		aspectSpinners[1] = new Spinner<Double>(-180, 180, 0.0);
		aspectSpinners[2] = new Spinner<Double>(-180, 180, 0.0);
		
		GridPane grid = new GridPane();
		grid.addRow(0, new Text("Latitude:"), aspectSliders[0], aspectSpinners[0]);
		grid.addRow(1, new Text("Longitude:"), aspectSliders[1], aspectSpinners[1]);
		grid.addRow(2, new Text("Orientation:"), aspectSliders[2], aspectSpinners[2]);
		
		for (int i = 0; i < 3; i ++) {
			final Slider sld = aspectSliders[i];
			final Spinner<Double> spn = aspectSpinners[i];
			GridPane.setHgrow(sld, Priority.ALWAYS);
			sld.setTooltip(new Tooltip("Change the aspect of the map"));
			sld.valueChangingProperty().addListener(
					(observable, then, now) -> {
						if (!now)	updateMap();
						if (spn.getValue() != sld.getValue())
							spn.getEditor().textProperty().set(Double.toString(sld.getValue()));
					});
			
			spn.setTooltip(new Tooltip("Change the aspect of the map"));
			spn.setPrefWidth(100);
			spn.setEditable(true);
			spn.getEditor().textProperty().addListener((ov, pv, nv) -> {	// link the Spinners
				if (spn.getEditor().textProperty().isEmpty().get())	return;
				try {
					Double.parseDouble(nv);
					spn.increment(0);	// forces the spinner to commit its value
					if (spn.getValue() != sld.getValue())
						sld.setValue(spn.getValue());
				} catch (NumberFormatException e) {
					spn.getEditor().textProperty().set(pv); //yeah, this is all pretty jank. JavaFX spinners are just weird by default
				}
			});
			spn.setOnKeyPressed((event) -> {
				System.out.println(event);
				if (event.getCode() == KeyCode.ENTER)
					updateMap();
			});
		}
		layout.getChildren().add(grid);
		
		layout.getChildren().add(new Separator());
		
		saver = new FileChooser();
		saver.setInitialDirectory(new File("output"));
		saver.setInitialFileName("myMap.svg");
		saver.setTitle("Save Map");
		saver.getExtensionFilters().addAll(
				new FileChooser.ExtensionFilter("SVG", "*.svg"));
		
		saveMap = new Button("Save Map...");
		saveMap.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				startFinalizingMap();
			}
		});
		saveMap.setTooltip(new Tooltip("Save the map with current settings."));
		stage.addEventHandler(KeyEvent.KEY_RELEASED, new EventHandler<KeyEvent>() {	// ctrl-S saves
			public void handle(KeyEvent event) {
				if (ctrlS.match(event))	saveMap.fire();
			}
		});
		
		HBox box = new HBox(5, saveMap);
		box.setAlignment(Pos.CENTER);
		layout.getChildren().add(box);
		
		viewer = new Canvas(IMG_WIDTH, IMG_WIDTH);
		
		new Thread(() -> {
			setInput("basic.svg", "input/basic.svg");
			updateMap();
		}).start();
		
		final HBox gui = new HBox(layout, viewer);
		gui.setAlignment(Pos.CENTER);
		gui.setSpacing(10);
		StackPane.setMargin(gui, new Insets(10));
		stage.setScene(new Scene(new StackPane(gui)));
		stage.show();
	}
	
	
	private void chooseInput() {
		final File f = inputChooser.showOpenDialog(stage);
		if (f != null) {
			new Thread(() -> {
				setInput(f.getName(), f.getAbsolutePath());
				updateMap();
			}).start();
		}
	}
	
	
	private void setInput(String name, String filename) {
		changeInput.setDisable(true);
		saveMap.setDisable(true);
		
		try {
			input = loadSVG(filename);
			inputLabel.setText(name);
		} catch (IllegalArgumentException e) {
			final Alert alert = new Alert(Alert.AlertType.ERROR);
			alert.setHeaderText("Unreadable file!");
			alert.setContentText("We couldn't read "+filename+". It may be corrupt or an unreadable format.");
		} catch (IOException e) {
			final Alert alert = new Alert(Alert.AlertType.ERROR);
			alert.setHeaderText("File not found!");
			alert.setContentText("Couldn't find "+filename+".");
		}
		
		changeInput.setDisable(false);
		saveMap.setDisable(false);
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
						currentShape.add(new double[] {x, y});
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
		return input;
	}
	
	
	private void setAxisByPreset(String preset) {
		if (preset.equals("Antipode")) {
			for (Slider s: aspectSliders)
				s.setValue(-s.getValue());
			aspectSliders[1].setValue((-aspectSliders[1].getValue()+360)%360-180);
			return;
		}
		if (preset.equals("Random")) {
			aspectSliders[0].setValue(Math.toDegrees(Math.asin(Math.random()*2-1)));
			aspectSliders[1].setValue(Math.random()*360-180);
			aspectSliders[2].setValue(Math.random()*360-180);
			return;
		}
		for (int i = 0; i < AXES.length; i ++) {
			if (AXES[i].equals(preset)) {
				for (int j = 0; j < 3; j ++)
					aspectSliders[j].setValue(DEF_ASPECTS[j][i]);
				break;
			}
		}
	}
	
	
	private void updateMap() {
		updateMap(false);
	}
	
	private void updateMap(boolean fast) {
		new Thread(new Task<Void>() {
			protected Void call() {
				try {
					drawImage(map(fast), viewer);
				} catch (Exception e) {
					e.printStackTrace();
				}
				return null;
			}
		}).start();
	}
	
	
	private void startFinalizingMap() {
		saveMap.setDisable(true);
		
		int p = 0;
		while (p < PROJ_ARR.length &&
				!PROJ_ARR[p].equals(projectionChooser.getValue()))
			p ++;
		
		ProgressBarDialog pBar = new ProgressBarDialog();
		pBar.show();
		new Thread(() -> {
			List<List<double[]>> map = map(0, pBar);
			Platform.runLater(() -> saveToSVG(map, pBar));
		}).start();
	}
	
	
	private void drawImage(List<List<double[]>> img, Canvas c) {
		GraphicsContext g = c.getGraphicsContext2D();
		g.clearRect(0, 0, c.getWidth(), c.getHeight());
		for (List<double[]> closedCurve: img) {
			g.beginPath();
			for (int i = 0; i < closedCurve.size(); i ++) {
				double[] p = closedCurve.get(i);
				if (i == 0)	g.moveTo(p[0], p[1]);
				else		g.lineTo(p[0], p[1]);
			}
			g.closePath();
			g.stroke();
		}
	}
	
	
	private void saveToSVG(List<List<double[]>> curves, ProgressBarDialog pBar) {	// call from the main thread!
		pBar.close();
		final File f = saver.showSaveDialog(stage);
		if (f == null) {
			saveMap.setDisable(false);
			return;
		}
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(f));
			
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
			Alert alert = new Alert(Alert.AlertType.ERROR);
			alert.setHeaderText("File not found!");
			alert.setContentText("Could not find "+f+".");
		}
		saveMap.setDisable(false);
	}
	
	
	private List<List<double[]>> map(boolean fast) {
		if (fast)
			return map(DEF_MAX_VTX/10, null);
		else
			return map(DEF_MAX_VTX, null);
	}
	
	private List<List<double[]>> map(int maxVtx,
			ProgressBarDialog pbar) {
		final Projection proj = projectionChooser.getValue();
		final double[] pole = new double[3];
		for (int i = 0; i < 3; i ++)
			pole[i] = Math.toRadians(aspectSliders[i].getValue());
		
		int step = maxVtx==0 ? 1 : numVtx/maxVtx+1;
		List<List<double[]>> output = new LinkedList<List<double[]>>();
		
		int i = 0;
		for (List<double[]> curve0: input) {
			if (curve0.size() < step*3)	continue;
			
			List<double[]> curve1 = new ArrayList<double[]>(curve0.size()/step);
			for (int j = 0; j < curve0.size(); j += step) {
				double[] radCoords = convCoordsToMathy(curve0.get(j));
				double[] plnCoords = proj.project(radCoords[0],radCoords[1],pole);
				curve1.add(convCoordsToImg(plnCoords));
			}
			output.add(curve1);
			
			if (pbar != null) {
				i ++;
				pbar.setProgress((double)i/input.size());
			}
		}
		
		return output;
	}
	
	
	private double[] convCoordsToMathy(double[] coords) { // changes svg coordinates to radians
		final double NORTHMOST = 1.459095;
		final double SOUTHMOST = -1.4868809;
		final double EASTMOST = -Math.PI;
		final double WESTMOST = Math.PI;
		return new double[] {linMap(coords[1], minY,maxY, SOUTHMOST,NORTHMOST),
				linMap(coords[0], minX,maxX, EASTMOST,WESTMOST)};
	}
	
	
	private double[] convCoordsToImg(double[] coords) { // changes [-1,1] coordinates to image coordinates
		return new double[] {linMap(coords[0], -Math.PI,Math.PI, 0,IMG_WIDTH),
				linMap(coords[1], -Math.PI,Math.PI, IMG_WIDTH,0)};
	}
	
	
	private double[] convCoordsToSVG(double[] coords) { // changes [-1,1] coordinates to image coordinates
		return new double[] {linMap(coords[0], 0,IMG_WIDTH, -100,100),
				linMap(coords[1], IMG_WIDTH,0, -100,100)};
	}
	
	
	private static final double linMap(double x, double a0, double a1, double b0, double b1) {
		return (x-a0)*(b1-b0)/(a1-a0) + b0;
	}
	
	
	private static final boolean isDigit(char c) {
		return c >= '-' && c <= '9';
	}

}
