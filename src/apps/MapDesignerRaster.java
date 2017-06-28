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
import java.lang.reflect.Array;
import java.util.Optional;
import javax.imageio.ImageIO;

import dialogs.MapConfigurationDialog;
import dialogs.ProgressBarDialog;
import javafx.application.Application;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.concurrent.Task;
import javafx.embed.swing.SwingFXUtils;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.MenuButton;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Separator;
import javafx.scene.control.Slider;
import javafx.scene.control.Spinner;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.image.PixelReader;
import javafx.scene.image.WritableImage;
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
 * An application to make raster oblique aspects of map projections
 * 
 * @author Justin Kunimune
 */
public class MapDesignerRaster extends Application {

	private static final int CONT_WIDTH = 300;
	private static final int IMG_WIDTH = 500;
	
	
	private static final KeyCombination ctrlO = new KeyCodeCombination(KeyCode.O, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination ctrlS = new KeyCodeCombination(KeyCode.S, KeyCodeCombination.CONTROL_DOWN);
	
	
	private static final Projection[] PROJ_ARR = { Projection.MERCATOR, Projection.EQUIRECTANGULAR, Projection.HOBODYER,
			Projection.GALL, Projection.STEREOGRAPHIC, Projection.POLAR, Projection.E_A_AZIMUTH,
			Projection.ORTHOGRAPHIC, Projection.GNOMONIC, Projection.LAMBERT_CONIC, Projection.E_D_CONIC,
			Projection.ALBERS, Projection.LEE, Projection.TETRAGRAPH, Projection.AUTHAGRAPH, Projection.SINUSOIDAL,
			Projection.AITOFF, Projection.MOLLWEIDE, Projection.TOBLER, Projection.VAN_DER_GRINTEN, Projection.ROBINSON,
			Projection.WINKEL_TRIPEL, Projection.PEIRCE_QUINCUNCIAL, Projection.GUYOU, Projection.LEMONS,
			Projection.MAGNIFIER, Projection.EXPERIMENT };
	
	private static final String[] AXES = { "Standard", "Transverse", "Center of Mass", "Jerusalem", "Point Nemo",
			"Longest Line", "Longest Line Transverse", "Cylindrical", "Conic", "Tetrahedral", "Quincuncial", "Antipode", "Random" };
	private static final double[][] DEF_ASPECTS = { { 90, 0, 29.9792, 31.7833, 48.8767, -28.5217, -46.4883, -35, -10, 47, 60 },
			{ 0, 0, 31.1344, 35.216, 56.6067, 141.451, 16.5305, -13.6064, 65, -173, -6 },
			{ 0, 0, -32, -35, -45, 161.5, 137, 145, -150, 138, -10 } };
	
	
	private Stage stage;
	private FileChooser inputChooser, saver;
	private Text inputLabel;
	private Button changeInput;
	private ComboBox<Projection> projectionChooser;
	private Text projectionDesc;
	private Slider[] aspectSliders;
	private Spinner<Double>[] aspectSpinners;
	private Button update, saveMap;
	private Image input;
	private ImageView output;
	
	
	
	public static final void main(String[] args) {
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
				new FileChooser.ExtensionFilter("JPG", "*.jpg; *.jpeg; *.jpe; *.jfif"),
				new FileChooser.ExtensionFilter("PNG", "*.png"));
		
		changeInput = new Button("Choose input...");
		changeInput.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				chooseInput();
			}
		});
		changeInput.setTooltip(new Tooltip(
				"Choose the image to determine your map's color scheme"));
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
			}
		});
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
		
		update = new Button("Update Map");
		update.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				updateMap();
			}
		});
		update.setTooltip(new Tooltip(
				"Update the current map with your parameters."));
		update.setDefaultButton(true);
		
		saver = new FileChooser();
		saver.setInitialDirectory(new File("output"));
		saver.setInitialFileName("myMap.jpg");
		saver.setTitle("Save Map");
		saver.getExtensionFilters().addAll(
				new FileChooser.ExtensionFilter("JPG", "*.jpg"),
				new FileChooser.ExtensionFilter("PNG", "*.png"));
		
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
		
		HBox box = new HBox(5, update, saveMap);
		box.setAlignment(Pos.CENTER);
		layout.getChildren().add(box);
		
		output = new ImageView();
		output.setFitWidth(IMG_WIDTH);
		output.setFitHeight(IMG_WIDTH);
		output.setPreserveRatio(true);
		
		new Thread(() -> {
			setInput("basic.jpg", "input/basic.jpg");
			update.fire();
		}).start();
		
		final HBox gui = new HBox(layout, output);
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
			}).start();
		}
	}
	
	
	private void setInput(String name, String filename) {
		changeInput.setDisable(true);
		update.setDisable(true);
		saveMap.setDisable(true);
		
		try {
			input = new Image("file:"+filename);
			inputLabel.setText(name);
		} catch (IllegalArgumentException e) {
			final Alert alert = new Alert(Alert.AlertType.ERROR);
			alert.setHeaderText("File not found!");
			alert.setContentText("Couldn't find "+filename+".");
		}
		
		changeInput.setDisable(false);
		update.setDisable(false);
		saveMap.setDisable(false);
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
		update.setDisable(true);
		new Thread(new Task<Void>() {
			protected Void call() {
				output.setImage(map());
				update.setDisable(false);
				return null;
			}
		}).start();
	}
	
	
	private void startFinalizingMap() {
		Projection p = projectionChooser.getValue();
		Dialog<Thread> dialog =
				new MapConfigurationDialog(p.getAspectRatio(), this);
		Optional<Thread> mapMaker = dialog.showAndWait();
		if (mapMaker.isPresent())	mapMaker.get().start();
		else						return;
	}
	
	
	public void saveImage(Image img, ProgressBarDialog pBar) {	// call from the main thread!
		pBar.close();
		
		final File f = saver.showSaveDialog(stage);
		if (f != null) {
			new Thread(() -> {
				try {
					saveMap.setDisable(true);
					ImageIO.write(SwingFXUtils.fromFXImage(img,null), "jpg", f);
					saveMap.setDisable(false);
				} catch (IOException e) {}
			}).start();
		}
	}
	
	
	public Image map() {
		final double a = projectionChooser.getValue().getAspectRatio();
		return map(IMG_WIDTH, (int)(IMG_WIDTH/a), 1);
	}
	
	public Image map(int outputWidth, int outputHeight, int smoothing) {
		return map(outputWidth,outputHeight,smoothing, null);
	}
	
	public Image map(int outputWidth, int outputHeight, int smoothing,
			ProgressBarDialog pbar) {
		final Projection proj = projectionChooser.getValue();
		final double[] pole = new double[3];
		for (int i = 0; i < 3; i ++)
			pole[i] = Math.toRadians(aspectSliders[i].getValue());
		final PixelReader ref = input.getPixelReader();
		final int[] refDims = {(int)input.getWidth(), (int)input.getHeight()};
		final int[] outDims = {outputWidth, outputHeight};
		
		WritableImage img = new WritableImage(outputWidth, outputHeight);
		
		for (int x = 0; x < outputWidth; x ++) {
			for (int y = 0; y < outputHeight; y ++) {
				int[] colors = new int[smoothing*smoothing];
				int i = 0;
				for (double dx = 0.5/smoothing; dx < 1; dx += 1.0/smoothing) {
					for (double dy = .5/smoothing; dy < 1; dy += 1.0/smoothing) {
						colors[i] = getArgb(x+dx, y+dy,
								proj,pole,ref,refDims,outDims);
						i ++;
					}
				}
				img.getPixelWriter().setArgb(x, y, blend(colors));
			}
			if (pbar != null)	pbar.setProgress((double)(x+1)/outputWidth);
		}
		
		return img;
	}
	
	
	public static int getArgb(double x, double y, Projection proj, double[] pole,
			PixelReader ref, int[] refDims, int[] outDims) {
		final double[] coords = proj.inverse(
				2.*x/outDims[0]-1, 1-2.*y/outDims[1], pole);
		if (coords != null)
			return getColorAt(coords, ref, refDims);
		else
			return 0;
	}
	
	
	public static int getColorAt(double[] coords,
			PixelReader ref, int[] refDims) { // returns the color of any coordinate on earth
		double x = 1/2.0 + coords[1]/(2*Math.PI);
		x = (x - Math.floor(x)) * refDims[0];
		
		double y = refDims[1]/2.0 - coords[0]*refDims[1]/(Math.PI);
		if (y < 0)
			y = 0;
		else if (y >= refDims[1])
			y = refDims[1] - 1;
		
		return ref.getArgb((int)x, (int)y);
	}
	
	
	public static final int blend(int[] colors) {
		int a_tot = 0;
		int r_tot = 0;
		int g_tot = 0;
		int b_tot = 0;
		for (int argb: colors) {
			double a = ((argb >> 24) & 0xFF);
			a_tot += a;
			r_tot += a*((argb >> 16) & 0xFF);
			g_tot += a*((argb >> 8) & 0xFF);
			b_tot += a*((argb >> 0) & 0xFF);
		}
		if (a_tot == 0)	return 0;
		else
			return (a_tot/colors.length<<24) +
					(r_tot/a_tot<<16) + (g_tot/a_tot<<8) + (b_tot/a_tot);
	}

}