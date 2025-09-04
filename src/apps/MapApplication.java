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

import dialogs.ProgressDialog;
import dialogs.ProjectionSelectionDialog;
import image.SavableImage;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.property.ReadOnlyBooleanProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.concurrent.Task;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Control;
import javafx.scene.control.Label;
import javafx.scene.control.MenuButton;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Slider;
import javafx.scene.control.Spinner;
import javafx.scene.control.SpinnerValueFactory;
import javafx.scene.control.SpinnerValueFactory.DoubleSpinnerValueFactory;
import javafx.scene.control.Tooltip;
import javafx.scene.control.cell.ComboBoxListCell;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
import javafx.scene.input.KeyEvent;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.util.converter.DoubleStringConverter;
import maps.ArbitraryPseudocylindrical;
import maps.Azimuthal;
import maps.Conic;
import maps.Cylindrical;
import maps.Danseiji;
import maps.Elastic;
import maps.EqualEarth;
import maps.Gyorffy;
import maps.HammerRetroazimuthal;
import maps.Lenticular;
import maps.Sargent;
import maps.Misc;
import maps.Octahedral;
import maps.Polyhedral;
import maps.Projection;
import maps.Pseudocylindrical;
import maps.Snyder;
import maps.Tobler;
import maps.WinkelTripel;
import utils.Flag;
import utils.MutableDouble;
import utils.Procedure;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.function.BooleanSupplier;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;

import static java.lang.Math.asin;
import static java.lang.Math.random;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;
import static utils.Math2.round;


/**
 * A base class for all GUI applications that deal with maps
 * 
 * @author jkunimune
 */
public abstract class MapApplication extends Application {

	protected static final int GUI_WIDTH = 300;
	protected static final int IMG_SIZE = 480;
	protected static final int V_SPACE = 7;
	protected static final int H_SPACE = 6;
	protected static final int MARGIN = 15;
	protected static final int SPINNER_WIDTH = 90;
	
	private static final KeyCombination CTRL_O = new KeyCodeCombination(KeyCode.O, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination CTRL_S = new KeyCodeCombination(KeyCode.S, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination CTRL_ENTER = new KeyCodeCombination(KeyCode.ENTER, KeyCodeCombination.CONTROL_DOWN);
	
	
	public static final Projection[] FEATURED_PROJECTIONS = { Cylindrical.MERCATOR, // the set of featured projections for the ComboBox
			Cylindrical.EQUIRECTANGULAR, Cylindrical.EQUAL_AREA, Cylindrical.MILLER,
			Azimuthal.STEREOGRAPHIC, Azimuthal.EQUIDISTANT, Azimuthal.EQUAL_AREA, Azimuthal.GNOMONIC,
			Azimuthal.PERSPECTIVE, Conic.LAMBERT, Conic.EQUIDISTANT, Conic.ALBERS,
			Pseudocylindrical.SINUSOIDAL, Pseudocylindrical.MOLLWEIDE, Tobler.TOBLER,
			Lenticular.VAN_DER_GRINTEN, ArbitraryPseudocylindrical.ROBINSON,
			WinkelTripel.WINKEL_TRIPEL, EqualEarth.EQUAL_EARTH, Octahedral.KEYES_SIMPLIFIED,
			Polyhedral.LEE_TETRAHEDRAL_RECTANGULAR, Misc.PEIRCE_QUINCUNCIAL, Misc.LEMONS };

	public static final String[] PROJECTION_CATEGORIES = { "Cylindrical", "Azimuthal", "Conic",
			"Tetrahedral", "Polyhedral", "Pseudocylindrical", "Lenticular", "Mesh-based", "Other" }; //the overarching categories by which I organise my projections
	public static final Projection[][] ALL_PROJECTIONS = { // every projection I have programmed
			{ Cylindrical.EQUAL_AREA, Cylindrical.EQUIRECTANGULAR, Cylindrical.GALL_ORTHOGRAPHIC,
			  Cylindrical.GALL_STEREOGRAPHIC, Cylindrical.HOBO_DYER, Cylindrical.LAMBERT,
			  Cylindrical.MERCATOR, Cylindrical.MILLER, Cylindrical.CENTRAL,
			  Cylindrical.PLATE_CARREE },
			{ Azimuthal.EQUAL_AREA, Azimuthal.EQUIDISTANT, Azimuthal.GNOMONIC, Azimuthal.MAGNIFIER,
			  Azimuthal.ORTHOGRAPHIC, Azimuthal.PERSPECTIVE, Azimuthal.STEREOGRAPHIC },
			{ Conic.ALBERS, Misc.BRAUN_CONIC, Conic.LAMBERT, Conic.EQUIDISTANT },
			{ Polyhedral.AUTHAGRAPH, Polyhedral.AUTHAGRAPH_ALT,
			  Polyhedral.LEE_TETRAHEDRAL_RECTANGULAR, Polyhedral.LEE_TETRAHEDRAL_TRIANGULAR,
			  Polyhedral.VAN_LEEUWEN, Polyhedral.ACTUAUTHAGRAPH, Polyhedral.TETRAGRAPH },
			{ Octahedral.CONFORMAL_CAHILL_BUTTERFLY, Octahedral.CONFORMAL_CAHILL_OCTANT,
			  Octahedral.CAHILL_CONCIALDI, Octahedral.KEYES_STANDARD, Octahedral.KEYES_BUTTERFLY,
			  Octahedral.KEYES_SIMPLIFIED, Octahedral.KEYES_OCTANT, Polyhedral.DYMAXION,
			  Octahedral.WATERMAN_BUTTERFLY, Octahedral.WATERMAN_SIMPLIFIED,
			  Octahedral.WATERMAN_OCTANT },
			{ Pseudocylindrical.ECKERT_IV, EqualEarth.EQUAL_EARTH,
			  Pseudocylindrical.HOMOLOSINE_INTERRUPTED, Pseudocylindrical.HOMOLOSINE,
			  Gyorffy.B, Gyorffy.D, Pseudocylindrical.KAVRAYSKIY_VII, Pseudocylindrical.MOLLWEIDE,
			  ArbitraryPseudocylindrical.ROBINSON, Pseudocylindrical.SINUSOIDAL,
			  Tobler.TOBLER, Pseudocylindrical.WAGNER_II, Pseudocylindrical.WAGNER_V },
			{ Lenticular.AITOFF, Lenticular.POLYCONIC, Lenticular.EISENLOHR, Lenticular.CHINA,
			  Gyorffy.E, Gyorffy.F, Lenticular.HAMMER, Lenticular.LAGRANGE, Lenticular.STREBE_95,
			  Lenticular.VAN_DER_GRINTEN, Lenticular.WAGNER_VIII, WinkelTripel.WINKEL_TRIPEL },
			{ Elastic.ELASTIC_I, Elastic.ELASTIC_II, Elastic.ELASTIC_III,
			  Sargent.LIQUID_EARTH, Sargent.SOLID_EARTH,
			  Danseiji.DANSEIJI_N, Danseiji.DANSEIJI_I, Danseiji.DANSEIJI_II, Danseiji.DANSEIJI_III,
			  Danseiji.DANSEIJI_IV, Danseiji.DANSEIJI_V, Danseiji.DANSEIJI_VI },
			{ Misc.BONNE, Misc.CASSINI, Snyder.GS50, Misc.GUYOU, HammerRetroazimuthal.FULL,
			  HammerRetroazimuthal.FRONT, HammerRetroazimuthal.BACK, Misc.LEMONS,
			  Misc.PEIRCE_QUINCUNCIAL, Misc.SPIRAL, Misc.T_SHIRT, Misc.TWO_POINT_EQUIDISTANT,
			  Misc.FLAT_EARTH },
	};
	
	private static final String[] ASPECT_NAMES = {
			"Normal", "Transverse", "Atlantis", "Jerusalem", "Point Nemo", "Longest Line",
			"Cylindrical", "Tetrahedral", "Antipode", "Random" };
	private static final double[][] ASPECT_VALUES = { //the aspect presets (in degrees)
			{ 90.,  0.,  -4., 31.78, 48.88, -28.52,  35.,   57. },
			{  0., 90.,  65., 35.22, 56.61, 141.45, 166.5,-175.5 },
			{  0.,-90.,-147.,-35.,  -45.,    30.,  -145.,  154. } };
	
	
	private final Map<ButtonType, Button> buttons = new HashMap<ButtonType, Button>();
	
	private final String name;
	private Stage root;
	private ComboBox<Projection> projectionChooser;
	private GridPane paramGrid;
	private Text[] paramLabels;
	private Slider[] paramSliders;
	private Spinner<Double>[] paramSpinners;
	private double[] currentParams;
	
	private final Flag isChanging = new Flag();
	private final Flag suppressListeners = new Flag(); //a flag to prevent events from triggering projection setting
	
	
	
	public MapApplication(String name) {
		super();
		this.name = name;
	}
	
	
	@Override
	public void start(Stage root) {
		this.root = root;
		this.root.setTitle(this.name);
		this.root.setScene(new Scene(new StackPane(makeWidgets())));
		this.root.show();
		
		this.suppressListeners.set();
		this.projectionChooser.setValue(FEATURED_PROJECTIONS[0]);
		this.suppressListeners.clear();
	}
	
	
	protected abstract Node makeWidgets();
	
	
	/**
	 * Create a set of widgets to select an input image.
	 * @param allowedExtensions - The List of possible file types for the input
	 * @param inputSetter - The method to be called when an input file is loaded.
	 * 		This method will be called from a nongui thread.
	 * @return The full formatted Region.
	 */
	protected Region buildInputSelector(
			FileChooser.ExtensionFilter[] allowedExtensions,
			FileChooser.ExtensionFilter defaultExtension,
			Function<File, Task<Void>> inputSetter) {
		final Label label = new Label("Current input:");
		final Text inputLabel = new Text("Basic"+defaultExtension.getExtensions().get(0).substring(1)); //this line kind of cheats, since it assumes the first input will be called "Basic", but I couldn't figure out a good way for this to update when the subclass programmatically sets the input
		
		final FileChooser inputChooser = new FileChooser();
		inputChooser.setInitialDirectory(new File("input"));
		inputChooser.setTitle("Choose an input map");
		inputChooser.getExtensionFilters().addAll(allowedExtensions);
		inputChooser.setSelectedExtensionFilter(defaultExtension);
		
		final Button loadButton = new Button("Choose input\u2026");
		loadButton.setTooltip(new Tooltip(
				"Choose the image to determine the style of your map"));
		loadButton.setOnAction((event) -> {
				File file;
				try {
					file = inputChooser.showOpenDialog(root); //have the user choose a file
				} catch (IllegalArgumentException e) {
					inputChooser.setInitialDirectory(new File(".")); //use the current directory if ./input/ doesn't exist for some reason
					file = inputChooser.showOpenDialog(root);
				}
				if (file != null) {
					inputChooser.setInitialDirectory(file.getParentFile()); //remember this directory for next time
					Task<Void> loadTask = inputSetter.apply(file);
					disableWhile(loadTask.runningProperty(), ButtonType.LOAD_INPUT,
							ButtonType.UPDATE_MAP, ButtonType.SAVE_MAP, ButtonType.SAVE_GRAPH);
					new Thread(loadTask).start(); //start loading
					inputLabel.setText(file.getName());
				}
			});
		loadButton.setTooltip(new Tooltip(
				"Change the input image"));
		root.addEventHandler(KeyEvent.KEY_PRESSED, (event) -> {	// ctrl-O opens
				if (CTRL_O.match(event)) {
					loadButton.requestFocus();
					loadButton.fire();
				}
			});
		
		buttons.put(ButtonType.LOAD_INPUT, loadButton);
		VBox output = new VBox(V_SPACE, new HBox(H_SPACE, label, inputLabel), loadButton);
		output.setAlignment(Pos.CENTER);
		return output;
	}
	
	
	/**
	 * Create a set of widgets to choose a Projection.
	 * @param projectionSetter - The default projection, before the user chooses anything.
	 * @return The full formatted Region.
	 */
	protected Region buildProjectionSelector(Procedure projectionSetter) {
		final Label label = new Label("Projection:");
		projectionChooser =
				new ComboBox<Projection>(FXCollections.observableArrayList(FEATURED_PROJECTIONS));
		projectionChooser.getItems().add(Projection.NULL_PROJECTION);
		
		final Text description = new Text();
		description.setWrappingWidth(GUI_WIDTH);
		
		projectionChooser.valueProperty().addListener((observable, old, now) -> {
				projectionChooser.setButtonCell(new ComboBoxListCell<Projection>()); //This makes it properly display values not in the featured list
				
				final boolean suppressedListeners = suppressListeners.isSet(); //save this value, because revealParameters()...
				if (projectionChooser.getValue() == Projection.NULL_PROJECTION) {
					chooseProjectionFromExpandedList(old); //<aside>NULL_PROJECTION is the "More" button. It triggers the expanded list</aside>
				}
				else {
					description.setText(projectionChooser.getValue().getDescription());
					revealParameters(projectionChooser.getValue()); //...clears suppressListeners. That's fine,
					if (!suppressedListeners) //because suppressListeners is only needed here for that one case.
						projectionSetter.execute();
				}
			});
		
		HBox comboRow = new HBox(H_SPACE, label, projectionChooser);
		comboRow.setAlignment(Pos.CENTER_LEFT);
		return new VBox(V_SPACE, comboRow, description);
	}
	
	
	/**
	 * Create a set of widgets to choose an aspect either from a preset or numbers.
	 * Also bind aspectArr to the sliders.
	 * @return The full formatted Region.
	 */
	protected Region buildAspectSelector(double[] aspectArr, Procedure aspectSetter) {
		final MenuButton presetChooser = new MenuButton("Aspect Presets");
		presetChooser.setTooltip(new Tooltip(
				"Set the aspect sliders based on a preset"));
		
		final String[] labels = { "Latitude:", "Longitude:", "Ctr. Meridian:" };
		final Slider[] sliders = new Slider[] {
				new Slider(-90, 90, 0.),
				new Slider(-180, 180, 0.),
				new Slider(-180, 180, 0.) };
		final Spinner<Double> spin0 = new Spinner<Double>(-90, 90, 0.); //yes, this is awkward. Java gets weird about arrays with generic types
		@SuppressWarnings("unchecked")
		final Spinner<Double>[] spinners = (Spinner<Double>[]) Array.newInstance(spin0.getClass(), 3);
		spinners[0] = spin0;
		spinners[1] = new Spinner<Double>(-180, 180, 0.);
		spinners[2] = new Spinner<Double>(-180, 180, 0.);
		
		for (int i = 0; i < 3; i ++) {
			aspectArr[i] = toRadians(sliders[i].getValue());
			link(sliders[i], spinners[i], i, aspectArr, Math::toRadians,
					aspectSetter, isChanging, suppressListeners);
		}
		setAspectByPreset("Normal", sliders);
		
		for (String preset: ASPECT_NAMES) {
			MenuItem m = new MenuItem(preset);
			m.setOnAction((event) -> {
					setAspectByPreset(((MenuItem) event.getSource()).getText(),
							sliders);
					for (int i = 0; i < 3; i ++)
						aspectArr[i] = toRadians(sliders[i].getValue());
					if (!suppressListeners.isSet())
						aspectSetter.execute();
				});
			presetChooser.getItems().add(m);
		}
		
		final GridPane grid = new GridPane();
		grid.setVgap(V_SPACE);
		grid.setHgap(H_SPACE);
		grid.getColumnConstraints().addAll(
				new ColumnConstraints(Control.USE_COMPUTED_SIZE),
				new ColumnConstraints(), new ColumnConstraints(SPINNER_WIDTH));
		for (int i = 0; i < 3; i ++) {
			GridPane.setHgrow(sliders[i], Priority.ALWAYS);
			sliders[i].setTooltip(new Tooltip("Change the aspect of the map"));
			spinners[i].setTooltip(new Tooltip("Change the aspect of the map"));
			spinners[i].setEditable(true);
			grid.addRow(i, new Text(labels[i]), sliders[i], spinners[i]);
		}
		
		VBox all = new VBox(V_SPACE, presetChooser, grid);
		all.setAlignment(Pos.CENTER);
		all.managedProperty().bind(all.visibleProperty()); //make it hide when it is hiding
		return all;
	}
	
	
	/**
	 * Create a grid of sliders and spinners not unlike the aspectSelector.
	 * @param parameterSetter - The function to execute when the parameters change.
	 * @return The full formatted Region.
	 */
	@SuppressWarnings("unchecked")
	protected Region buildParameterSelector(Procedure parameterSetter) {
		currentParams = new double[4];
		paramLabels = new Text[4];
		paramSliders = new Slider[4]; // I don't think any projection has more than four parameters
		final Spinner<Double> spin0 = new Spinner<Double>(0.,0.,0.); //yes, this is awkward. Java gets weird about arrays with generic types
		paramSpinners = (Spinner<Double>[])Array.newInstance(spin0.getClass(), 4);
		paramSpinners[0] = spin0;
		
		for (int i = 0; i < 4; i ++) {
			paramLabels[i] = new Text();
			paramSliders[i] = new Slider();
			if (i != 0)
				paramSpinners[i] = new Spinner<Double>(0.,0.,0.);
			link(paramSliders[i], paramSpinners[i], i, currentParams, (d)->d, parameterSetter,
					isChanging, suppressListeners);
		}
		
		for (int i = 0; i < 4; i ++) {
			GridPane.setHgrow(paramSliders[i], Priority.ALWAYS);
			paramSpinners[i].setEditable(true);
		}
		
		paramGrid = new GridPane();
		paramGrid.setVgap(V_SPACE);
		paramGrid.setHgap(H_SPACE);
		paramGrid.getColumnConstraints().addAll(
				new ColumnConstraints(Control.USE_COMPUTED_SIZE),
				new ColumnConstraints(), new ColumnConstraints(SPINNER_WIDTH));
		return paramGrid;
	}
	
	
	/**
	 * Create a block of niche options - specifically antimeridian cropping and whether
	 * there should be a graticule
	 * @param extendMap - The mutable boolean value to which to bind the "Crop at antimeridian" CheckBox.
	 * @param graticule - The mutable double value to which to bind the "Graticule" Spinner.
	 * @return The full formatted Region.
	 */
	protected Region buildOptionPane(Flag extendMap, MutableDouble graticule) {
		final CheckBox extendBox = new CheckBox("Extend map area"); //the CheckBox for whether there should be shown imagery beyond += 180 degrees
		extendBox.setSelected(extendMap.isSet());
		extendBox.setTooltip(new Tooltip("Show points beyond the normal boundaries of the map (so some locations are shown twice)"));
		extendBox.selectedProperty().addListener((observable, old, now) -> extendMap.set(now));
		
		final ObservableList<Double> factorsOf90 = FXCollections.observableArrayList();
		for (double f = 1; f <= 90; f += 0.5)
			if (90%f == 0)
				factorsOf90.add(f);
		final Spinner<Double> gratSpinner = new Spinner<Double>(factorsOf90); //spinner for the graticule value
		gratSpinner.getValueFactory().setConverter(new DoubleStringConverter());
		gratSpinner.getValueFactory().setValue(15.);
		gratSpinner.setDisable(true);
		gratSpinner.setEditable(true);
		gratSpinner.setTooltip(new Tooltip("The spacing in degrees between shown parallels and meridians."));
		gratSpinner.setPrefWidth(SPINNER_WIDTH);
		gratSpinner.valueProperty().addListener((observable, old, now) -> {
				graticule.set(now); //which is tied to the mutable graticule spacing variable
			});
		
		final CheckBox gratBox = new CheckBox("Graticule:"); //the CheckBox for whether there should be a graticule
		gratBox.setSelected(false);
		gratBox.setTooltip(new Tooltip("Overlay a mesh of parallels and meridians."));
		gratBox.selectedProperty().addListener((observable, old, now) -> {
				if (now)
					graticule.set(gratSpinner.getValue()); //set the value of graticule appropriately when checked
				else
					graticule.set(0); //when not checked, represent "no graticule" as a spacing of 0
				gratSpinner.setDisable(!now); //disable the graticule Spinner when appropriate
			});
		
		final HBox gratRow = new HBox(H_SPACE, gratBox, gratSpinner);
		gratRow.setAlignment(Pos.CENTER_LEFT);
		return new VBox(V_SPACE, extendBox, gratRow);
	}
	
	
	/**
	 * Create a default button that will update the map.
	 * @return The button.
	 */
	protected Region buildUpdateButton(String text, Supplier<Task<SavableImage>> mapUpdater) {
		Button updateButton = new Button(text);
		updateButton.setOnAction((event) -> {
			Task<SavableImage> updateTask = mapUpdater.get();
			disableWhile(updateTask.runningProperty(),
					ButtonType.UPDATE_MAP, ButtonType.SAVE_GRAPH, ButtonType.SAVE_MAP);
			new Thread(updateTask).start();
		});
		updateButton.setTooltip(new Tooltip(
				"Update the current map with your parameters."));
		
		updateButton.setDefaultButton(true);
		root.addEventHandler(KeyEvent.KEY_PRESSED, (event) -> {
			if (CTRL_ENTER.match(event)) {
				updateButton.requestFocus();
				updateButton.fire();
			}
		});
		
		this.buttons.put(ButtonType.UPDATE_MAP, updateButton);
		return updateButton;
	}
	
	
	/**
	 * Build a button that will save something
	 * @param savee - The name of the thing being saved.
	 * @param allowedExtensions - The allowed file formats that can be saved.
	 * @param defaultExtension - The default file format to be saved.
	 * @param mapVerifier - A callback to run just before the saving happens that returns true if it should commence.
	 * @param mapCalculator - The callback that saves the thing.
	 * @return The button, ready to be pressed.
	 */
	protected Region buildSaveButton(String savee,
	                                 FileChooser.ExtensionFilter[] allowedExtensions,
	                                 FileChooser.ExtensionFilter defaultExtension,
	                                 BooleanSupplier mapVerifier, Supplier<Task<SavableImage>> mapCalculator) {
		FileChooser saver = new FileChooser();
		saver.setInitialDirectory(new File("output"));
		saver.setInitialFileName("my"+savee+defaultExtension.getExtensions().get(0).substring(1));
		saver.setTitle("Save "+savee);
		saver.getExtensionFilters().addAll(allowedExtensions);
		saver.setSelectedExtensionFilter(defaultExtension);
		try {
			if (!saver.getInitialDirectory().exists())
				saver.getInitialDirectory().mkdirs();
		} catch (SecurityException ignored) {}
		
		final Button saveButton = new Button("Save "+savee+"\u2026");
		final ButtonType buttonType =
				savee.equals("map") ? ButtonType.SAVE_MAP : ButtonType.SAVE_GRAPH;
		this.buttons.put(buttonType, saveButton);
		saveButton.setOnAction((actionEvent) -> { //when this save button is pressed...
			File file;
			try {
				file = saver.showSaveDialog(root); //let the user pick a file
			} catch(IllegalArgumentException e) {
				saver.setInitialDirectory(new File("."));
				file = saver.showSaveDialog(root);
			}
			if (file != null) { //progress only if a File was chosen
				saver.setInitialDirectory(file.getParentFile()); //remember this directory for next time
				final File f = file;
				
				if (mapVerifier.getAsBoolean()) { //if the optional verification process verifies it (TODO: do this before choosing a file)
					Task<SavableImage> task = mapCalculator.get(); //create the Task
					ProgressDialog<SavableImage> pBar = new ProgressDialog<SavableImage>(task); //and track its progress all the while
					pBar.show();
					task.setOnSucceeded((succeedEvent) -> { //save the result to disk when it finishes
						pBar.unbindAndSetProgress(-1);
						pBar.unbindAndSetHeaderText("Saving to disk\u2026");
						try {
							task.getValue().save(f);
						} catch (IOException e) {
							showError("Failure!",
									"Could not access "+f.getAbsolutePath()+". It's possible that another program has it open.");
						}
						pBar.close();
//						// TODO: show a popup here with a button to open the folder
					});
					disableWhile(task.runningProperty(), buttonType);
					new Thread(task).start(); //now execute!
				}
			}
		});
		saveButton.setTooltip(new Tooltip("Save the "+savee+" with current settings"));
		
		// bind ctrl+S to save
		root.addEventHandler(KeyEvent.KEY_PRESSED, (event) -> {
			if (CTRL_S.match(event)) {
				saveButton.requestFocus();
				saveButton.fire();
			}
		});
		
		return saveButton;
	}
	
	
	protected Projection getProjection() {
		return projectionChooser.getValue();
	}
	
	
	protected boolean getParamsChanging() { //are the aspect or parameters actively changing?
		return isChanging.isSet();
	}
	
	
	protected void loadParameters() {
		try {
			getProjection().initialize(currentParams);
		} catch (IllegalArgumentException e) {
			showError("Failed to load projection", e.getLocalizedMessage());
		}
	}
	
	
	protected void disableWhile(ReadOnlyBooleanProperty condition, ButtonType... buttons) {
		for (ButtonType bt: buttons)
			if (this.buttons.containsKey(bt))
				this.buttons.get(bt).disableProperty().bind(condition);
	}
	
	
	private void chooseProjectionFromExpandedList(Projection lastProjection) {
		final ProjectionSelectionDialog selectDialog = new ProjectionSelectionDialog();
		
		do {
			final Optional<Projection> result = selectDialog.showAndWait();
			
			if (result.isPresent()) {
				if (result.get() == Projection.NULL_PROJECTION) {
					showError("No projection chosen", "Please select a projection.");
				}
				else {
					Platform.runLater(() -> projectionChooser.setValue(result.get()));
					break;
				}
			}
			else {
				Platform.runLater(() -> projectionChooser.setValue(lastProjection)); //I need these runLater()s because JavaFX is dumb and throws an obscure error on recursive edits to ComboBoxes
				break;
			}
		} while (true);
	}
	
	
	private void setAspectByPreset(String presetName,
			Slider[] sliders) {
		this.suppressListeners.set();
		if (presetName.equals("Antipode")) {
			sliders[0].setValue(-sliders[0].getValue());
			sliders[1].setValue((sliders[1].getValue()+360)%360-180);
			sliders[2].setValue(-sliders[2].getValue());
		}
		else if (presetName.equals("Random")) {
			sliders[0].setValue(toDegrees(asin(random()*2-1)));
			sliders[1].setValue(random()*360-180);
			sliders[2].setValue(random()*360-180);
		}
		else {
			for (int i = 0; i < ASPECT_NAMES.length; i ++) {
				if (ASPECT_NAMES[i].equals(presetName)) {
					for (int j = 0; j < 3; j ++)
						sliders[j].setValue(ASPECT_VALUES[j][i]);
					break;
				}
			}
		}
		
		this.suppressListeners.clear();
	}
	
	
	private void revealParameters(Projection proj) {
		this.suppressListeners.set();
		final String[] paramNames = proj.getParameterNames();
		final double[][] paramValues = proj.getParameterValues();
		paramGrid.getChildren().clear();
		for (int i = 0; i < proj.getNumParameters(); i ++) {
			paramLabels[i].setText(paramNames[i]+":");
			paramSliders[i].setMin(paramValues[i][0]);
			paramSliders[i].setMax(paramValues[i][1]);
			paramSliders[i].setValue(paramValues[i][2]);
			final SpinnerValueFactory.DoubleSpinnerValueFactory paramSpinVF =
					(DoubleSpinnerValueFactory) paramSpinners[i].getValueFactory();
			paramSpinVF.setMin(paramValues[i][0]);
			paramSpinVF.setMax(paramValues[i][1]);
			paramSpinVF.setValue(paramValues[i][2]);
			final Tooltip tt = new Tooltip(
					"Change the "+paramNames[i]+" of the map (default is " + paramValues[i][2]+")");
			paramSliders[i].setTooltip(tt);
			paramSpinners[i].setTooltip(tt);
			paramGrid.addRow(i, paramLabels[i], paramSliders[i], paramSpinners[i]);
		}
		this.suppressListeners.clear();
	}
	
	
	private static void link(Slider sld, Spinner<Double> spn, int i, double[] doubles,
		DoubleUnaryOperator converter, Procedure callback, Flag isChanging, Flag suppressListeners) {
		sld.valueChangingProperty().addListener((observable, prev, now) -> { //link spinner to slider
				isChanging.set(now);
				if (!now) {
					if (spn.getValue() != sld.getValue())
						spn.getValueFactory().setValue(sld.getValue());
					doubles[i] = converter.applyAsDouble(round(sld.getValue(),3));
					if (!suppressListeners.isSet())
						callback.execute();
				}
			});
		sld.valueProperty().addListener((observable, prev, now) -> {
				if (spn.getValue() != sld.getValue())
					spn.getValueFactory().setValue(sld.getValue());
				doubles[i] = converter.applyAsDouble(round(now.doubleValue(),3));
				if (!suppressListeners.isSet())
					callback.execute();
			});
		
		spn.valueProperty().addListener((observable, prev, now) -> { //link slider to spinner
				if (spn.getValue() != sld.getValue())
					sld.setValue(spn.getValue());
			});
		
		spn.focusedProperty().addListener((observable, prev, now) -> { //make spinner act rationally
				if (!now) 	spn.increment(0);
			});
	}
	
	
	protected static void showError(String header, String message) { //a simple error handling thing
		final Alert alert = new Alert(Alert.AlertType.ERROR);
		alert.setHeaderText(header);
		Text contentText = new Text(message);
		contentText.setWrappingWidth(350);
		StackPane contentPane = new StackPane(contentText);
		contentPane.setPadding(new Insets(10, 10, 10, 10));
		contentPane.setAlignment(Pos.CENTER);
		alert.getDialogPane().setContent(contentPane);
		alert.showAndWait();
	}
	
	
	
	protected enum ButtonType {
		LOAD_INPUT, UPDATE_MAP, SAVE_MAP, SAVE_GRAPH
	}
}
