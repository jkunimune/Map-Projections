/**
 * 
 */
package apps;

import java.io.File;
import java.lang.reflect.Array;
import java.util.function.BiConsumer;
import java.util.function.BooleanSupplier;
import java.util.function.Consumer;
import java.util.function.DoubleUnaryOperator;

import dialogs.ProgressBarDialog;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
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
import maps.Cylindrical;
import maps.Projection;
import utils.Math2;
import utils.Procedure;


/**
 * A base class for all GUI applications that deal with maps
 * 
 * @author jkunimune
 */
public abstract class MapApplication extends Application {

	protected static final int GUI_WIDTH = 350;
	protected static final int IMG_WIDTH = 500;
	
	private static final KeyCombination ctrlO = new KeyCodeCombination(KeyCode.O, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination ctrlS = new KeyCodeCombination(KeyCode.S, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination ctrlEnter = new KeyCodeCombination(KeyCode.ENTER, KeyCodeCombination.CONTROL_DOWN);
	
	private static final String[] ASPECT_NAMES = { "Standard", "Transverse", "Cassini",
			"Atlantis", "AuthaGraph", "Jerusalem", "Point Nemo", "Longest Line",
			"Cylindrical", "Tetrahedral", "Antipode", "Random" };
	private static final double[][] ASPECT_VALS = { //the aspect presets (in degrees)
			{ 90., 0.,  0.,  -4.,  77., 31.78, 48.88, -28.52,-35.,    57. },
			{  0., 0., 90.,  65., 143., 35.22, 56.61, 141.45,-13.61,-176. },
			{  0., 0.,-90.,-147.,  17.,-35.,  -45.,    30.,  145.,   154. } };
	
	
	final private String name;
	private Stage root;
	private ComboBox<Projection> projectionChooser;
	private GridPane paramGrid;
	private Label[] paramLabels;
	private Slider[] paramSliders;
	private Spinner<Double>[] paramSpinners;
	private double[] currentParams;
	
	private Flag isChanging = new Flag();
	private Flag suppressListeners = new Flag();
	
	
	
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
		this.projectionChooser.setValue(Cylindrical.MERCATOR);
		this.suppressListeners.clear();
	}
	
	
	protected abstract Node makeWidgets();
	
	
	/**
	 * Create a set of widgets to select an input image
	 * @param allowedExtensions The List of possible file types for the input
	 * @param setInput The method to be called when an input file is loaded.
	 * 		This method will be called from a nongui thread.
	 * @return the full formatted Region
	 */
	protected Region buildInputSelector(
			FileChooser.ExtensionFilter[] allowedExtensions,
			FileChooser.ExtensionFilter defaultExtension,
			Consumer<File> inputSetter) {
		final Label label = new Label("Current input:");
		final Text inputLabel = new Text("None");
		
		final FileChooser inputChooser = new FileChooser(); //TODO: remember last directory
		inputChooser.setInitialDirectory(new File("input"));
		inputChooser.setTitle("Choose an input map");
		inputChooser.getExtensionFilters().addAll(allowedExtensions);
		inputChooser.setSelectedExtensionFilter(defaultExtension);
		
		final Button loadButton = new Button("Choose input...");
		loadButton.setTooltip(new Tooltip(
				"Choose the image to determine the style of your map"));
		loadButton.setOnAction((event) -> {
				File file;
				try {
					file = inputChooser.showOpenDialog(root);
				} catch (IllegalArgumentException e) {
					inputChooser.setInitialDirectory(new File("."));
					file = inputChooser.showOpenDialog(root);
				}
				if (file != null) {
					final File f = file;
					trigger(loadButton, () -> inputSetter.accept(f));
					inputLabel.setText(f.getName());
				}
			});
		loadButton.setTooltip(new Tooltip(
				"Change the input image"));
		root.addEventHandler(KeyEvent.KEY_PRESSED, (event) -> {	// ctrl-O opens
				if (ctrlO.match(event)) {
					loadButton.requestFocus();
					loadButton.fire();
				}
			});
		
		VBox output = new VBox(5, new HBox(3, label, inputLabel), loadButton);
		output.setAlignment(Pos.CENTER);
		return output;
	}
	
	
	/**
	 * Create a set of widgets to choose a Projection
	 * @param projections The Projections from which the user may choose
	 * @param defProj The default projection, before the user chooses anything
	 * @return the full formatted Region
	 */
	protected Region buildProjectionSelector(Projection[] projections, Procedure projectionSetter) {
		final Label label = new Label("Projection:");
		projectionChooser =
				new ComboBox<Projection>(FXCollections.observableArrayList(projections));
		projectionChooser.setPrefWidth(210);
		
		final Text description = new Text();
		description.setWrappingWidth(GUI_WIDTH);
		
		projectionChooser.setOnAction((event) -> {
				final boolean suppressedListeners = suppressListeners.isSet(); //save this value, because
				description.setText(projectionChooser.getValue().getDescription()); //revealParameters
				revealParameters(projectionChooser.getValue()); //clears suppressListeners. That's fine,
				if (!suppressedListeners) //because suppressListeners is only needed here for that one case.
					projectionSetter.execute();
			});
		
		return new VBox(5, new HBox(3, label, projectionChooser), description);
	}
	
	
	/**
	 * Create a set of widgets to choose an aspect either from a preset or numbers
	 * Also bind aspectArr to the sliders
	 * @return the full formatted Region
	 */
	protected Region buildAspectSelector(double[] aspectArr, Procedure aspectSetter) {
		final MenuButton presetChooser = new MenuButton("Aspect Presets");
		presetChooser.setTooltip(new Tooltip(
				"Set the aspect sliders based on a preset"));
		
		final String[] labels = { "Latitude:", "Longitude:", "Ctr. Meridian:" };
		final Slider[] sliders = new Slider[] {
				new Slider(-90, 90, 0.), //TODO: can we call setAspectByPreset("Standard") instead of this?
				new Slider(-180, 180, 0.),
				new Slider(-180, 180, 0.) };
		final Spinner<Double> spin0 = new Spinner<Double>(-90, 90, 0.); //yes, this is awkward. Java gets weird about arrays with generic types
		@SuppressWarnings("unchecked")
		final Spinner<Double>[] spinners = (Spinner<Double>[]) Array.newInstance(spin0.getClass(), 3);
		spinners[0] = spin0;
		spinners[1] = new Spinner<Double>(-180, 180, 0.);
		spinners[2] = new Spinner<Double>(-180, 180, 0.);
		
		for (int i = 0; i < 3; i ++) {
			aspectArr[i] = Math.toRadians(sliders[i].getValue());
			link(sliders[i], spinners[i], i, aspectArr, Math::toRadians,
					aspectSetter, isChanging, suppressListeners);
		}
		setAspectByPreset("Standard", sliders, spinners);
		
		for (String preset: ASPECT_NAMES) {
			MenuItem m = new MenuItem(preset);
			m.setOnAction((event) -> {
					setAspectByPreset(((MenuItem) event.getSource()).getText(),
							sliders, spinners);
					for (int i = 0; i < 3; i ++)
						aspectArr[i] = Math.toRadians(sliders[i].getValue());
					if (!suppressListeners.isSet())
						aspectSetter.execute();
				});
			presetChooser.getItems().add(m);
		}
		
		final GridPane grid = new GridPane();
		grid.setVgap(5);
		grid.setHgap(3);
		grid.getColumnConstraints().addAll(
				new ColumnConstraints(92,Control.USE_COMPUTED_SIZE,Control.USE_COMPUTED_SIZE),
				new ColumnConstraints(), new ColumnConstraints(92));
		for (int i = 0; i < 3; i ++) {
			GridPane.setHgrow(sliders[i], Priority.ALWAYS);
			sliders[i].setTooltip(new Tooltip("Change the aspect of the map"));
			spinners[i].setTooltip(new Tooltip("Change the aspect of the map"));
			spinners[i].setEditable(true);
			grid.addRow(i, new Label(labels[i]), sliders[i], spinners[i]);
		}
		
		VBox all = new VBox(5, presetChooser, grid);
		all.setAlignment(Pos.CENTER);
		return all;
	}
	
	
	/**
	 * Create a grid of sliders and spinners not unlike the aspectSelector
	 * @param parameterSetter The function to execute when the parameters change
	 * @return the full formatted Region
	 */
	@SuppressWarnings("unchecked")
	protected Node buildParameterSelector(Procedure parameterSetter) {
		currentParams = new double[4];
		paramLabels = new Label[4];
		paramSliders = new Slider[4]; // I don't think any projection has more than four parameters
		final Spinner<Double> spin0 = new Spinner<Double>(0.,0.,0.); //yes, this is awkward. Java gets weird about arrays with generic types
		paramSpinners = (Spinner<Double>[]) Array.newInstance(spin0.getClass(), 4);
		paramSpinners[0] = spin0;
		
		for (int i = 0; i < 4; i ++) {
			paramLabels[i] = new Label();
			paramSliders[i] = new Slider();
			if (i != 0)
				paramSpinners[i] = new Spinner<Double>(0.,0.,0.);
			link(paramSliders[i], paramSpinners[i], i, currentParams, (d)->d, parameterSetter, isChanging, suppressListeners);
		}
		
		for (int i = 0; i < 4; i ++) {
			GridPane.setHgrow(paramSliders[i], Priority.ALWAYS);
			paramSpinners[i].setEditable(true);
		}
		
		paramGrid = new GridPane();
		paramGrid.setVgap(5);
		paramGrid.setHgap(3);
		paramGrid.getColumnConstraints().addAll(
				new ColumnConstraints(92,Control.USE_COMPUTED_SIZE,Control.USE_COMPUTED_SIZE),
				new ColumnConstraints(), new ColumnConstraints(92));
		return paramGrid;
	}
	
	
	/**
	 * Create a default button that will update the map
	 * @return the button
	 */
	protected Button buildUpdateButton(Runnable mapUpdater) {
		final Button updateButton = new Button("Update map");
		updateButton.setOnAction((event) -> {
				trigger(updateButton, mapUpdater);
			});
		updateButton.setTooltip(new Tooltip(
				"Update the current map with your parameters"));
		
		updateButton.setDefaultButton(true);
		root.addEventHandler(KeyEvent.KEY_PRESSED, (event) -> {
				if (ctrlEnter.match(event)) {
					updateButton.requestFocus();
					updateButton.fire();
				}
			});
		
		return updateButton;
	}
	
	
	/**
	 * Build a button that will save something
	 * @param bindCtrlS Should ctrl+S trigger this button?
	 * @param savee The name of the thing being saved
	 * @param allowedExtensions The allowed file formats that can be saved
	 * @param defaultExtension The default file format to be saved
	 * @param settingCollector A callback to run just before the saving happens that returns true if it should commence
	 * @param calculateAndSaver The callback that saves the thing
	 * @return the button, ready to be pressed
	 */
	protected Button buildSaveButton(boolean bindCtrlS, String savee,
			FileChooser.ExtensionFilter[] allowedExtensions,
			FileChooser.ExtensionFilter defaultExtension,
			BooleanSupplier settingCollector,
			BiConsumer<File, ProgressBarDialog> calculateAndSaver) {
		final FileChooser saver = new FileChooser();
		saver.setInitialDirectory(new File("output"));
		saver.setInitialFileName("my"+savee+defaultExtension.getExtensions().get(0).substring(1));
		saver.setTitle("Save "+savee);
		saver.getExtensionFilters().addAll(allowedExtensions);
		saver.setSelectedExtensionFilter(defaultExtension);
		try {
			if (!saver.getInitialDirectory().exists())
				saver.getInitialDirectory().mkdirs();
		} catch (SecurityException e) {}
		
		final Button saveButton = new Button("Save "+savee+"...");
		saveButton.setOnAction((event) -> {
				File file;
				try {
					file = saver.showSaveDialog(root);
				} catch(IllegalArgumentException e) {
					saver.setInitialDirectory(new File("."));
					file = saver.showSaveDialog(root);
				}
				if (file != null) {
					final File f = file;
					if (settingCollector.getAsBoolean()) {
						final ProgressBarDialog pBar = new ProgressBarDialog();
						pBar.setContentText("Finalizing "+savee+"...");
						pBar.show();
						trigger(saveButton, () -> {
								calculateAndSaver.accept(f, pBar);
								Platform.runLater(pBar::close);
							});
					}
				}
			});
		saveButton.setTooltip(new Tooltip("Save the "+savee+" with current settings"));
		
		if (bindCtrlS) // ctrl+S saves
			root.addEventHandler(KeyEvent.KEY_PRESSED, (event) -> {
				if (ctrlS.match(event)) {
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
		getProjection().setParameters(currentParams);
	}
	
	
	private void setAspectByPreset(String presetName,
			Slider[] sliders, Spinner<Double>[] spinners) {
		this.suppressListeners.set();
		if (presetName.equals("Antipode")) {
			sliders[0].setValue(-sliders[0].getValue());
			sliders[1].setValue((sliders[1].getValue()+360)%360-180);
		}
		else if (presetName.equals("Random")) {
			sliders[0].setValue(Math.toDegrees(Math.asin(Math.random()*2-1)));
			sliders[1].setValue(Math.random()*360-180);
			sliders[2].setValue(Math.random()*360-180);
		}
		else {
			for (int i = 0; i < ASPECT_NAMES.length; i ++) {
				if (ASPECT_NAMES[i].equals(presetName)) {
					for (int j = 0; j < 3; j ++)
						sliders[j].setValue(ASPECT_VALS[j][i]);
					break;
				}
			}
		}
		
		for (int i = 0; i < 3; i ++)
			spinners[i].getEditor().textProperty().set(
					Double.toString(sliders[i].getValue()));
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
					"Change the "+paramNames[i]+" of the map (default is "+
					paramValues[i][2]+")");
			paramSliders[i].setTooltip(tt);
			paramSpinners[i].setTooltip(tt);
			paramGrid.addRow(i, paramLabels[i], paramSliders[i], paramSpinners[i]);
		}
		this.suppressListeners.clear();
	}
	
	
	private static void trigger(Button btn, Runnable task) {
		btn.setDisable(true);
		new Thread(() -> {
			try {
				task.run();
			} catch (Exception e) {
				e.printStackTrace();
			} finally {
				btn.setDisable(false);
			}
		}).start();
	}
	
	
	private static void link(Slider sld, Spinner<Double> spn, int i, double[] doubles,
		DoubleUnaryOperator converter, Procedure callback, Flag isChanging, Flag suppressListeners) {
		sld.valueChangingProperty().addListener((observable, prev, now) -> { //link spinner to slider
				isChanging.set(now);
				if (!now) {
					if (spn.getValue() != sld.getValue())
						spn.getValueFactory().setValue(sld.getValue());
					doubles[i] = converter.applyAsDouble(Math2.round(sld.getValue(),3));
					if (!suppressListeners.isSet())
						callback.execute();
				}
			});
		sld.valueProperty().addListener((observable, prev, now) -> {
				if (spn.getValue() != sld.getValue())
					spn.getValueFactory().setValue(sld.getValue());
				doubles[i] = converter.applyAsDouble(Math2.round(now.doubleValue(),3));
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
	
	
	protected static void showError(String header, String message) { //a simple thread-safe error handling thing
		Platform.runLater(() -> {
				final Alert alert = new Alert(Alert.AlertType.ERROR);
				alert.setHeaderText(header);
				alert.setContentText(message);
				alert.showAndWait();
			});
	}
	
	
	
	/**
	 * Because Java apparently doesn't already have a mutable Bullion
	 * @author jkunimune
	 */
	private class Flag {
		private boolean set;
		
		public Flag() {
			this(false);
		}
		
		public Flag(boolean set) {
			this.set = set;
		}
		
		public boolean isSet() {
			return this.set;
		}
		
		public void set(boolean set) {
			this.set = set;
		}
		
		public void set() {
			this.set = true;
		}
		
		public void clear() {
			this.set = false;
		}
		
		public String toString() {
			return "Flag("+this.set+")";
		}
	}

}
