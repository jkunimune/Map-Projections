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
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.Scene;
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
import maps.Projection;
import util.Procedure;


/**
 * A base class for all applications that deal with maps
 * 
 * @author jkunimune
 */
public abstract class MapApplication extends Application {

	protected static final int GUI_WIDTH = 350;
	protected static final int IMG_WIDTH = 500;
	
	private static final KeyCombination ctrlO = new KeyCodeCombination(KeyCode.O, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination ctrlS = new KeyCodeCombination(KeyCode.S, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination ctrlEnter = new KeyCodeCombination(KeyCode.ENTER, KeyCodeCombination.CONTROL_DOWN);
	
	private static final String[] ASPECT_NAMES = { "Standard", "Transverse", "Center of Mass", "Jerusalem", "Point Nemo",
			"Longest Line", "Longest Line Transverse", "Cylindrical", "Conic", "Tetrahedral", "Quincuncial", "Antipode", "Random" };
	private static final double[][] ASPECT_VALS = {
			{ 90., 0., 29.98, 31.78, 48.88, -28.52,-46.4883,-35.,  -10.,  47., 60. },
			{  0., 0., 31.13, 35.22, 56.61, 141.45, 16.5305,-13.61, 65.,-173., -6. },
			{  0., 0.,-32.,  -35.,  -45.,   161.5, 137.,    145., -150., 138.,-10. } };
	
	
	final private String name;
	private Stage root;
	private ComboBox<Projection> projectionChooser;
	private GridPane paramGrid;
	private Label[] paramLabels;
	private Slider[] paramSliders;
	private Spinner<Double>[] paramSpinners;
	private double[] currentParams;
	
	
	
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
		
		final FileChooser inputChooser = new FileChooser();
		inputChooser.setInitialDirectory(new File("input"));
		inputChooser.setTitle("Choose an input map");
		inputChooser.getExtensionFilters().addAll(allowedExtensions);
		inputChooser.setSelectedExtensionFilter(defaultExtension);
		
		final Button loadButton = new Button("Choose input...");
		loadButton.setTooltip(new Tooltip(
				"Choose the image to determine the style of your map"));
		loadButton.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				final File f = inputChooser.showOpenDialog(root);
				if (f != null) {
					trigger(loadButton, () -> inputSetter.accept(f));
					inputLabel.setText(f.getName());
				}
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
	protected Region buildProjectionSelector(Projection[] projections, Projection defProj, Procedure projectionSetter) {
		final Label label = new Label("Projection:");
		projectionChooser =
				new ComboBox<Projection>(FXCollections.observableArrayList(projections));
		projectionChooser.setPrefWidth(210);
		projectionChooser.setValue(defProj);
		
		final Text description = new Text(projectionChooser.getValue().getDescription());
		description.setWrappingWidth(GUI_WIDTH);
		
		projectionChooser.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				description.setText(projectionChooser.getValue().getDescription());
				revealParameters(projectionChooser.getValue());
				projectionSetter.execute();
			}
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
				new Slider(-90, 90, 90), //TODO: can we call setAspectByPreset("Standard") instead of this?
				new Slider(-180,180,0),
				new Slider(-180,180,0) };
		final Spinner<Double> spin0 = new Spinner<Double>(-90, 90, 90.0); //yes, this is awkward. Java gets weird about arrays with generic types
		@SuppressWarnings("unchecked")
		final Spinner<Double>[] spinners = (Spinner<Double>[]) Array.newInstance(spin0.getClass(), 3);
		spinners[0] = spin0;
		spinners[1] = new Spinner<Double>(-180, 180, 0.0);
		spinners[2] = new Spinner<Double>(-180, 180, 0.0);
		
		for (int i = 0; i < 3; i ++)
			aspectArr[i] = Math.toRadians(sliders[i].getValue());
		
		for (String preset: ASPECT_NAMES) {
			MenuItem m = new MenuItem(preset);
			m.setOnAction(new EventHandler<ActionEvent>() {
				public void handle(ActionEvent event) {
					setAspectByPreset(((MenuItem) event.getSource()).getText(),
							sliders, spinners);
					for (int i = 0; i < 3; i ++)
						aspectArr[i] = Math.toRadians(sliders[i].getValue());
					aspectSetter.execute();
				}
			});
			presetChooser.getItems().add(m);
		}
		
		link(sliders, spinners, aspectArr, Math::toRadians, aspectSetter);
		
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
		currentParams = new double[3];
		paramLabels = new Label[] {
				new Label(), new Label(), new Label() };
		paramSliders = new Slider[] {
				new Slider(), new Slider(), new Slider() }; // I don't think any projection has more than three parameters
		
		final Spinner<Double> spin0 = new Spinner<Double>(0.,0.,0.); //yes, this is awkward. Java gets weird about arrays with generic types
		paramSpinners = (Spinner<Double>[]) Array.newInstance(spin0.getClass(), 3);
		paramSpinners[0] = spin0;
		paramSpinners[1] = new Spinner<Double>(0.,0.,0.);
		paramSpinners[2] = new Spinner<Double>(0.,0.,0.);
		
		link(paramSliders, paramSpinners, currentParams, (d)->d, parameterSetter);
		
		for (int i = 0; i < 3; i ++) {
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
		updateButton.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				trigger(updateButton, mapUpdater);
			}
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
	 * @param bindCtrlS true iff this button should have ^S bound to it
	 * @param savee The name of the thing being saved
	 * @param allowedExtensions The List of possible file types for the output
	 * @param startSaving The method that will be called when the button is
	 * 		pressed. This will generally be called from the main thread. maybe.
	 * @return the button
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
		
		final Button saveButton = new Button("Save "+savee+"...");
		saveButton.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				final File f = saver.showSaveDialog(root);
				if (f != null) {
					if (settingCollector.getAsBoolean()) {
						final ProgressBarDialog pBar = new ProgressBarDialog();
						pBar.setContentText("Finalizing "+savee+"...");
						pBar.show();
						trigger(saveButton, () -> {
								calculateAndSaver.accept(f, pBar);
								Platform.runLater(pBar::close); });
					}
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
	
	
	protected double[] getParams() {
		return currentParams;
	}
	
	
	protected static final void trigger(Button btn, Runnable task) {
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
	
	
	private void setAspectByPreset(String presetName,
			Slider[] sliders, Spinner<Double>[] spinners) {
		if (presetName.equals("Antipode")) {
			sliders[0].setValue(-sliders[0].getValue());
			sliders[1].setValue((sliders[1].getValue()+180)%360);
			sliders[2].setValue(-sliders[2].getValue());
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
	}
	
	
	private void revealParameters(Projection proj) {
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
	}
	
	
	private static void link(Slider[] sliders, Spinner<Double>[] spinners, double[] doubles,
			DoubleUnaryOperator converter, Procedure callback) {
		for (int i = 0; i < doubles.length; i ++) {
			final Slider sld = sliders[i];
			final Spinner<Double> spn = spinners[i];
			final int I = i;
			
			sld.valueChangingProperty().addListener((observable, prev, now) -> { //link spinner to slider
					if (prev) {
						if (spn.getValue() != sld.getValue())
							spn.getValueFactory().setValue(sld.getValue());
						doubles[I] = converter.applyAsDouble(sld.getValue());
						callback.execute();
					}
				});
			
			spn.valueProperty().addListener((observable, prev, now) -> { //link slider to spinner
					if (spn.getValue() != sld.getValue())
						sld.setValue(spn.getValue());
					doubles[I] = converter.applyAsDouble(spn.getValue());
					callback.execute();
				});
			
			spn.focusedProperty().addListener((observable, prev, now) -> { //make spinner act rationally
					if (!now) 	spn.increment(0);
				});
		}
	}

}
