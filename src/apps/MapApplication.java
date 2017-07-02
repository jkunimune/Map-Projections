/**
 * 
 */
package apps;

import java.io.File;
import java.lang.reflect.Array;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
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
import javafx.scene.control.Label;
import javafx.scene.control.MenuButton;
import javafx.scene.control.MenuItem;
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

	protected static final int CONT_WIDTH = 300;
	protected static final int IMG_WIDTH = 500;
	
	private static final KeyCombination ctrlO = new KeyCodeCombination(KeyCode.O, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination ctrlS = new KeyCodeCombination(KeyCode.S, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination ctrlEnter = new KeyCodeCombination(KeyCode.ENTER, KeyCodeCombination.CONTROL_DOWN);
	
	private static final String[] ASPECT_NAMES = { "Standard", "Transverse", "Center of Mass", "Jerusalem", "Point Nemo",
			"Longest Line", "Longest Line Transverse", "Cylindrical", "Conic", "Tetrahedral", "Quincuncial", "Antipode", "Random" };
	private static final double[][] ASPECT_VALS = {
			{ 90., 0., 29.9792, 31.7833, 48.8767, -28.5217,-46.4883,-35.,    -10.,  47., 60. },
			{  0., 0., 31.1344, 35.2160, 56.6067, 141.451,  16.5305,-13.6064, 65.,-173., -6. },
			{  0., 0.,-32.,    -35.,    -45.,     161.5,   137.,    145.,   -150., 138.,-10. } };
	
	
	protected Stage root;
	final protected String name;
	private ComboBox<Projection> projectionChooser;
	
	
	
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
			Consumer<File> inputSetter) {
		final Label label = new Label("Current input:");
		final Text inputLabel = new Text("None");
		
		final FileChooser inputChooser = new FileChooser();
		inputChooser.setInitialDirectory(new File("input"));
		inputChooser.setTitle("Choose an input map");
		inputChooser.getExtensionFilters().addAll(allowedExtensions);
		
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
		root.addEventHandler(KeyEvent.KEY_PRESSED, new EventHandler<KeyEvent>() {	// ctrl-O opens
			public void handle(KeyEvent event) {
				if (ctrlO.match(event))	loadButton.fire();
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
		description.setWrappingWidth(CONT_WIDTH);
		
		projectionChooser.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				description.setText(projectionChooser.getValue().getDescription());
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
		
		final Slider[] sliders = new Slider[] {
				new Slider(-90, 90, 90),
				new Slider(-180,180,0),
				new Slider(-180,180,0)
		};
		copyToArray(sliders, aspectArr);
		
		final Spinner<Double> latSpinner = new Spinner<Double>(-90, 90, 90.0); //yes, this is awkward. Java gets weird about arrays with generic types
		@SuppressWarnings("unchecked")
		final Spinner<Double>[] spinners = (Spinner<Double>[]) Array.newInstance(latSpinner.getClass(), 3);
		spinners[0] = latSpinner;
		spinners[1] = new Spinner<Double>(-180, 180, 0.0);
		spinners[2] = new Spinner<Double>(-180, 180, 0.0);
		
		for (String preset: ASPECT_NAMES) {
			MenuItem m = new MenuItem(preset);
			m.setOnAction(new EventHandler<ActionEvent>() {
				public void handle(ActionEvent event) {
					setAspectByPreset(((MenuItem) event.getSource()).getText(),
							sliders, spinners);
					copyToArray(sliders, aspectArr);
					aspectSetter.execute();
				}
			});
			presetChooser.getItems().add(m);
		}
		
		for (int i = 0; i < 3; i ++) {
			final Slider sld = sliders[i];
			final Spinner<Double> spn = spinners[i];
			GridPane.setHgrow(sld, Priority.ALWAYS);
			sld.setTooltip(new Tooltip("Change the aspect of the map"));
			sld.valueChangingProperty().addListener(
					(observable, then, now) -> {
						if (spn.getValue() != sld.getValue()) //TODO the boxes need to change when the dropdown fires
							spn.getEditor().textProperty().set(Double.toString(sld.getValue()));
						copyToArray(sliders, aspectArr);
						aspectSetter.execute();
					});
			
			spn.setTooltip(new Tooltip("Change the aspect of the map"));
			spn.setPrefWidth(100);
			spn.setEditable(true);
			spn.getEditor().textProperty().addListener((ov, pv, nv) -> {	// link the Spinners
				if (spn.getEditor().textProperty().isEmpty().get()) 	return;
				if (nv.equals("-") || nv.equals(".")) 					return;
				try {
					Double.parseDouble(nv);
					spn.increment(0);	// forces the spinner to commit its value
					if (spn.getValue() != sld.getValue())
						sld.setValue(spn.getValue());
					copyToArray(sliders, aspectArr);
					aspectSetter.execute();
				} catch (NumberFormatException e) {
					spn.getEditor().textProperty().set(pv); //yeah, this is all pretty jank. JavaFX spinners are just weird by default
				}
			});
		}
		
		final GridPane grid = new GridPane();
		grid.addRow(0, new Text("Latitude:"), sliders[0], spinners[0]);
		grid.addRow(1, new Text("Longitude:"), sliders[1], spinners[1]);
		grid.addRow(2, new Text("Ctr. Meridian:"), sliders[2], spinners[2]);
		
		VBox all = new VBox(5, presetChooser, grid);
		all.setAlignment(Pos.CENTER);
		return all;
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
				"Update the current map with your parameters."));
		
		updateButton.setDefaultButton(true);
		root.addEventHandler(KeyEvent.KEY_PRESSED, (KeyEvent event) -> {
			if (ctrlEnter.match(event))	updateButton.fire(); });
		
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
			Procedure parameterCollector,
			BiConsumer<File, ProgressBarDialog> calculateAndSaver) {
		final FileChooser saver = new FileChooser();
		saver.setInitialDirectory(new File("output"));
		saver.setInitialFileName("my"+savee+allowedExtensions[0].getExtensions().get(0));
		saver.setTitle("Save "+savee);
		saver.getExtensionFilters().addAll(allowedExtensions);
		
		final Button saveButton = new Button("Save "+savee+"...");
		saveButton.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				final File f = saver.showOpenDialog(root);
				if (f != null) {
					parameterCollector.execute();
					final ProgressBarDialog pBar = new ProgressBarDialog();
					pBar.show();
					trigger(saveButton, () -> {
							calculateAndSaver.accept(f, pBar);
							Platform.runLater(pBar::close); });
				}
			}
		});
		saveButton.setTooltip(new Tooltip("Save the "+savee+" with current settings."));
		
		if (bindCtrlS) // ctrl+S saves
			root.addEventHandler(KeyEvent.KEY_PRESSED, (KeyEvent event) -> {
				if (ctrlS.match(event))	saveButton.fire(); });
		
		return saveButton;
	}
	
	
	protected void setAspectByPreset(String presetName,
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
	
	
	protected Projection getProjection() {
		return projectionChooser.getValue();
	}
	
	
	private static void copyToArray(Slider[] sliders, double[] output) {
		assert sliders.length == 3; //ugh I can never remember how to accually activate these things. I'll just leave it here as a self-commenting thing
		for (int i = 0; i < 3; i ++)
			output[i] = Math.toRadians(sliders[i].getValue());
	}

}
