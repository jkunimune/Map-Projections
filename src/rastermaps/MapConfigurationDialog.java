package rastermaps;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Pos;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Dialog;
import javafx.scene.control.DialogPane;
import javafx.scene.control.Label;
import javafx.scene.control.Spinner;
import javafx.scene.control.ButtonBar.ButtonData;
import javafx.scene.image.Image;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.VBox;

public class MapConfigurationDialog extends Dialog<Thread> {

	public final double DEF_SIZE = 1000;
	
	
	private final double defaultRatio;
	private final VBox gui;
	private final CheckBox maintainRatio;
	private final Spinner<Integer> widthBox, heightBox;
	private final ComboBox<String> smoothBox;
	private final MapProjections parent;
	
	private boolean realEdit; // is a real edit happening, or is it just me?
	
	public MapConfigurationDialog(double defAsp, MapProjections mp) {
		this.defaultRatio = defAsp;
		this.parent = mp;
		DialogPane pane = this.getDialogPane();
		
		this.maintainRatio = new CheckBox("Maintain aspect ratio");	// instantiate the components
		this.maintainRatio.setSelected(true);
		
		this.widthBox = new Spinner<Integer>(1,12000,
				10*(int)Math.round(DEF_SIZE*Math.sqrt(defaultRatio)/10));
		this.widthBox.setEditable(true);
		this.widthBox.setMaxWidth(Double.MAX_VALUE);
		
		this.heightBox = new Spinner<Integer>(1,12000,
				10*(int)Math.round(this.widthBox.getValue()/defaultRatio/10));
		this.heightBox.setEditable(true);
		this.widthBox.setMaxWidth(Double.MAX_VALUE);
		
		this.widthBox.getEditor().textProperty().addListener((ov, pv, nv) -> {	// link the Spinners
			if (widthBox.getEditor().textProperty().isEmpty().get())	return;
			try {
				widthBox.increment(0);	// forces the spinner to commit its value
			} catch (NumberFormatException e) {
				widthBox.getEditor().textProperty().set(pv);
				widthBox.increment(0);
			}
			if (realEdit) {
				realEdit = false;
				if (maintainRatio.isSelected()) {
					int prefHeight = (int)Math.floor(widthBox.getValue()/defaultRatio);
					heightBox.getValueFactory().setValue(prefHeight);
				}
				realEdit = true;
			}
		});
		
		this.heightBox.getEditor().textProperty().addListener((ov, pv, nv) -> {
			if (heightBox.getEditor().textProperty().isEmpty().get())	return;
			try {
				heightBox.increment(0);	// forces the spinner to commit its value
			} catch (NumberFormatException e) {
				heightBox.getEditor().textProperty().set(pv);
				heightBox.increment(0);
			}
			if (realEdit) {
				realEdit = false;
				
				if (maintainRatio.isSelected()) {
					int prefWidth = (int)Math.ceil(heightBox.getValue()*defaultRatio);
					widthBox.getValueFactory().setValue(prefWidth);
				}
				realEdit = true;
			}
		});
		
		ObservableList<String> items = FXCollections.observableArrayList(
				"None","Low","High");
		this.smoothBox = new ComboBox<String>(items);
		this.smoothBox.setValue("Low");
		this.smoothBox.setMaxWidth(Double.MAX_VALUE);
		
		this.gui = new VBox();
		this.gui.setSpacing(20);
		
		pane.contentTextProperty().addListener((arg0) -> {	// set it to refresh the gui when... the content texts?
			this.updateGUI();
		});
		this.setTitle("Save options");	// set the words on it
		pane.setHeaderText("Please configure the final image");
		pane.getButtonTypes().addAll(new ButtonType[] { ButtonType.OK, ButtonType.CANCEL });	// add buttons
		this.updateGUI();
		
		this.setResultConverter((btn) -> {	// set the return value
			ButtonData bData = btn == null ? null : btn.getButtonData();
			if (bData == ButtonData.OK_DONE) {
				ProgressBarDialog pBar = new ProgressBarDialog();
				pBar.show();
				
				return new Thread(() -> {
					Image map = parent.map(
							widthBox.getValue(),
							heightBox.getValue(),
							getSmoothing(),
							pBar);
					Platform.runLater(() -> parent.saveImage(map, pBar));
				});
			}
			else
				return null;
		});
		
		realEdit = true;
	}
	
	
	private void updateGUI() {
		this.gui.getChildren().clear();
		this.gui.getChildren().add(this.maintainRatio);
		
		GridPane grid = new GridPane();
		grid.setHgap(10);
		grid.setMaxWidth(Double.MAX_VALUE);
		grid.setAlignment(Pos.CENTER_LEFT);
		grid.getChildren().clear();
		grid.add(new Label("Width:"), 0, 0);
		grid.add(this.widthBox, 1, 0);
		grid.add(new Label("Height:"), 0, 1);
		grid.add(this.heightBox, 1, 1);
		grid.add(new Label("Smoothing:"), 0, 2);
		grid.add(this.smoothBox, 1, 2);
		this.gui.getChildren().add(grid);
		
		this.getDialogPane().setContent(this.gui);
	}
	
	
	private int getSmoothing() {
		if (smoothBox.getValue().equals("None"))
			return 1;
		else if (smoothBox.getValue().equals("Low"))
			return 2;
		else if (smoothBox.getValue().equals("High"))
			return 3;
		else
			return 0;
	}

}
