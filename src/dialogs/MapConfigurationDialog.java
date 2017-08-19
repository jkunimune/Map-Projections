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
package dialogs;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Pos;
import javafx.scene.control.ButtonBar.ButtonData;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Dialog;
import javafx.scene.control.DialogPane;
import javafx.scene.control.Label;
import javafx.scene.control.Spinner;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.VBox;

public class MapConfigurationDialog extends Dialog<Boolean> {

	public final double DEF_SIZE = 1000;
	
	
	private final double defaultRatio;
	private final VBox gui;
	private final CheckBox maintainRatio;
	private final Spinner<Integer> widthBox, heightBox;
	private final ComboBox<String> smoothBox;
	
	private boolean realEdit; // is a real edit happening, or is it just me?
	
	public MapConfigurationDialog(double defAsp) {
		this.defaultRatio = defAsp;
		DialogPane pane = this.getDialogPane();
		
		this.maintainRatio = new CheckBox("Maintain aspect ratio");	// instantiate the components
		this.maintainRatio.setSelected(true);
		
		this.widthBox = new Spinner<Integer>(5,10000,
				10*(int)Math.round(DEF_SIZE*Math.sqrt(defaultRatio)/10));
		this.widthBox.setEditable(true);
		this.widthBox.setMaxWidth(Double.MAX_VALUE);
		
		this.heightBox = new Spinner<Integer>(5,10000,
				10*(int)Math.round(this.widthBox.getValue()/defaultRatio/10));
		this.heightBox.setEditable(true);
		this.widthBox.setMaxWidth(Double.MAX_VALUE);
		
		this.widthBox.valueProperty().addListener((observable, prev, now) -> {	// link the Spinners
				if (realEdit && maintainRatio.isSelected()) {
					realEdit = false;
					int prefHeight = (int)Math.round(widthBox.getValue()/defaultRatio);
					heightBox.getValueFactory().setValue(prefHeight);
					realEdit = true;
				}
			});
		this.heightBox.valueProperty().addListener((observable, prev, now) -> {
				if (realEdit && maintainRatio.isSelected()) {
					realEdit = false;
					int prefWidth = (int)Math.round(heightBox.getValue()*defaultRatio);
					widthBox.getValueFactory().setValue(prefWidth);
					realEdit = true;
				}
			});
		
		this.widthBox.focusedProperty().addListener((observable, prev, now) -> { //make the spinners commit their
				if (!now) 	widthBox.increment(0);
			});
		this.heightBox.focusedProperty().addListener((observable, prev,now) -> { //values when focus is lost
				if (!now) 	heightBox.increment(0);
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
		
		this.setResultConverter((btn) -> {
				if (btn != null && btn.getButtonData() == ButtonData.OK_DONE)
					return true;
				else
					return false;
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
	
	
	public int[] getDims() {
		return new int[] { widthBox.getValue(), heightBox.getValue() };
	}
	
	
	public int getSmoothing() {
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
