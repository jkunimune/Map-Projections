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

import static java.lang.Math.round;
import static java.lang.Math.sqrt;

public class MapConfigurationDialog extends Dialog<Boolean> {

	public final int MIN_SIZE = 16, MAX_SIZE = 60000;
	public final double DEF_SIZE = 1200;
	
	
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
		
		this.widthBox = new Spinner<Integer>(MIN_SIZE, MAX_SIZE,
				10*(int)round(DEF_SIZE*sqrt(defaultRatio)/10));
		this.widthBox.setEditable(true);
		this.widthBox.setMaxWidth(Double.MAX_VALUE);
		
		this.heightBox = new Spinner<Integer>(MIN_SIZE, MAX_SIZE,
				(int)round(this.widthBox.getValue()/defaultRatio));
		this.heightBox.setEditable(true);
		this.widthBox.setMaxWidth(Double.MAX_VALUE);
		
		this.widthBox.valueProperty().addListener((observable, prev, now) -> {	// link the Spinners
				if (realEdit && maintainRatio.isSelected()) {
					realEdit = false;
					int prefHeight = (int)round(widthBox.getValue()/defaultRatio);
					heightBox.getValueFactory().setValue(prefHeight);
					realEdit = true;
				}
			});
		this.heightBox.valueProperty().addListener((observable, prev, now) -> {
				if (realEdit && maintainRatio.isSelected()) {
					realEdit = false;
					int prefWidth = (int)round(heightBox.getValue()*defaultRatio);
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
		
		ObservableList<String> items = FXCollections.observableArrayList("None","Low","High");
		this.smoothBox = new ComboBox<String>(items);
		this.smoothBox.setValue("Low");
		this.smoothBox.setMaxWidth(Double.MAX_VALUE);
		
		this.gui = new VBox(20);
		
		pane.contentTextProperty().addListener((arg0) -> this.updateGUI());	// set it to refresh the gui when... the content texts?
		this.setTitle("Save options");	// set the words on it
		pane.setHeaderText("Please configure the final image");
		pane.getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);	// add buttons
		this.updateGUI();
		
		this.setResultConverter(
				(btn) -> (btn != null && btn.getButtonData() == ButtonData.OK_DONE));
		
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
		grid.addRow(0, new Label("Width:"), this.widthBox);
		grid.addRow(1, new Label("Height:"), this.heightBox);
		grid.addRow(2, new Label("Smoothing:"), this.smoothBox);
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
