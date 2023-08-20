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

import java.util.HashMap;
import java.util.Map;

import apps.MapApplication;
import javafx.scene.Node;
import javafx.scene.control.ButtonBar.ButtonData;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Dialog;
import javafx.scene.control.DialogPane;
import javafx.scene.control.Label;
import javafx.scene.control.TreeCell;
import javafx.scene.control.TreeItem;
import javafx.scene.control.TreeView;
import javafx.scene.control.cell.TextFieldTreeCell;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.scene.text.TextFlow;
import maps.Projection;


/**
 * A dialog for perusing all of the map projections I have to offer in an ordered drop-down-type
 * thing.
 * 
 * @author jkunimune
 */
public class ProjectionSelectionDialog extends Dialog<Projection> {
	
	
	private static final double MENU_WIDTH = 250;
	private static final double TEXT_WIDTH = 300;
	
	
	private final Map<TreeItem<String>, Projection> projMap;
	private final TreeView<String> menu;
	private final TextFlow flow;
	private final GridPane text;
	
	
	
	public ProjectionSelectionDialog() {
		projMap = new HashMap<TreeItem<String>, Projection>();
		
		final TreeItem<String> root = new TreeItem<String>();
		menu = new TreeView<String>(root);
		menu.setShowRoot(false); //create and configure the TreeView of options
		menu.setPrefWidth(MENU_WIDTH);
		
		flow = new TextFlow(); //create and configure the description area
		flow.setPrefWidth(TEXT_WIDTH);
		text = new GridPane();
		text.setHgap(10);
		
		menu.getSelectionModel().selectedItemProperty().addListener((observable, old, now) -> {
				if (projMap.containsKey(now)) //selection callback to describe each projection
					describe(projMap.get(now));
				else if (now != null) {
					describe(null);
				}
			});
		menu.setCellFactory((tView) -> { //factoring cells to detect double-clicks
				final TreeCell<String> cell = new TextFieldTreeCell<String>();
				cell.setOnMouseClicked((event) -> { //on double click, close dialog
						if (event.getClickCount() >= 2 && projMap.containsKey(cell.getTreeItem())) {
							this.setResult(projMap.get(cell.getTreeItem()));
						}
					});
				return cell;
			});
		
		String[] categories = MapApplication.PROJECTION_CATEGORIES;
		Projection[][] projections = MapApplication.ALL_PROJECTIONS;
		for (int i = 0; i < categories.length; i ++) { //finally, populate the TreeView
			final TreeItem<String> header = new TreeItem<String>(categories[i]);
			root.getChildren().add(header);
			for (int j = 0; j < projections[i].length; j ++) {
				final TreeItem<String> leaf = new TreeItem<String>(projections[i][j].getName());
				projMap.put(leaf, projections[i][j]);
				header.getChildren().add(leaf);
			}
		}
		
		this.setTitle("Projection selection"); //set general properties for the dialog
		final DialogPane pane = this.getDialogPane();
		pane.setHeaderText("Choose a projection from the list below.");
		pane.getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL); //add buttons
		pane.setContent(new HBox(10, menu, new VBox(10, flow, text)));
		
		this.setResultConverter((btn) -> { //how to return a result:
				if (btn != null && btn.getButtonData() == ButtonData.OK_DONE) {
					final TreeItem<String> selection =  menu.getSelectionModel().getSelectedItem();
					return projMap.getOrDefault(selection, Projection.NULL_PROJECTION); //return the corresponding projection
				} //or NULL_PROJECTION if the user never chose anything
				else {
					return null;
				}
			});
	}
	
	
	private void describe(Projection p) {
		if (p == null) {
			flow.getChildren().clear();
			text.getChildren().clear();
			return;
		}
		
		final Text head = new Text(p.getName()+"\n");
		head.setFont(Font.font(head.getFont().getFamily(), FontWeight.BOLD, 18));
		final Text body = new Text(p.getDescription());
		flow.getChildren().setAll(head, body);
		
		text.getChildren().clear();
		text.addRow(0, new Label("Geometry:"), new Label(p.getType().getName()));
		text.addRow(1, new Label("Property:"), new Label(p.getProperty().getName()));
		text.addRow(2, new Label("Uninterrupted:"), new Label(p.isContinuous() ? "Yes" : "No"));
		text.addRow(3, new Label("Shows entire world:"), new Label(p.isFinite() ? "Yes" : "No"));
		text.addRow(4, new Label("Closed-form solution:"), new Label(p.isSolveable() ? "Yes" : "No"));
		text.addRow(5, new Label("Closed-form inverse:"), new Label(p.isInvertable() ? "Yes" : "No"));
		for (Node label: text.getChildren())
			((Label)label).setFont(body.getFont());
	}
	
}
