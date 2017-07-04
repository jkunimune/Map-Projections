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
import javafx.application.Platform;
import javafx.geometry.Pos;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Dialog;
import javafx.scene.control.DialogPane;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;

public class ProgressBarDialog extends Dialog<Void> {

	public final double DEF_SIZE = 1000;
	
	
	private final HBox box;
	private final ProgressBar bar;
	private final Label words;
	
	public ProgressBarDialog() {
		DialogPane pane = this.getDialogPane();
		
		this.bar = new ProgressBar();
		this.bar.setPrefWidth(300);
		this.bar.setTooltip(new Tooltip("Oh, did you want a cancel button? Well, TOO BAD."));
		
		this.words = new Label();
		this.words.setMinWidth(40);
		this.words.setAlignment(Pos.BASELINE_RIGHT);
		
		this.box = new HBox(5);
		
		this.setTitle("");	// set the words on it
		pane.setHeaderText("Please wait.");
		pane.getButtonTypes().addAll(new ButtonType[] { ButtonType.CLOSE });
		((Button) pane.lookupButton(ButtonType.CLOSE)).setText("Run in background");	// you can't close it
		this.resetBar();
		
		this.setResultConverter((btn) -> {	// set the return value
			return null;
		});
	}
	
	
	private void resetBar() {
		this.bar.setProgress(0);
		this.words.setText("0.0%");
		this.box.getChildren().clear();
		this.box.getChildren().add(this.bar);
		this.box.getChildren().add(this.words);
		HBox.setHgrow(this.bar, Priority.ALWAYS);
		this.getDialogPane().setContent(this.box);
	}
	
	
	public void setProgress(double p) {
		Platform.runLater(() -> {
				if (p >= 0)
					this.words.setText((Math.round(1000*p)/10.0)+"%");
				else {
					this.bar.setProgress(-1);
					this.words.setText("...");
				}
			});
		if (p >= 0)
			this.bar.setProgress(p);
	}

}
