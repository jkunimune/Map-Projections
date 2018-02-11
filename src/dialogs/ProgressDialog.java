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
import javafx.concurrent.Worker;
import javafx.geometry.Pos;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonBar.ButtonData;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Dialog;
import javafx.scene.control.DialogPane;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;

public class ProgressDialog<V> extends Dialog<V> {
	
	public final double DEF_SIZE = 1000;
	
	
	private final Worker<V> worker;
	private final ProgressBar bar;
	private final Label words;
	
	
	public static <V> void trackProgress(Worker<V> worker) { //track the worker with a nicely formatted prog
		ProgressDialog<V> dialog = new ProgressDialog<V>(worker);
		worker.runningProperty().addListener((wasRunning, isRunning, observable) -> {
			if (isRunning)
				dialog.show(); //show this automatically when the Worker runs
			else
				dialog.close(); //close this automatically when the Worker stops running
		});
	}
	
	public ProgressDialog(Worker<V> worker) {
		this.worker = worker;
			
		DialogPane pane = this.getDialogPane();
		pane.setHeaderText("Please wait.");
		pane.headerTextProperty().bind(worker.messageProperty()); //the Worker's message is the header text;
		pane.getButtonTypes().add(ButtonType.CLOSE);
		((Button) pane.lookupButton(ButtonType.CLOSE)).setTooltip(new Tooltip("Please don't deactivate me!"));
		
		this.words = new Label();
		this.words.setMinWidth(50); //the words are the progress in text form;
		this.words.setAlignment(Pos.BASELINE_RIGHT);
		
		this.bar = new ProgressBar();
		this.bar.setPrefWidth(300);
		this.bar.progressProperty().addListener((old, now, observable) -> { //update the percentage text whenever the progress changes{
			if (now.doubleValue() >= 0) {
				String percentage = (Math.round(1000*now.doubleValue())/10.0)+"%";
				this.words.setText(percentage);
				this.setTitle("Please wait - "+percentage);
			}
			else {
				this.words.setText("...");
				this.setTitle("Please wait.");
			}
		});
		this.bar.progressProperty().bind(worker.progressProperty()); //the ProgressBar tracks the Worker's progress, naturally;
		
		pane.setContent(new HBox(5, bar, words));
		HBox.setHgrow(bar, Priority.ALWAYS);
		
		this.setResultConverter((btn) -> {	// set the return value
			if (btn.getButtonData() == ButtonData.CANCEL_CLOSE)
				this.worker.cancel();
			return this.worker.getValue();
		});
	}
	
	
	public void unbindAndSetProgress(double p) {
		this.bar.progressProperty().unbind();
		this.bar.setProgress(p);
	}
	
	
	public void unbindAndSetHeaderText(String text) {
		this.headerTextProperty().unbind();
		this.setHeaderText(text);
	}
	
}
