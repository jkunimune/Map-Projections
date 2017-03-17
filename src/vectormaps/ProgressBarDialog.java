package vectormaps;
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
		
		pane.contentTextProperty().addListener((arg0) -> {	// set it to refresh the gui when... the content texts?
			this.resetBar();
		});
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
		this.bar.setProgress(p);
		Platform.runLater(() ->
				this.words.setText((Math.round(1000*p)/10.0)+"%"));
	}

}
