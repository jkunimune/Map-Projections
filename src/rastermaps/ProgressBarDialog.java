package rastermaps;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Dialog;
import javafx.scene.control.DialogPane;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.Tooltip;

public class ProgressBarDialog extends Dialog<Void> {

	public final double DEF_SIZE = 1000;
	
	
	private final ProgressBar bar;
	
	public ProgressBarDialog() {
		DialogPane pane = this.getDialogPane();
		
		this.bar = new ProgressBar();
		this.bar.setPrefWidth(300);
		this.bar.setTooltip(new Tooltip("Oh, did you want a cancel button? Well, TOO BAD."));
		
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
		this.getDialogPane().setContent(this.bar);
	}
	
	
	public void setProgress(double p) {
		this.bar.setProgress(p);
	}

}
