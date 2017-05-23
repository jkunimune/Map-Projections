package mapAnalyzer;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart.Series;
import javafx.stage.Stage;

public class MapOptimizer extends Application {
	
	
	private static final String[] EXISTING_PROJECTIONS = { "Hobo-Dyer", "Winkel Tripel", "Robinson", "Van der Grinten",
			"Mercator" };
	private ScatterChart<Number, Number> chart;
	
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage stage) throws Exception {
		final Series<Number, Number> oldMaps = analyzeAll(EXISTING_PROJECTIONS);
		final Series<Number, Number> hyperMaps = optimizeHyperelliptical();
		final Series<Number, Number> roundMaps = optimizeElliptihypercosine();
		
		chart = new ScatterChart<Number, Number>(new NumberAxis("Size distortion", 0, 2, 0.5),
				new NumberAxis("Shape distortion", 0, Math.PI / 2, 0.5));
		chart.getData().add(oldMaps);
		chart.getData().add(hyperMaps);
		chart.getData().add(roundMaps);
		
		stage.setTitle("Map Projections");
		stage.setScene(new Scene(chart));
		stage.show();
	}
	
	
	private static Series<Number, Number> analyzeAll(String... projs) { //analyze and plot the specified preexisting map projections
		return new Series<Number, Number>();
	}
	
	
	private static Series<Number, Number> optimizeHyperelliptical() { //optimize and plot some hyperelliptical maps
		return new Series<Number, Number>();
	}
	
	
	private static Series<Number, Number> optimizeElliptihypercosine() { //optimize and plot some elliptical-hypebolic-cosine maps
		return new Series<Number, Number>();
	}
	
}
