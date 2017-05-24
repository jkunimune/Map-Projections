package mapAnalyzer;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.chart.XYChart.Series;
import javafx.stage.Stage;
import mapAnalyzer.MapAnalyzer.Projection;
import util.Stat;

/**
 * An application to compare and optimize map projections
 * 
 * @author Justin Kunimune
 */
public class MapOptimizer extends Application {
	
	
	private static final String[] EXISTING_PROJECTIONS = { "Hobo-Dyer", "Winkel Tripel", "Robinson", "Van der Grinten",
			"Mercator" };
	private ScatterChart<Number, Number> chart;
	
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage stage) throws Exception {
		double[][][] globe = MapAnalyzer.globe(0.02);
		final Series<Number, Number> oldMaps = analyzeAll(globe, EXISTING_PROJECTIONS);
		final Series<Number, Number> hyperMaps = optimizeHyperelliptical();
		final Series<Number, Number> roundMaps = optimizeElliptihypercosine();
		
		chart = new ScatterChart<Number, Number>(new NumberAxis("Size distortion", 0, 1, 0.1),
				new NumberAxis("Shape distortion", 0, 0.5, 0.1));
		chart.getData().add(oldMaps);
		chart.getData().add(hyperMaps);
		chart.getData().add(roundMaps);
		
		stage.setTitle("Map Projections");
		stage.setScene(new Scene(chart));
		stage.show();
	}
	
	
	private static Series<Number, Number> analyzeAll(double[][][] points,
			String... projs) { //analyze and plot the specified preexisting map projections
		Series<Number, Number> output = new Series<Number, Number>();
		output.setName("Basic Projections");
		
		for (String name: projs)
			output.getData().add(plot(points, MapAnalyzer.projFromName(name)));
		
		return output;
	}
	
	
	private static Series<Number, Number> optimizeHyperelliptical() { //optimize and plot some hyperelliptical maps
		return new Series<Number, Number>();
	}
	
	
	private static Series<Number, Number> optimizeElliptihypercosine() { //optimize and plot some elliptical-hypebolic-cosine maps
		return new Series<Number, Number>();
	}
	
	
	private static Data<Number, Number> plot(double[][][] pts, Projection proj) {
		double[][][] distortion = MapAnalyzer.calculateDistortion(pts, proj);
		return new Data<Number, Number>(
				Stat.stdDev(distortion[0]),
				Stat.average(distortion[1]));
	}

}
