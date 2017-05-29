package mapAnalyzer;

import java.util.Arrays;
import java.util.function.BinaryOperator;

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
			"Mercator" }; //"Pierce Quincuncial" };
	private static final double[] WEIGHTS = { .25, .43, .67, 1.0, 1.5, 2.3, 4.0 };
	private ScatterChart<Number, Number> chart;
	
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage stage) throws Exception {
		final long startTime = System.currentTimeMillis();
		
		chart = new ScatterChart<Number, Number>(
				new NumberAxis("Size distortion", 0, 1, 0.2),
				new NumberAxis("Shape distortion", 0, 1, 0.2));
		double[][][] globe = MapAnalyzer.globe(0.02);
		
		chart.getData().add(analyzeAll(globe, EXISTING_PROJECTIONS));
		chart.getData().add(optimizeHyperelliptical(globe));
		chart.getData().add(optimizeElliptabola(globe));
		
		System.out.println("Total time elapsed: "+
				(System.currentTimeMillis()-startTime)/1000.+"s");
		
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
	
	
	private static Series<Number, Number> optimizeFamily(BinaryOperator<double[]> projectionFam, String name,
			double[][] bounds, double[][][] points) { // optimize and plot some maps of a given family maps
		final double[][] currentBest = new double[WEIGHTS.length][3+bounds.length]; //the 0-3 cols are the min distortions for each weight, the other cols are the values of k and n that caused that
		for (int k = 0; k < WEIGHTS.length; k ++)
			currentBest[k][0] = Integer.MAX_VALUE;
		
		final double[] params = new double[bounds.length];
		for (int i = 0; i < params.length; i ++)	params[i] = bounds[i][0]; //initialize params
		
		while (true) {
			int i;
			for (i = 0; i < params.length; i ++) { //iterate the parameters
				if (params[i] < bounds[i][1]) {
					for (int j = 0; j < i; j ++)
						params[j] = bounds[j][0];
					params[i] += bounds[i][2];
					break;
				}
			}
			if (i == params.length)	break; //if you made it through the for loop without breaking (finding a parameter to increment), you're done!
			System.out.println(Arrays.toString(params));
			
			double[][][] dist = MapAnalyzer.calculateDistortion(points,
					(coords) -> projectionFam.apply(coords, params));
			final double avgSizeD = Stat.stdDev(dist[0]);
			final double avgShapeD = Stat.mean(dist[1]);
			for (int k = 0; k < WEIGHTS.length; k ++) {
				final double avgDist = Math.pow(avgSizeD,1.5) + WEIGHTS[k]*Math.pow(avgShapeD,1.5);
				if (avgDist < currentBest[k][0]) {
					currentBest[k][0] = avgDist;
					currentBest[k][1] = avgSizeD;
					currentBest[k][2] = avgShapeD;
					System.arraycopy(params,0, currentBest[k],3,params.length);
				}
			}
		}
		
		final Series<Number, Number> output = new Series<Number, Number>();
		output.setName(name);
		
		System.out.println("We got the best "+name+" projections using:");
		for (double[] best: currentBest) {
			System.out.print("\t");
			for (int i = 0; i < params.length; i ++)
				System.out.print("t"+i+"="+best[3+i]+"; ");
			System.out.println();
			output.getData().add(new Data<Number, Number>(best[1], best[2]));
		}
		return output;
	}
	
	
	private static Series<Number, Number> optimizeHyperelliptical(double[][][] points) { //optimize and plot some hyperelliptical maps
		return optimizeFamily(MapOptimizer::hyperelliptical, "Hyperelliptic",
				new double[][] {{2.5,5,.125}, {0.5,1.75,.125}, {1.0,1.5,.0625}}, points);
	}
	
	
	private static Series<Number, Number> optimizeElliptabola(double[][][] points) { //optimize and plot some elliptical-hypebolic-cosine maps
		return optimizeFamily(MapOptimizer::elliptabola, "Elliptabola",
				new double[][] {{.25,1,.125}, {.75,1.75,.125}, {1.25,2,.125}, {.25,1.5,.125}}, points);
	}
	
	
	private static Data<Number, Number> plot(double[][][] pts, Projection proj) {
		double[][][] distortion = MapAnalyzer.calculateDistortion(pts, proj);
		return new Data<Number, Number>(
				Stat.stdDev(distortion[0]),
				Stat.mean(distortion[1]));
	}
	
	
	private static final double[] hyperelliptical(double[] coords, double[] params) { //a hyperelliptic map projection with hyperellipse order k and lattitudinal spacind described by x^n/n
		final double lat = coords[0], lon = coords[1];
		final double k = params[0], n = params[1], a = params[2];
		
		return new double[] {
				Math.pow(1 - Math.pow(Math.abs(lat/(Math.PI/2)), k),1/k)*lon,
				(1-Math.pow(1-Math.abs(lat/(Math.PI/2)), n))/Math.sqrt(n)*Math.signum(lat)*Math.PI/2*a};
	}
	
	
	private static double[] elliptabola(double[] coords, double[] params) {
		final double lat = coords[0], lon = coords[1];
		final double k1 = params[0], k2 = params[1], a = params[2], p2 = params[3];
		
		final double f1p2 = Math.pow(1-Math.pow(1-Math.abs(lat)/(Math.PI/2), k1),2);
		final double f2l2 = Math.pow(1-Math.pow(1-Math.abs(lon)/Math.PI, k2),2);
		final double ax = f1p2/(p2*p2);
		final double bx = 2*f1p2/(p2)+1/(f2l2);
		final double cx = f1p2 - 1;
		final double ay = f2l2;
		final double by = p2/Math.sqrt(f1p2);
		final double cy = -p2 - f2l2;
		return new double[] {
				Math.PI*Math.sqrt((-bx+Math.sqrt(bx*bx-4*ax*cx))/(2*ax))*Math.signum(lon),
				Math.PI*(-by+Math.sqrt(by*by-4*ay*cy))/(2*ay)*Math.signum(lat)/a
		};
	}

}
