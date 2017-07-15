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
package apps;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import javax.imageio.ImageIO;

import javafx.application.Application;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.LineChart.SortingPolicy;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart.Data;
import javafx.scene.chart.XYChart.Series;
import javafx.stage.Stage;
import maps.Projection;

/**
 * An application to compare and optimize map projections
 * 
 * @author Justin Kunimune
 */
public class MapOptimizer extends Application {

	private static final Projection[] EXISTING_PROJECTIONS = { Projection.HOBO_DYER, Projection.ROBINSON,
			Projection.PLATE_CARREE, Projection.PEIRCE_QUINCUNCIAL };
	private static final double[] WEIGHTS = { .083, .20, .33, .50, .71, 1.0, 1.4, 2.0, 3.0, 5.0, 12. };
	private static final int NUM_DESCENT = 40;
	private LineChart<Number, Number> chart;
	
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage stage) throws Exception {
		final long startTime = System.currentTimeMillis();
		
		chart = new LineChart<Number, Number>(
				new NumberAxis("Size distortion", 0, 1, 0.2),
				new NumberAxis("Shape distortion", 0, 1, 0.2));
		chart.setCreateSymbols(true);
		chart.setAxisSortingPolicy(SortingPolicy.NONE);
		double[][][] globe = Projection.globe(0.01);
		PrintStream log = new PrintStream(new File("output/parameters.txt"));
		
		chart.getData().add(analyzeAll(globe, EXISTING_PROJECTIONS));
//		chart.getData().add(optimizeFamily(Projection.WINKEL_TRIPEL, globe, log));
//		chart.getData().add(optimizeFamily(Projection.TOBLER, globe, log));
//		chart.getData().add(optimizeFamily(Projection.HYPERELLIPOWER, globe, log));
//		chart.getData().add(optimizeFamily(Projection.TETRAPOWER, globe, log));
		chart.getData().add(optimizeFamily(Projection.TETRAFILLET, globe, log));
		
		System.out.println("Total time elapsed: "+
				(System.currentTimeMillis()-startTime)/1000.+"s");
		
		stage.setTitle("Map Projections");
		stage.setScene(new Scene(chart));
		
		ImageIO.write(
				SwingFXUtils.fromFXImage(chart.snapshot(null, null), null),
				"png", new File("output/graph - optimizer.png"));
		
		stage.show();
		log.close();
	}
	
	
	private static Series<Number, Number> analyzeAll(double[][][] points,
			Projection... projs) { //analyze and plot the specified preexisting map projections. These projections must not be parametrized
		System.out.println("Analyzing "+Arrays.toString(projs));
		Series<Number, Number> output = new Series<Number, Number>();
		output.setName("Basic Projections");
		
		for (Projection proj: projs)
			if (!proj.isParametrized())
				output.getData().add(plotDistortion(points, proj, new double[0]));
		
		return output;
	}
	
	
	private static Series<Number, Number> optimizeFamily(
			Projection proj, double[][][] points, PrintStream log) { // optimize and plot some maps of a given family maps
		System.out.println("Optimizing "+proj.getName());
		final double[][] currentBest = new double[WEIGHTS.length][3+proj.getNumParameters()]; //the 0-3 cols are the min distortions for each weight, the other cols are the values of k and n that caused that
		for (int k = 0; k < WEIGHTS.length; k ++)
			currentBest[k][0] = Integer.MAX_VALUE;
		
		final double[][] bounds = proj.getParameterValues();
		final double[] params = new double[proj.getNumParameters()];
		for (int i = 0; i < params.length; i ++) 	params[i] = bounds[i][0]; //initialize params
		
		while (true) { //start with brute force
			double[] distortions = proj.avgDistortion(points, params);
			System.out.println(Arrays.toString(params)+": "+Arrays.toString(distortions));
			for (int k = 0; k < WEIGHTS.length; k ++) {
				final double avgDist = weighDistortion(WEIGHTS[k], distortions);
				if (avgDist < currentBest[k][0]) {
					currentBest[k][0] = avgDist;
					currentBest[k][1] = distortions[0];
					currentBest[k][2] = distortions[1];
					System.arraycopy(params,0, currentBest[k],3, params.length);
				}
			}
			
			int i;
			for (i = 0; i < params.length; i ++) { //iterate the parameters
				final double step = (bounds[i][1]-bounds[i][0])/Math.floor(Math.pow(16, 1./params.length));
				if (params[i]+step < bounds[i][1]+1e-5) {
					for (int j = 0; j < i; j ++)
						params[j] = bounds[j][0];
					params[i] += step;
					break;
				}
			}
			if (i == params.length)	break; //if you made it through the for loop without breaking (finding a parameter to increment), you're done!
		}
		
		final double h = 1e-7;
		final double delX = 5e-2;
		for (int k = 0; k < WEIGHTS.length; k ++) { //now do gradient descent
			System.arraycopy(currentBest[k],3, params,0, params.length);
			System.out.println("Starting gradient descent with weight "+WEIGHTS[k]+" and initial parameters "+Arrays.toString(params));
			double fr0 = currentBest[k][0];
			double[] frd = new double[params.length];
			for (int i = 0; i < NUM_DESCENT; i ++) {
				if (i > 0)
					fr0 = weighDistortion(WEIGHTS[k],
							proj.avgDistortion(points, params)); //calculate the distortion here
				System.out.println(Arrays.toString(params)+" -> "+fr0);
				for (int j = 0; j < params.length; j ++) {
					params[j] += h;
					frd[j] = weighDistortion(WEIGHTS[k],
							proj.avgDistortion(points, params)); //and the distortion nearby
					params[j] -= h;
				}
				for (int j = 0; j < params.length; j ++)
					params[j] -= (frd[j]-fr0)/h*delX; //use that to approximate the gradient and go in that direction
			}
			System.arraycopy(params,0, currentBest[k],3, params.length);
			System.arraycopy(proj.avgDistortion(points, params), 0,
					currentBest[k], 1, 2);
		}
		
		final Series<Number, Number> output = new Series<Number, Number>();
		output.setName(proj.getName());
		
		log.println("We got the best "+proj.getName()+" projections using:");
		for (double[] best: currentBest) {
			log.print("\t");
			for (int i = 0; i < params.length; i ++)
				log.print("t"+i+"="+best[3+i]+"; ");
			log.println("\t("+best[1]+", "+best[2]+")");
			output.getData().add(new Data<Number, Number>(best[1], best[2]));
		}
		log.println();
		return output;
	}
	
	
	private static final double weighDistortion(double weight, double... distortions) {
		return Math.pow(distortions[0], 1.5) +
				weight*Math.pow(distortions[1], 1.5);
	}
	
	
	private static Data<Number, Number> plotDistortion(
			double[][][] pts, Projection proj, double[] params) {
		double[] distortion = proj.avgDistortion(pts, params);
		return new Data<Number, Number>(distortion[0], distortion[1]);
	}

}
