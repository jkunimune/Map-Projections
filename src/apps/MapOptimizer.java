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
import java.util.function.Function;

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
import maps.Arbitrary;
import maps.Cylindrical;
import maps.Misc;
import maps.Polyhedral;
import maps.Projection;
import maps.Tobler;
import maps.WinkelTripel;
import utils.linalg.Matrix;
import utils.linalg.Vector;

/**
 * An application to compare and optimize map projections
 * 
 * @author Justin Kunimune
 */
public class MapOptimizer extends Application {
	
	private static final Projection[] EXISTING_PROJECTIONS = { Cylindrical.BEHRMANN,
			Arbitrary.ROBINSON, Cylindrical.PLATE_CARREE, Cylindrical.GALL_STEREOGRAPHIC,
			Misc.PEIRCE_QUINCUNCIAL };
	private static final Projection[] PROJECTIONS_TO_OPTIMIZE = { Tobler.TOBLER,
			WinkelTripel.WINKEL_TRIPEL, Polyhedral.TETRAPOWER, Polyhedral.AUTHAPOWER };
	private static final double[] WEIGHTS = { 0., .125, .25, .375, .5, .625, .75, .875, 1. };
	private static final int NUM_BRUTE_FORCE = 30;
	private static final int NUM_BFGS_ITERATE = 7;
	private static final double GOLDSTEIN_C = 0.5;
	private static final double BACKTRACK_TAU = 0.5;
	private static final double BACKTRACK_ALF0 = 4;
	private static final double DEL_X = 0.05;
	private LineChart<Number, Number> chart;
	
	private static final double[][][] GLOBE = Projection.hemisphere(0.01);
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage stage) throws Exception {
		final long startTime = System.currentTimeMillis();
		
		chart = new LineChart<Number, Number>(new NumberAxis("Size distortion", 0, .6, 0.1),
				new NumberAxis("Shape distortion", 0, .6, 0.1));
		chart.setCreateSymbols(true);
		chart.setAxisSortingPolicy(SortingPolicy.NONE);
		
		PrintStream log = new PrintStream(new File("output/parameters.txt"));
		
		chart.getData().add(analyzeAll(EXISTING_PROJECTIONS));
		for (Projection p: PROJECTIONS_TO_OPTIMIZE)
			chart.getData().add(optimiseFamily(p, log));
		
		System.out.println("Total time elapsed: " + (System.currentTimeMillis() - startTime) / 60000. + "m");
		
		stage.setTitle("Map Projections");
		stage.setScene(new Scene(chart));
		
		ImageIO.write(SwingFXUtils.fromFXImage(chart.snapshot(null, null), null), "png",
				new File("output/graph - optimizer.png"));
		
		stage.show();
		log.close();
	}
	
	
	private static Series<Number, Number> analyzeAll(Projection... projs) { //analyze and plot the specified preexisting map projections.
		System.out.println("Analyzing " + Arrays.toString(projs));
		Series<Number, Number> output = new Series<Number, Number>(); //These projections must not be parametrized
		output.setName("Basic Projections");
		
		for (Projection proj : projs)
			if (!proj.isParametrized())
				output.getData().add(plotDistortion(proj, new double[0]));
		
		return output;
	}
	
	
	private static Data<Number, Number> plotDistortion(Projection proj, double[] params) {
		double[] distortion = proj.avgDistortion(GLOBE, proj.getDefaultParameters());
		return new Data<Number, Number>(distortion[0], distortion[1]);
	}
	
	
	private static final double weighDistortion(double[] distortions, double weight) {
		return distortions[0]*weight + distortions[1]*(1-weight);
	}
	
	
	private static Series<Number, Number> optimiseFamily(
			Projection proj, PrintStream log) { //optimize and plot some maps of a given family
		System.out.println("Optimizing " + proj.getName());
		
		final double[][] best = new double[WEIGHTS.length][proj.getNumParameters()];
		
		for (int k = 0; k < WEIGHTS.length; k ++) {
			final double weighFactor = WEIGHTS[k];
			double[] currentBest = bruteForceMinimise(
					(params) -> weighDistortion(proj.avgDistortion(GLOBE, params), weighFactor),
					proj.getParameterValues());
			best[k] = bfgsMinimise(
					(params) -> weighDistortion(proj.avgDistortion(GLOBE, params), weighFactor),
					currentBest);
		}
		
		final Series<Number, Number> output = new Series<Number, Number>();
		output.setName(proj.getName());
		
		log.println("We got the best " + proj.getName() + " projections using:"); //now log it
		for (double[] bestForWeight : best) { //for each weight
			log.print("\t");
			
			for (int i = 0; i < proj.getNumParameters(); i++)
				log.print("t" + i + "=" + bestForWeight[i] + "; "); //print the parameters used
			
			double[] distortion = proj.avgDistortion(GLOBE, bestForWeight);
			log.println("\t(" + distortion[0] + ", " + distortion[1] + ")"); //print the resulting distortion
			
			output.getData().add(new Data<Number, Number>(distortion[0], distortion[1])); //plot it
		}
		log.println();
		return output;
	}
	
	
	/**
	 * Returns the parameters that minimise the function, based on a simple brute-force
	 * parameter sweep.
	 * @param func - The function to minimise.
	 * @param bounds - Parameter limits for each argument.
	 * @return An array containing the best input to func that it found.
	 */
	private static double[] bruteForceMinimise(Function<double[], Double> func, double[][] bounds) {
		System.out.println("BF = [");
		final double[] params = new double[bounds.length];
		for (int i = 0; i < params.length; i++)
			params[i] = bounds[i][0]; // initialize params
		double bestValue = Double.POSITIVE_INFINITY;
		double[] bestParams = new double[params.length];
		
		while (true) { // run until you've exhausted the parameter space
			double avgDist = func.apply(params);
			if (avgDist < bestValue) {
				bestValue = avgDist;
				bestParams = params.clone();
			}
			for (int i = 0; i < params.length; i ++)
				System.out.print(params[i]+", ");
			System.out.println(avgDist+";");
			
			int i;
			for (i = 0; i <= params.length; i++) { // iterate the parameters
				if (i == params.length) {
					System.out.println("];");
					return bestParams; // if you made it through all the parameters without breaking, you're done!
				}
				
				final double step = (bounds[i][1] - bounds[i][0]) /
						Math.floor(Math.pow(NUM_BRUTE_FORCE, 1./params.length));
				if (params[i] + step < bounds[i][1] + 1e-5) {
					for (int j = 0; j < i; j ++)
						params[j] = bounds[j][0];
					params[i] += step;
					break;
				}
			}
		}
	}
	
	
	/**
	 * Calculates the set of parameters that minimises the function using BFGS optimisation with
	 * a backtracking line search.
	 * @param arrFunction - The function that takes a parameter array and returns a double value.
	 * @param x0 - The initial guess.
	 * @return The array of parameters that mimimise arrFunction.
	 */
	private static double[] bfgsMinimise(Function<double[], Double> arrFunction, double[] x0) { //The Broyden-Fletcher-Goldfarb-Shanno algorithm
		System.out.println("BFGS = [");
		final int n = x0.length;
		final Matrix I = Matrix.identity(n);
		final Function<Vector, Double> func = (vec) -> arrFunction.apply(vec.asArray());
		
		Vector xk = new Vector(x0); //initial variable values
		double fxk = func.apply(xk);
		Matrix Binv = hessian(func, xk, fxk).inverse();
		Vector gradFxk = grad(func, xk, fxk); //function at current location
		
		for (int k = 0; k < NUM_BFGS_ITERATE; k ++) { //(I'm not sure how to test for convergence here, so I'm just running a set number of iterations)
			Vector pk = Vector.fromMatrix(Binv.times(gradFxk)); //apply Newton's method for initial step direction
			pk = pk.times(-Math.signum(pk.dot(gradFxk))); //but make sure it points downhill
			
			double alfk = BACKTRACK_ALF0; //perform a backtracking line search to find the best alpha
			double fxkp1 = func.apply(xk.plus(pk.times(alfk)));
			while ((!Double.isFinite(fxkp1) || fxkp1 > fxk + alfk*pk.dot(gradFxk)*GOLDSTEIN_C)) {
				if (alfk <= 1e-5)
					return xk.asArray(); //a simple way to check for convergence: if xk gets ridiculously small, we're done here.
				alfk *= BACKTRACK_TAU;
				fxkp1 = func.apply(xk.plus(pk.times(alfk)));
			}
			
			Vector sk = pk.times(alfk); //iterate
			Vector xkp1 = xk.plus(sk);
			
			Vector gradFxkp1 = grad(func, xkp1, fxkp1); //compute new gradient
			Vector yk = gradFxkp1.minus(gradFxk); //and gradient change
			
			Matrix a = I.minus(sk.times(yk.T()).times(1/yk.dot(sk)));
			Matrix b = sk.times(sk.T()).times(1/yk.dot(sk));
			Binv = a.times(Binv).times(a.T()).plus(b); //update Binv
			
			xk = xkp1;
			fxk = fxkp1;
			gradFxk = gradFxkp1; //and save the gradient
		}
		System.out.println("];");
		return xk.asArray();
	}
	
	
	/**
	 * Calculates the gradient vector of f at x.
	 * @param f - The function to differentiate.
	 * @param x - The point at which to differentiate.
	 * @param fx - The value of f(x), to speed computations.
	 * @return The vector of partial derivatives of f at x.
	 */
	private static Vector grad(Function<Vector, Double> f, Vector x, double fx) {
		final int n = x.getLength();
		Vector gradF = new Vector(n); //compute the gradient
		
		for (int i = 0; i < n; i ++) {
			Vector xph = x.plus(Vector.unit(i,n).times(DEL_X));
			double fxph = f.apply(xph);
			gradF.setElement(i, (fxph-fx)/DEL_X);
		}
		for (double d: x.asArray())
			System.out.print(d+", ");
		System.out.println(fx+";");
		return gradF;
	}
	
	
	/**
	 * Computes the Hessian matrix of f at x.
	 * @param f - The function to differentiate.
	 * @param x - The point at which to differentiate.
	 * @param fx - The value of f(x), to aid in computation.
	 * @return The Jacobian of the gradient, a symmetric Matrix of second derivatives.
	 */
	private static Matrix hessian(Function<Vector, Double> f, Vector x, double fx) {
		final int n = x.getLength();
		
		double[] values = new double[(int)Math.pow(3, n-1)*2+1]; //points in array placed with ternary coordinates
		values[0] = fx;
		for (int i = 0; i < n; i ++) { //for each primary dimension
			for (int j = i; j < n; j ++) { //for each secondary dimension (skip a few to prevent redundant calculations)
				int k = (int)Math.pow(3, i) + (int)Math.pow(3, j); //calculate the ternary index
				Vector dx = Vector.unit(i, n).plus(Vector.unit(j, n)).times(DEL_X); //go a bit in both directions
				values[k] = f.apply(x.plus(dx)); //calculate and save
			}
			int k = (int)Math.pow(3, i); //do the same with just i, no j
			Vector dx = Vector.unit(i, n).times(DEL_X);
			values[k] = f.apply(x.plus(dx));
		}
		
		Matrix h = new Matrix(n, n);
		for (int i = 0; i < n; i ++) { //compute the derivatives and fill the matrix
			for (int j = i; j < n; j ++) {
				int dxi = (int)Math.pow(3, i);
				int dxj = (int)Math.pow(3, j);
				double dfdx0 = (values[dxi] - values[0])/DEL_X;
				double dfdx1 = (values[dxi+dxj] - values[dxj])/DEL_X;
				double d2fdx2 = (dfdx1 - dfdx0)/DEL_X;
				h.setElement(i, j, d2fdx2);
				h.setElement(j, i, d2fdx2);
			}
		}
		
		return h;
	}
	
}
