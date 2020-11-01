/**
 * MIT License
 * 
 * Copyright (c) 2020 Justin Kunimune
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

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

import maps.Gyorffy;
import maps.Lenticular;
import maps.Meshed;
import maps.Polyhedral;
import maps.Projection;
import maps.Pseudocylindrical;
import maps.Tobler;
import maps.WinkelTripel;

/**
 * A script to get some distortion quantiles for my paper.
 * 
 * @author Justin Kunimune
 */
public class MapEvaluator {
	
	public static final double STEP = 1e-2;
	
	public static final boolean ONLY_LAND = true;
	public static final double THRESHOLD = 3;
	public static final Projection[] PROJECTIONS = {
			Tobler.TOBLER, Pseudocylindrical.ECKERT_IV, Lenticular.WAGNER_VIII, Meshed.DANSEIJI_I,
			WinkelTripel.WINKEL_TRIPEL, Gyorffy.E, Meshed.DANSEIJI_II,
			Pseudocylindrical.HOMOLOSINE_INTERRUPTED, Meshed.DANSEIJI_III,
			Polyhedral.DYMAXION, Meshed.DANSEIJI_IV,
	};
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		BufferedImage img = ImageIO.read(new File("input/Sillouette.png")); // start by loading the mask
		boolean[][] land = new boolean[img.getHeight()][img.getWidth()];
		for (int i = 0; i < land.length; i ++)
			for (int j = 0; j < land[i].length; j ++)
				land[i][j] = (img.getRGB(j, i) & 0xFF) < 128;
		
		for (Projection projection : PROJECTIONS) {
			System.out.println(projection.getName());
			
			projection.setParameters(projection.getDefaultParameters());
			
			List<Double> arealDistortion = new ArrayList<Double>();
			List<Double> angularDistortion = new ArrayList<Double>();
			List<Double> weights = new ArrayList<Double>();
			
			double ɸStep = STEP;
			for (double ɸMin = -Math.PI/2; ɸMin < Math.PI/2; ɸMin += ɸStep) {
				double ɸMax = Math.min(ɸMin + ɸStep, Math.PI/2);
				int i = (int) ((-(ɸMin+ɸMax)/2 / Math.PI + .5) * land.length);
				
				double λStep = STEP/Math.cos((ɸMax + ɸMin)/2);
				for (double λMin = -Math.PI; λMin < Math.PI; λMin += λStep) {
					double λMax = Math.min(λMin + λStep, Math.PI);
					int j = (int) (((λMin+λMax)/2 / (2*Math.PI) + .5) * land[i].length);
					
					if (ONLY_LAND && !land[i][j]) // if we only care about land and this isn't land
						continue;
					
					double[] n = projection.project((ɸMin + 3*ɸMax)/4, (λMin + λMax)/2);
					double[] s = projection.project((3*ɸMin + ɸMax)/4, (λMin + λMax)/2);
					double[] e = projection.project((ɸMin + ɸMax)/2, (λMin + 3*λMax)/4);
					double[] w = projection.project((ɸMin + ɸMax)/2, (3*λMin + λMax)/4);
					
					double lengthScale = Math.sqrt(
							(e[0] - w[0])*(n[1] - s[1]) - (n[0] - s[0])*(e[1] - w[1]));
					double nonlinearity = Math.hypot(
							n[0] + s[0] - e[0] - w[0], n[1] + s[1] - e[1] - w[1])/lengthScale;
					if (nonlinearity < 2) { // only count distortion measurements where the map is smooth
						double dA = (Math.sin(ɸMax) - Math.sin(ɸMin))*(λMax - λMin);
						
						double dxdɸ = (n[0] - s[0])/((ɸMax - ɸMin)/2);
						double dxdλ = (e[0] - w[0])/((λMax - λMin)/2);
						double dydɸ = (n[1] - s[1])/((ɸMax - ɸMin)/2);
						double dydλ = (e[1] - w[1])/((λMax - λMin)/2);
						double dsdɸ = 1;
						double dsdλ = Math.cos((ɸMin + ɸMax)/2);
						
						double s1ps2 = Math.hypot(dxdλ/dsdλ + dydɸ/dsdɸ, dydλ/dsdλ - dxdɸ/dsdɸ);
						double s1ms2 = Math.hypot(dxdλ/dsdλ - dydɸ/dsdɸ, dydλ/dsdλ + dxdɸ/dsdɸ);
						angularDistortion.add((s1ps2 - s1ms2)/(s1ps2 + s1ms2));
						arealDistortion.add((dxdλ*dydɸ - dxdɸ*dydλ)/dsdɸ/dsdλ);
						weights.add(dA);
					}
					else {
//						System.out.printf("[%f,%f],[%f,%f],[%f,%f],[%f,%f]\n", n[0], n[1], e[0], e[1], s[0], s[1], w[0], w[1]);
//						System.out.println("Nonlinearity detected at "+ɸMin+"N, "+λMin+"E: "+nonlinearity);
					}
				}
			}
			
			double minScale = Double.POSITIVE_INFINITY, maxArea = Double.NEGATIVE_INFINITY;
			for (int i = 0; i < weights.size(); i ++) {
				if (arealDistortion.get(i) < minScale)
					minScale = arealDistortion.get(i);
				if (weights.get(i) > maxArea)
					maxArea = weights.get(i);
			}
			double shearedArea = 0, swollenArea = 0, totalArea = 0;
			for (int i = 0; i < weights.size(); i ++) {
				if (absoluteQuality(angularDistortion.get(i)) > THRESHOLD)
					shearedArea += weights.get(i);
				if (arealDistortion.get(i) - minScale > THRESHOLD)
					swollenArea += weights.get(i);
				totalArea += weights.get(i);
			}
			
			System.out.printf("Usable area: %.6f sr\n", totalArea);
			System.out.printf("Anglular:   %8.4f %%\n", 100*shearedArea/totalArea);
			System.out.printf("Areal:      %8.4f %%\n", 100*swollenArea/totalArea);
			System.out.printf("Precision:  %8.4f %%\n", 100*maxArea/totalArea);
			System.out.println();
		}
	}
	
	private static final double absoluteQuality(double x) {
		if (Math.abs(x) < 1)
			return 1/x;
		else
			return x;
	}

}
