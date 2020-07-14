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

import maps.Cylindrical;
import maps.Gyorffy;
import maps.Lenticular;
import maps.Meshed;
import maps.Polyhedral;
import maps.Projection;
import maps.Pseudocylindrical;
import maps.WinkelTripel;

/**
 * A script to get some distortion quantiles for my paper.
 * 
 * @author Justin Kunimune
 */
public class MapEvaluator {
	
	public static final double STEP = 1e-2;
	
	public static final double THRESHOLD = 2.0;
	public static final Projection[] PROJECTIONS = {
			Pseudocylindrical.SINUSOIDAL, Lenticular.VAN_DER_GRINTEN, Cylindrical.MERCATOR,
//			Pseudocylindrical.MOLLWEIDE, Pseudocylindrical.ECKERT_IV,
//			WinkelTripel.WINKEL_TRIPEL, Gyorffy.E, Gyorffy.F, Polyhedral.DYMAXION,
//			Meshed.DANSEIJI_I, Meshed.DANSEIJI_II, Meshed.DANSEIJI_III,
	};

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		for (Projection projection : PROJECTIONS) {
			System.out.println(projection.getName());
			
			double totalArea = 0;
			double shearedArea = 0, swollenArea = 0, unitArea = 0;
			double ɸStep = STEP;
			for (double ɸMin = -Math.PI/2; ɸMin < Math.PI/2; ɸMin += ɸStep) {
				double ɸMax = Math.min(ɸMin + ɸStep, Math.PI/2);
				double λStep = STEP/Math.cos((ɸMax + ɸMin)/2);
				for (double λMin = -Math.PI; λMin < Math.PI; λMin += λStep) {
					double λMax = Math.min(λMin + λStep, Math.PI);
					
					double dA = (Math.sin(ɸMax) - Math.sin(ɸMin))*(λMax - λMin);
					
					double[] n = projection.project((ɸMin + 3*ɸMax)/4, (λMin + λMax)/2);
					double[] s = projection.project((3*ɸMin + ɸMax)/4, (λMin + λMax)/2);
					double[] e = projection.project((ɸMin + ɸMax)/2, (λMin + 3*λMax)/4);
					double[] w = projection.project((ɸMin + ɸMax)/2, (3*λMin + λMax)/4);
					
					double lengthScale = Math.sqrt(Math.abs(
							(e[0] - w[0])*(n[1] - s[1]) - (n[0] - s[0])*(e[1] - w[1])));
					double nonlinearity = Math.hypot(
							n[0] + s[0] - e[0] - w[0], n[1] + s[1] - e[1] - w[1])/lengthScale;
					if (nonlinearity <= 1) { // only count distortion measurements where the map is smooth
						double dxdɸ = (n[0] - s[0])/((ɸMax - ɸMin)/2);
						double dxdλ = (e[0] - w[0])/((λMax - λMin)/2);
						double dydɸ = (n[1] - s[1])/((ɸMax - ɸMin)/2);
						double dydλ = (e[1] - w[1])/((λMax - λMin)/2);
						double dsdɸ = 1;
						double dsdλ = Math.cos((ɸMin + ɸMax)/2);
						
						double s1ps2 = Math.hypot(dxdλ/dsdλ + dydɸ/dsdɸ, dydλ/dsdλ - dxdɸ/dsdɸ);
						double s1ms2 = Math.hypot(dxdλ/dsdλ - dydɸ/dsdɸ, dydλ/dsdλ + dxdɸ/dsdɸ);
						double angularDistortion = (s1ps2 - s1ms2)/(s1ps2 + s1ms2);
						double arealDistortion = (dxdλ*dydɸ - dxdɸ*dydλ)/dsdɸ/dsdλ;
						
						if (angularDistortion < 1/THRESHOLD || angularDistortion > THRESHOLD)
							shearedArea += dA;
						if (arealDistortion < 1/THRESHOLD || arealDistortion > THRESHOLD)
							swollenArea += dA;
						if (unitArea < dA)
							unitArea = dA;
						totalArea += dA;
					}
					else {
						System.out.printf("[%f,%f],[%f,%f],[%f,%f],[%f,%f]\n", n[0], n[1], e[0], e[1], s[0], s[1], w[0], w[1]);
						System.out.println("Nonlinearity detected at "+ɸMin+"N, "+λMin+"E: "+nonlinearity);
					}
				}
			}
			
			System.out.printf("Usable area:%.6f sr\n", totalArea);
			System.out.printf("Anglular:   %8.4f %%\n", 100*shearedArea/totalArea);
			System.out.printf("Areal:      %8.4f %%\n", 100*swollenArea/totalArea);
			System.out.printf("Precision:  %8.4f %%\n", 100*unitArea/totalArea);
			System.out.println();
		}
	}

}
