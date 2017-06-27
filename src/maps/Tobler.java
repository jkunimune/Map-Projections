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
package maps;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;

/**
 * A class of values and functions used to approximate the Tobler projection
 * 
 * @author jkunimune
 */
public class Tobler {

	private static final double[][] TABLE = {
			{  -90,-85,-80,-75,-70,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-05, 00,
					05, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90 },
			{-5.000000,-4.921966,-4.782840,-4.610387,-4.410173,-4.184957,-3.940243,-3.678861,-3.398039,-3.103326,-2.794726,-2.475016,-2.144195,-1.802264,-1.454780,-1.096182,-0.734807,-0.367979, 0.000000,
						0.367979, 0.734807, 1.096182, 1.454780, 1.802264, 2.144195, 2.475016, 2.794726, 3.103326, 3.398039, 3.678861, 3.940243, 4.184957, 4.410173, 4.610387, 4.782840, 4.921966, 5.000000 },
			{ 0.000000, 0.022944, 0.030411, 0.035519, 0.039413, 0.042524, 0.045070, 0.047181, 0.048940, 0.050406, 0.051625, 0.052627, 0.053442, 0.054091, 0.054595, 0.054973, 0.055242, 0.055429, 0.055555,
						0.055429, 0.055242, 0.054973, 0.054595, 0.054091, 0.053442, 0.052627, 0.051625, 0.050406, 0.048940, 0.047181, 0.045070, 0.042524, 0.039413, 0.035519, 0.030411, 0.022944, 0.000000 }
		}; //these values were calculated with some kind of iterative approach in 1973.
	
	private static final SplineInterpolator SI = new SplineInterpolator(); //I don't understand why these can't be static functions, but fine
	
	
	
	public static final double xfacFromLat(double lat) {
		return SI.interpolate(TABLE[0], TABLE[2]).value(Math.toDegrees(lat)); //I also don't understand why I can't call interpolate into a final Object instead of calling it every time. It's weird. If I try to do that, then invoking this class terminates the thread, like there's something so awful about calling that method outside of any method that Java just dies
	}
	
	
	public static final double yFromLat(double lat) {
		return SI.interpolate(TABLE[0], TABLE[1]).value(Math.toDegrees(lat));
	}
	
	
	public static final double latFromY(double pdfe) {
		return Math.toRadians(SI.interpolate(TABLE[1], TABLE[0]).value(pdfe));
	}
	
	
	public static final double xfacFromY(double pdfe) {
		return SI.interpolate(TABLE[1], TABLE[2]).value(pdfe);
	}

}
