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
package image;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Path2D;
import java.util.List;

import static java.lang.Math.pow;

/**
 * A collection of (currently just one) methods for dealing with images and coordinates
 * 
 * @author jkunimune
 */
public class ImageUtils {
	
	public static int blend(int[] colors) {
		return blend(colors, 2.2);
	}
	
	
	public static int blend(int[] colors, double gamma) {
		int a_tot = 0;
		int r_tot = 0;
		int g_tot = 0;
		int b_tot = 0;
		for (int argb: colors) {
			int a = ((argb >> 24)&0xFF);
			a_tot += a;
			r_tot += a*(pow((argb>>16)&0xFF, gamma));
			g_tot += a*(pow((argb>> 8)&0xFF, gamma));
			b_tot += a*(pow((argb>> 0)&0xFF, gamma));
		}
		if (a_tot == 0)	return 0;
		else
			return (a_tot/colors.length << 24) +
					((int)pow((double)r_tot/a_tot, 1/gamma) << 16) +
					((int)pow((double)g_tot/a_tot, 1/gamma) << 8) +
					((int)pow((double)b_tot/a_tot, 1/gamma) << 0);
	}
	
	
	public static void drawSVGPath(List<Path.Command> path, Color stroke, float strokeWidth, boolean antialias, Graphics2D g) {
		g.setStroke(new BasicStroke(strokeWidth, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND));
		g.setColor(stroke);
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				antialias ? RenderingHints.VALUE_ANTIALIAS_ON : RenderingHints.VALUE_ANTIALIAS_OFF);
		drawSVGPath(path, g);
	}

	public static void drawSVGPath(List<Path.Command> path, Graphics2D g) {
		Path2D awtPath = new Path2D.Double(Path2D.WIND_NON_ZERO, path.size());
		for (Path.Command svgCmd: path) {
			switch (svgCmd.type) {
			case 'M':
				awtPath.moveTo(svgCmd.args[0], svgCmd.args[1]);
				break;
			case 'L':
				awtPath.lineTo(svgCmd.args[0], svgCmd.args[1]);
				break;
			case 'Z':
				awtPath.closePath();
			}
		}
		g.draw(awtPath);
	}
}
