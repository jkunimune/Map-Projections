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

import java.awt.Transparency;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import static java.lang.Math.PI;
import static java.lang.Math.floor;

/**
 * An input equirectangular map based on a raster image file
 * 
 * @author jkunimune
 */
public class PixelMap {
	
	private final BufferedImage pixels;
	private final WritableRaster alphaPixels;
	
	
	public PixelMap(File f) throws IOException {
		pixels = ImageIO.read(f);
		alphaPixels = pixels.getAlphaRaster();
	}


	public int getWidth() {
		return this.pixels.getWidth();
	}


	public int getHeight() {
		return this.pixels.getHeight();
	}
	
	
	public int getArgb(double lat, double lon) {
		double x = 0.5 + lon/(2*PI);
		x = (x - floor(x)) * pixels.getWidth();
		
		double y = pixels.getHeight()*(.5 - lat/PI);
		if (y < 0)
			y = 0;
		else if (y >= pixels.getHeight())
			y = pixels.getHeight() - 1;

		int alpha;
		if (pixels.getTransparency() != Transparency.OPAQUE)
			alpha = alphaPixels.getPixel((int) x, (int) y, new int[1])[0];
		else
			alpha = 0xFF;
		int rgb = pixels.getRGB((int) x, (int) y);
		
		return (alpha << 24) | rgb;
	}
}
