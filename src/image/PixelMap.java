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

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

/**
 * An input equirectangular map based on a raster image file
 * 
 * @author jkunimune
 */
public class PixelMap {
	
	private BufferedImage pixels;
	
	
	public PixelMap(File f) throws IOException {
		pixels = ImageIO.read(f);
	}
	
	
	public int getArgb(double lat, double lon) {
		double x = 0.5 + lon/(2*Math.PI);
		x = (x - Math.floor(x)) * pixels.getWidth();
		
		double y = pixels.getHeight()*(.5 - lat/Math.PI);
		if (y < 0)
			y = 0;
		else if (y >= pixels.getHeight())
			y = pixels.getHeight() - 1;
		
		return (0xFF000000) | pixels.getRGB((int) x, (int) y);
	}
}
