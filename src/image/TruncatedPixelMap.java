/*
 * MIT License
 *
 * Copyright (c) 2022 Justin Kunimune
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

import java.io.File;
import java.io.IOException;

import static java.lang.Math.PI;

/**
 * a PixelMap that only loads one octant of the globe, to save memory
 */
public class TruncatedPixelMap extends PixelMap {
	public TruncatedPixelMap(File f) throws IOException {
		super(f);
	}

	public int getWidth() {
		return super.getWidth()*4;
	}


	public int getHeight() {
		return super.getHeight()*2;
	}


	public int getArgb(double lat, double lon) {
		if (lat < 0 || lon < 0 || lon > PI/2)
			return 0;
		else
			return super.getArgb(2*lat - PI/2, 4*lon - PI);
	}
}
