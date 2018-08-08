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

import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import javafx.embed.swing.SwingFXUtils;
import javafx.scene.image.Image;

/**
 * An image that can be saved. It provides a convenient interface for both vector and raster images.
 * 
 * @author jkunimune
 */
public interface SavableImage {
	
	public abstract void save(File file) throws IOException;
	
	public static SavableImage savable(RenderedImage img) {
		return new SavableImage() {
			public void save(File file) throws IOException {
				String filename = file.getName();
				String extension = filename.contains(".") ?
						filename.substring(filename.lastIndexOf(".")+1) : "";
				ImageIO.write(img, extension, file);
			}
		};
	}
	
	public static SavableImage savable(Image img) {
		return new SavableImage() {
			public void save(File file) throws IOException {
				String filename = file.getName();
				String extension = filename.contains(".") ?
						filename.substring(filename.lastIndexOf(".")+1) : "";
				ImageIO.write(SwingFXUtils.fromFXImage(img, null), extension, file);
			}
		};
	}
	
}
