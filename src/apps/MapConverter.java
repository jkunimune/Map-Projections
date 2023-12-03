/**
 * MIT License
 * <p>
 * Copyright (c) 2022 Justin Kunimune
 * <p>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * <p>
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package apps;

import image.PixelMap;
import image.SavableImage;
import maps.Polyhedral;
import maps.Projection;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Stream;

import static java.lang.Math.sqrt;

/**
 * a script that automatically converts a bunch of inputs into a particular projection.
 */
public class MapConverter {
	public static void main(String[] args) throws IOException {

		// the directory with the equirectangular maps to convert
		String directory = "C:\\Users\\justi\\Downloads\\new_dymaxion_maps";
		// the desired map projection
		Projection projection = Polyhedral.DYMAXION;

		Stream<Path> paths = Files.walk(Paths.get(directory));
		Iterable<Path> pathIterable = paths::iterator;

		// iterate thru the directory
		for (Path inputPath: pathIterable) {
			// look for images that are not dymaxion projections
			if (inputPath.toString().endsWith(".jpg") ||
			    inputPath.toString().endsWith(".jpeg") ||
			    inputPath.toString().endsWith(".tif") ||
			    inputPath.toString().endsWith(".tiff") ||
			    inputPath.toString().endsWith(".png") &&
			    !inputPath.toString().endsWith(".dymaxion.png")) {
				System.out.println(inputPath);
				PixelMap inputImage = new PixelMap(inputPath.toFile());

				// reduce the area by 2 to avoid pixelation
				double area = inputImage.getWidth() * inputImage.getHeight() / 2.;
				// and calculate the new dimensions according to Dymaxion's aspect ratio
				int width = (int) sqrt(area * projection.getAspectRatio());
				int height = (int) sqrt(area / projection.getAspectRatio());

				// generate the new map
				BufferedImage outputImage = MapDesignerRaster.calculate(
					  width, height, 2, inputImage, projection,
					  null, true, 0,
					  null, null, null);

				// update the filename and save to disk
				String outputPath = inputPath.toString();
				outputPath = outputPath.replace(".jpg", ".png");
				outputPath = outputPath.replace(".jpeg", ".png");
				outputPath = outputPath.replace(".tif", ".png");
				outputPath = outputPath.replace(".tiff", ".png");
				outputPath = outputPath.replace(".png", ".dymaxion.png");
				SavableImage.savable(outputImage).save(new File(outputPath));
			}
		}

	}
}
