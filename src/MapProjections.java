import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

/**
 * 
 */

/**
 * @author Justin Kunimune
 *
 */
public class MapProjections {
	private static final int QUINCUNCIAL = 0;
	private static final int MERCATOR = 1;
	
	
	
	
	public static void main(String[] args) throws IOException {
		BufferedImage input = ImageIO.read(new File("realistic.jpg"));
		
		BufferedImage output = new BufferedImage(540,540,BufferedImage.TYPE_INT_RGB);
		
		map(input,output,MERCATOR);
		
		saveImage(output);
	}
	
	
	
	
	/* PROJECTION METHODS: Return RGB at a given pixel based on a reference map and a unique projection method */
	public static int quincuncial(final double lat0, final double lon0, final double orientation,
			                      double x, double y, BufferedImage ref) { // a tessalatable square map
		return 0;
	}
	
	
	/* PROJECTION METHODS: Return RGB at a given pixel based on a reference map and a unique projection method */
	public static int mercator(final double lat0, final double lon0, final double orientation,
			                      double x, double y, BufferedImage ref) { // a popular shape-preserving map
		return 0;
	}
	
	
	public static void map(BufferedImage input, BufferedImage output, int projection) {
		for (int x = 0; x < output.getWidth(); x ++) {
			for (int y = 0; y < output.getHeight(); y ++) {
				switch (projection) {
				case QUINCUNCIAL:
					output.setRGB(x, y, quincuncial(0,0,0,x,y,input));
					break;
				case MERCATOR:
					output.setRGB(x, y, mercator(0,0,0,x,y,input));
					break;
				}
			}
		}
	}
	
	
	private static void saveImage(BufferedImage img) {
		try {
			File outputFile = new File("output.jpg");
			ImageIO.write(img, "jpg", outputFile);
		} catch (IOException e) {}
	}
	
	
	public static final byte[] valuesFromRGB(int RGB) { // converts an int to three bytes
		byte[] output = new byte[3];
		output[0] = (byte)((RGB & 0xFF0000) >> 16); // red
		output[1] = (byte)((RGB & 0x00FF00) >> 8); // green
		output[2] = (byte)(RGB & 0x0000FF); // blue
		return output;
	}
	
	
	public static final int RGBFromValues(byte r, byte g, byte b) { // converts three bytes to an int
		return (r<<16) + (g<<8) + b;
	}
}
