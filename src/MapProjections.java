import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Scanner;

import javax.imageio.ImageIO;

/**
 * 
 */

/**
 * @author Justin Kunimune
 *
 */
public class MapProjections {
	private static final int QUINCUNCIAL = 1;
	private static final int EQUIRECTANGULAR = 0;
	private static final int MERCATOR = 2;
	
	
	
	
	public static void main(String[] args) throws IOException {
		Scanner in = new Scanner(System.in);
		String response;
		
		System.out.println("Realistic or Simplified? (r or s)");
		response = in.nextLine();
		BufferedImage input;
		if (response.toLowerCase().charAt(0) == 'r') // reads equirectangular from file
			input = ImageIO.read(new File("input/realistic.jpg"));
		else
			input = ImageIO.read(new File("input/simplified.jpg"));
		
		System.out.println("And what aspect ratio would you like? (Please enter a decimal)");
		response = in.nextLine();
		double x2y = Double.parseDouble(response);
		System.out.println("Pixel width?");
		response = in.nextLine();
		int w = Integer.parseInt(response);
		BufferedImage output = new BufferedImage(w,(int)(w/x2y),BufferedImage.TYPE_INT_RGB);
		
		System.out.println("Finally, pick a projection:\n");
		System.out.println(EQUIRECTANGULAR+" --- Equirectangular");
		System.out.println(QUINCUNCIAL    +" --- Peirce Quincuncial");
		System.out.println(MERCATOR       +" --- Mercator");
		response = in.nextLine();
		int projection = Integer.parseInt(response);
		map(input,output,projection);
		
		saveImage(output);
		
		in.close();
		
		System.out.println("Done!");
	}
	
	
	
	
	/* PROJECTION METHODS: Return RGB at a given pixel based on a reference map and a unique projection method */
	public static int quincuncial(final double lat0, final double lon0, final double orientation,
			                      final int width, final int height, int x, int y, BufferedImage ref) { // a tessalatable square map
		return 0;
	}
	
	
	public static int equirectangular(final double lat0, final double lon0, final double orientation,
			                          final int width, final int height, int x, int y, BufferedImage ref) { // a popular shape-preserving map
		
		return ref.getRGB((int)((double)ref.getWidth()*x/width), (int)((double)ref.getHeight()*y/height));
	}
	
	
	public static int mercator(final double lat0, final double lon0, final double orientation,
		                       final int width, final int height, int x, int y, BufferedImage ref) { // a popular shape-preserving map
		y -= height/2;
		double phi = Math.atan(Math.sinh((double)y/width*2*Math.PI)) + Math.PI/2; // latitude from 0 to pi
		return ref.getRGB((int)((double)x/width*ref.getWidth()), (int)(phi/Math.PI*ref.getHeight()));
	}
	
	
	public static void map(BufferedImage input, BufferedImage output, int projection) {
		final int width = output.getWidth();
		final int height = output.getHeight();
		
		for (int x = 0; x < output.getWidth(); x ++) {
			for (int y = 0; y < output.getHeight(); y ++) {
				switch (projection) {
				case QUINCUNCIAL:
					output.setRGB(x, y, quincuncial(0,0,0,width,height,x,y,input));
					break;
				case EQUIRECTANGULAR:
					output.setRGB(x, y, equirectangular(0,0,0,width,height,x,y,input));
					break;
				case MERCATOR:
					output.setRGB(x, y, mercator(0,0,0,width,height,x,y,input));
					break;
				default:
					System.err.println("Justin, you forgot to add a projection to the switch case!");
				}
			}
		}
	}
	
	
	private static void saveImage(BufferedImage img) {
		try {
			File outputFile = new File("output/myMap.jpg");
			ImageIO.write(img, "jpg", outputFile);
		} catch (IOException e) {}
	}
}
