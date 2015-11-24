import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Scanner;

import javax.imageio.ImageIO;

import org.apache.commons.math3.complex.Complex;

/**
 * 
 */

/**
 * @author Justin Kunimune
 *
 */
public class MapProjections {
	private static final int QUINCUNCIAL = 3;
	private static final int EQUIRECTANGULAR = 1;
	private static final int MERCATOR = 2;
	private static final int POLAR = 4;
	
	
	
	
	public static void main(String[] args) throws IOException {
		Scanner in = new Scanner(System.in);
		String response;
		
		System.out.println("Realistic or Simplified? (r or s)");
		response = in.nextLine();
		if (response.length() == 0)
			response = "c";
		BufferedImage input;
		if (response.toLowerCase().charAt(0) == 'r') // reads equirectangular from file
			input = ImageIO.read(new File("input/realistic.jpg"));
		else if (response.toLowerCase().charAt(0) == 's')
			input = ImageIO.read(new File("input/simplified.jpg"));
		else
			input = ImageIO.read(new File("input/color.png"));
		
		System.out.println("And what aspect ratio would you like? (Please enter as a decimal)");
		response = in.nextLine();
		if (response.length() == 0)
			response = "1";
		double x2y = Double.parseDouble(response);
		System.out.println("Pixel width?");
		response = in.nextLine();
		if (response.length() == 0)
			response = "540";
		int w = Integer.parseInt(response);
		BufferedImage output = new BufferedImage(w,(int)(w/x2y),BufferedImage.TYPE_INT_RGB);
		
		System.out.println("Finally, pick a projection:\n");
		System.out.println(EQUIRECTANGULAR+" --- Equirectangular");
		System.out.println(MERCATOR       +" --- Mercator");
		System.out.println(POLAR    +" --- Polar");
		System.out.println(QUINCUNCIAL    +" --- Peirce Quincuncial");
		
		response = in.nextLine();
		if (response.length() == 0)
			response = "3";
		int projection = Integer.parseInt(response);
		System.out.println("Wait...");
		map(input,output,projection);
		
		saveImage(output);
		
		in.close();
		
		System.out.println("Done!");
	}
	
	
	
	
	/* PROJECTION METHODS: Return RGB at a given pixel based on a reference map and a unique projection method */
	public static int quincuncial(final double lat0, final double lon0, final double orientation,
			                      final int width, final int height, int x, int y, BufferedImage ref) { // a tessalatable square map
		final double dx = .01;
		
		double phi, lambda; // the correct phi and lambda combination is the one where func goes to 0.
		double phiMin = -Math.PI/2;
		double phiMax = Math.PI/2;
		double lamMin = -Math.PI;
		double lamMax = Math.PI;
		
		double[] error = new double[3];
		int i = 0;
		do {
			phi = (phiMax+phiMin)/2.0;
			lambda = (lamMax+lamMin)/2.0;
			
			error[0] = func(phi, lambda, 2.0*x/width-1, 2.0*y/height-1);
			error[1] = func(phi+dx, lambda, 2.0*x/width-1, 2.0*y/height-1);
			error[2] = func(phi, lambda+dx, 2.0*x/width-1, 2.0*y/height-1);
			
			if (error[0] > error[1]) // phi is too small
				phiMin = phi;
			else // phi is too big
				phiMax = phi;
			if (error[0] > error[2]) // lambda is too small
				lamMin = lambda;
			else // lambda is too big
				lamMax = lambda; 

			i ++;
		} while (phiMax-phiMin > Math.PI/ref.getHeight() && lamMax-lamMin > 2*Math.PI/ref.getWidth() && i < 256);
		
		return getColor(phi,lambda,ref);
	}
	
	
	public static int equirectangular(final double lat0, final double lon0, final double orientation,
			                          final int width, final int height, int x, int y, BufferedImage ref) { // a popular shape-preserving map
		
		return getColor((double)y/height*Math.PI - Math.PI/2, (double)x/width*2*Math.PI, ref);
	}
	
	
	public static int mercator(final double lat0, final double lon0, final double orientation,
		                       final int width, final int height, int x, int y, BufferedImage ref) { // a popular shape-preserving map
		y -= height/2;
		double phi = Math.atan(Math.sinh((double)y/width*2*Math.PI)); // latitude from 0 to pi
		return getColor(phi, (double)x/width*2*Math.PI, ref);
	}
	
	
	public static int polar(final double lat0, final double lon0, final double orientation,
		                       final int width, final int height, int x, int y, BufferedImage ref) { // the projection used on the UN flag
		double phi = 2*Math.PI*Math.hypot((double)x/width-.5, (double)y/height-.5) - Math.PI/2;
		if (Math.abs(phi) < Math.PI/2)
			return getColor(phi, Math.atan2(y-height/2.0, x-width/2.0), ref);
		else
			return 0;
	}
	
	
	public static int getColor(double lattitude, double longitude, BufferedImage ref) { // returns the color of any coordinate on earth
		if (lattitude >= Math.PI/2)
			return ref.getRGB(0, 0);
		return ref.getRGB((int)(longitude*ref.getWidth()/(2*Math.PI) +ref.getWidth())%ref.getWidth(), (int)((lattitude*ref.getHeight()/Math.PI)+ref.getHeight()/2.0));
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
				case POLAR:
					output.setRGB(x, y, polar(0,0,0,width,height,x,y,input));
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
	
	
	/* PRECONDITION: -1 <= x,y <= 1 */
	private static double func(double p, double l, double x, double y) { // a super complicated function that calculates stuff for peirce quincuncial
		final double rt2 = Math.sqrt(2);
		Complex X = new Complex(rt2*Math.tan(p/2)*Math.cos(l), rt2*Math.tan(p/2)*Math.sin(l));
		Complex Z = new Complex(x,y); // the point on the map
//	    		           X*sqrt(-X^2/2 + 1)/2
		Complex radical1 = X.multiply(X.multiply(X).divide(-2).add(1).pow(.5)).divide(2);
//			       	     asin(X/sqrt(2))/sqrt(2)
		Complex arcSin = X.divide(rt2).asin().divide(rt2);
//				           X*sqrt(-X^2/4 + 1)
		Complex radical2 = X.multiply(X.multiply(X).divide(-4).add(1).pow(.5));
//				           2*Z + pi*sqrt(2)/8
		Complex finalBit = Z.multiply(2).subtract(Math.PI*rt2/8);
		System.out.println("F("+X+","+Z+") = "+(radical1.subtract(arcSin).subtract(radical2).subtract(finalBit)).abs());
		return (radical1.subtract(arcSin).subtract(radical2).subtract(finalBit)).abs();
	}
}