import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Scanner;

import javax.imageio.ImageIO;

import ellipticFunctions.Jacobi;
import mfc.field.Complex;

/**
 * 
 */

/**
 * @author Justin Kunimune
 *
 */
public class MapProjections {
	private static final int QUINCUNCIAL = 4;
	private static final int EQUIRECTANGULAR = 1;
	private static final int MERCATOR = 2;
	private static final int POLAR = 3;
	
	
	
	
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
		Complex u = new Complex(3.7116*x/width, 3.7116*y/height-1.8558); // don't ask me where 3.7116 come from because I have no idea
		Complex k = new Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus stuff
		Complex ans = Jacobi.cn(u, k);
		double p = 2*Math.atan(ans.abs());
		double theta = Math.atan2(ans.getRe(), ans.getIm());
		double lambda = p-Math.PI/2;
		return getColor(lambda,theta,ref);
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
		Number X = new Number(rt2*Math.tan(p/2), l);
		Number Z = new Number(x,y,true); // the point on the map
//	    		           X*sqrt(-X^2/2 + 1)/2
		Number radical1 = X.times(X.sqrd().over(-2).plus(1)).over(2);
//			       	     asin(X/sqrt(2))/sqrt(2)
		Number arcSin = Number.asin(X.over(rt2)).over(rt2);
//				           X*sqrt(-X^2/4 + 1)
		Number radical2 = X.times(Number.sqrt(X.sqrd().over(-4).plus(1)));
//				           2*Z - pi*sqrt(2)/8
		Number finalBit = Z.times(2).minus(Math.PI*rt2/8);
		//System.out.println("F("+X+","+Z+") = "+(radical1.subtract(arcSin).subtract(radical2).subtract(finalBit)).abs());
		return Number.abs(radical1.minus(radical2).minus(arcSin).minus(finalBit));
	}
}