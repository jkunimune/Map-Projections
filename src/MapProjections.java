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
	private static final String MAP_TYPES = "abcmpsgtw";
	private static final String[] FILE = {"altitude","blackandwhite","color","mars","political","satellite",
			"milkyway","terrain","snowy"};
	private static final String AXES = "cstmjnl";
	private static final String[] AXIS_NAMES = {"Custom","Standard","Transverse","Center of Mass","Jerusalem","Point Nemo","Longest Straight Line"};
	private static final double[] lats = {90,0,29.9792,31.7833,48.8767,-28.5217};
	private static final double[] lons = {0,0,31.1344,35.216,56.6067,141.451};
	private static final double[] thts = {0,0,-35,-35,-45,-18.5};
	
	
	
	
	public static void main(String[] args) throws IOException {
		Scanner in = new Scanner(System.in);
		String response;
		System.out.println("Welcome to the map configurer. You will be asked for a seiries of values. Leaving the field blank at any time will set values to default.");
		
		BufferedImage input;
		int w;
		double x2y;
		double latD, lonD, thtD;
		int projection;
		
		System.out.println("First, enter a file name, or choose a preset map style:");
		for (int i = 0; i < FILE.length; i ++)
			System.out.println(MAP_TYPES.charAt(i)+" --- "+FILE[i]);
		response = in.nextLine();
		if (response.length() == 0)
			response = "s";
		final int index = MAP_TYPES.indexOf(response); // checks for presets
		if (index >= 0) // reads equirectangular from file
			input = ImageIO.read(new File("input/"+FILE[index]+".jpg"));
		else if (response.indexOf(".") >= 0)
			input = ImageIO.read(new File("input/"+response));
		else
			input = ImageIO.read(new File("input/"+response+".jpg"));
		
		System.out.println("And what aspect ratio would you like? (Please enter as a decimal)");
		response = in.nextLine();
		if (response.length() == 0)
			response = "1";
		x2y = Double.parseDouble(response);
		System.out.println("Pixel width?");
		response = in.nextLine();
		if (response.length() == 0)
			response = "540";
		w = Integer.parseInt(response);
		BufferedImage output = new BufferedImage(w,(int)(w/x2y),BufferedImage.TYPE_INT_RGB);
		
		System.out.println("Would you like to use a preset axis, or custom?");
		for (int i = 0; i < AXIS_NAMES.length; i ++)
			System.out.println(AXES.charAt(i)+" --- "+AXIS_NAMES[i]);
		response = in.nextLine();
		if (response.length() == 0)
			response = "s";
		int i = AXES.indexOf(response);
		if (i > 0) { // if it is a preset
			latD = lats[i-1];
			lonD = lons[i-1];
			thtD = thts[i-1];
		}
		else {
			System.out.println("What is the latitude of your desired axis? [-90, 90]");
			response = in.nextLine();
			if (response.length() == 0)
				response = "90";
			latD = Double.parseDouble(response);
			System.out.println("Longitude? [-180, 180]");
			response = in.nextLine();
			if (response.length() == 0)
				response = "0";
			lonD = Double.parseDouble(response);
			System.out.println("What about your orientation? [-180, 180]");
			response = in.nextLine();
			if (response.length() == 0)
				response = "0";
			thtD = Double.parseDouble(response);
		}
		
		System.out.println("Finally, pick a projection:");
		System.out.println(EQUIRECTANGULAR+" --- Equirectangular");
		System.out.println(MERCATOR       +" --- Mercator");
		System.out.println(POLAR    +" --- Polar");
		System.out.println(QUINCUNCIAL    +" --- Peirce Quincuncial");
		
		response = in.nextLine();
		if (response.length() == 0)
			response = Integer.toString(QUINCUNCIAL);
		projection = Integer.parseInt(response);
		System.out.println("Wait...");
		map(input,output,projection,latD,lonD,thtD);
		
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
		return getColor(lat0, lon0, orientation, lambda, theta, ref);
	}
	
	
	public static int equirectangular(final double lat0, final double lon0, final double orientation,
			                          final int width, final int height, int x, int y, BufferedImage ref) { // a basic scale
		return getColor(lat0,lon0,orientation, (double)y/height*Math.PI - Math.PI/2, (double)x/width*2*Math.PI, ref);
	}
	
	
	public static int mercator(final double lat0, final double lon0, final double orientation,
		                       final int width, final int height, int x, int y, BufferedImage ref) { // a popular shape-preserving map
		y -= height/2;
		double phi = Math.atan(Math.sinh((double)y/width*2*Math.PI)); // latitude from 0 to pi
		return getColor(lat0,lon0,orientation, phi, (double)x/width*2*Math.PI, ref);
	}
	
	
	public static int polar(final double lat0, final double lon0, final double orientation,
		                       final int width, final int height, int x, int y, BufferedImage ref) { // the projection used on the UN flag
		double phi = 2*Math.PI*Math.hypot((double)x/width-.5, (double)y/height-.5) - Math.PI/2;
		if (Math.abs(phi) < Math.PI/2)
			return getColor(lat0,lon0,orientation, phi, Math.atan2(y-height/2.0, -x+width/2.0), ref);
		else
			return 0;
	}
	/*END PROJECTION METHODS*/
	
	
	public static int getColor(final double lat0, final double lon0, final double orientation,
			                   double lat1, double lon1, BufferedImage ref) { // returns the color of any coordinate on earth		
		lon1 += orientation;
		double latitude = Math.asin(Math.sin(lat0)*Math.sin(lat1) + Math.cos(lat0)*Math.cos(lon1)*Math.cos(lat1));
		double longitude;
		if (lat0  >= Math.PI/2)
			longitude = lon1+Math.PI;
		else if (lat0 <= -Math.PI/2)
			longitude = -lon1;
		else if (Math.sin(lon1) < 0)
			longitude = lon0 + Math.acos(Math.sin(lat1)/Math.cos(lat0)/Math.cos(latitude)-Math.tan(lat0)*Math.tan(latitude));
		else
			longitude = lon0 - Math.acos(Math.sin(lat1)/Math.cos(lat0)/Math.cos(latitude)-Math.tan(lat0)*Math.tan(latitude));
		
		int x = (int)(longitude*ref.getWidth()/(2*Math.PI));
		int y = (int)((latitude*ref.getHeight()/Math.PI)+ref.getHeight()/2.0);
		
		while (x < 0)
			x += ref.getWidth();
		x %= ref.getWidth();
		if (y < 0)
			y = 0;
		else if (y >= ref.getHeight())
			y = ref.getHeight()-1;
				
		return ref.getRGB(x, y);
	}
	
	
	public static void map(BufferedImage input, BufferedImage output, int projection, double latD, double lonD, double thtD) {
		final int width = output.getWidth();
		final int height = output.getHeight();
		final double lat0 = Math.toRadians(latD);
		final double lon0 = Math.toRadians(lonD);
		final double tht0 = Math.toRadians(thtD+180);
		
		for (int x = 0; x < output.getWidth(); x ++) {
			for (int y = 0; y < output.getHeight(); y ++) {
				switch (projection) {
				case QUINCUNCIAL:
					output.setRGB(x, y, quincuncial(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case EQUIRECTANGULAR:
					output.setRGB(x, y, equirectangular(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case MERCATOR:
					output.setRGB(x, y, mercator(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case POLAR:
					output.setRGB(x, y, polar(lat0,lon0,tht0,width,height,x,y,input));
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