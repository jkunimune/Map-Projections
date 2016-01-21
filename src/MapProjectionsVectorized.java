import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

import org.apache.commons.math3.complex.Complex;

/**
 * 
 */

/**
 * @author Justin Kunimune
 *
 */
public class MapProjectionsVectorized {
	private static final int QUINCUNCIAL = 11;
	private static final int EQUIRECTANGULAR = 1;
	private static final int MERCATOR = 2;
	private static final int POLAR = 5;
	private static final int GALL = 3;
	private static final int SINUSOIDAL = 12;
	private static final int STEREOGRAPHIC = 6;
	private static final int ORTHOGONAL = 8;
	private static final int LEMONS = 13;
	private static final int EA_AZIMUTH = 7;
	private static final int EA_CYLINDER = 4;
	private static final int CONICAL = 10;
	private static final int GNOMIC = 9;
	private static final int QUINSHIFT = 14;
	
	private static final String AXES = "cstmjnlx123";
	private static final String[] AXIS_NAMES = {"Custom","Standard","Transverse","Center of Mass","Jerusalem",
			"Point Nemo","Longest Line","Longest Line Transverse","Cylindrical","Conical","Quincuncial"};
	private static final double[] lats = {90,0,29.9792,31.7833,48.8767,-28.5217,-46.4883,-35,10,59};
	private static final double[] lons = {0,0,31.1344,35.216,56.6067,141.451,16.5305,-13.6064,-115,19};
	private static final double[] thts = {0,180,-32,-35,-45,107,137,145,150,50};
	
	
	
	
	public static void main(String[] args) {
		/*Complex one = new Complex(1);
		Complex rt2 = new Complex(Math.sqrt(.5));
		Complex two = new Complex(2);
		Complex piover2 = new Complex(Math.PI/2);
		Complex piover4 = new Complex(Math.PI/4);
		
		System.out.println("F(pi/4,rt2/2) = "+F(piover4,rt2));
		System.out.println("F(1,rt2/2) = "+F(one,rt2));
		System.out.println("F(pi/2,rt2/2) = "+F(piover2,rt2));
		System.out.println("F(2,rt2/2) = "+F(two,rt2));
		System.out.println("F(rt2/2,rt2/2) = "+F(rt2,rt2));*/
		
		Scanner in = new Scanner(System.in);
		String response;
		System.out.println("Welcome to the map configurer. You will be asked for a series of values. Leaving the field blank at any time will set values to default.");
		
		BufferedReader input = null;
		BufferedWriter output = null;
		double latD, lonD, thtD;
		int projection;
		
		try {
			input = new BufferedReader(new FileReader("input/vectorized.svg"));
			output = new BufferedWriter(new FileWriter("output/vectorized.svg"));
		}
		catch (FileNotFoundException e) {
			System.err.println("Where is the freaking file?");
		}
		catch (IOException e) {
			System.err.println("IOException?! What does that even mean?");
		}
				
		while (true) {
			try {
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
				
				break;
			} catch (Error e) {
				System.out.println("I don't like that. Enter something else.");
			}
		}
		
		while (true) {
			try {
				System.out.println("Finally, pick a projection:");
				System.out.println(EQUIRECTANGULAR+" --- Equirectangular");
				System.out.println(MERCATOR       +" --- Mercator");
				System.out.println(GALL           +" --- Gall Stereographic");
				System.out.println(EA_CYLINDER    +" --- Cylindrical Equal-Area");
				System.out.println(POLAR          +" --- Polar");
				System.out.println(STEREOGRAPHIC  +" --- Stereographic");
				System.out.println(EA_AZIMUTH     +" --- Azimuthal Equal-Area");
				System.out.println(ORTHOGONAL     +" --- Orthogonal");
				System.out.println(GNOMIC         +" --- Gnomic");
				System.out.println(CONICAL        +" --- Lambert Conic");
				System.out.println(QUINCUNCIAL    +" --- Peirce Quincuncial");
				System.out.println(SINUSOIDAL     +" --- Sinusoidal");
				System.out.println(LEMONS         +" --- BURN LIFE'S HOUSE DOWN");
				System.out.println(QUINSHIFT      +" --- Shifted Quincuncial");
				
				response = in.nextLine();
				if (response.length() == 0)
					response = Integer.toString(QUINCUNCIAL);
				projection = Integer.parseInt(response);
				
				break;
			} catch (Error e) {
				System.out.println("I don't like that response. Enter something else.");
			}
		}
		System.out.println("Wait...");
		parse(input,output,projection,latD,lonD,thtD);
		
		try {
			output.close();
		} catch (IOException e) {
			System.err.println("Seriously?");
		}
		in.close();
		
		System.out.println("Done!");
	}
	
	
	
	
	/* PROJECTION METHODS: Return map coordinates based on a given latitude and longitude and a projection algorithm */
	public static double[] equirectangular(final double lat0, final double lon0, final double orientation,
                                      double lat, double lon) { // a basic scale
		double[] output = new double[2];
		output[0] = lon*2000/(2*Math.PI);
		output[1] = lat*1000/Math.PI;
		return output;
	}
	
	
	public static double[] polar(final double lat0, final double lon0, final double orientation,
            double lat, double lon) { // UN flag
		double[] output = new double[2];
		output[0] = (lat+Math.PI/2)*1000/Math.PI * Math.cos(lon);
		output[1] = (lat+Math.PI/2)*1000/Math.PI * Math.sin(lon);
		return output;
	}
	
	
	public static double[] quincuncial(final double lat0, final double lon0, final double orientation,
            double lat, double lon) { // awesomeness
		double[] output = new double[2];
				
		final double wMag = Math.tan(lat/2+Math.PI/4);
		final Complex w = new Complex(wMag*Math.cos(lon), wMag*Math.sin(lon));
		final Complex k = new Complex(Math.sqrt(0.5));
		Complex z = F(w.acos(),k);
		if (z.isInfinite() || z.isNaN() || z.abs() > 10)
			z = new Complex(0);
				
		output[0] = -z.getReal()*1000;
		output[1] = z.getImaginary()*1000;
		return output;
		}
		
		
	public static double[] shiftquin(final double lat0, final double lon0, final double orientation,
            double lat, double lon) { // more awesomeness
		double[] output = new double[2];
				
		if (lat >= 0) {
			final double wMag = Math.tan(-lat/2+Math.PI/4);
			final Complex w = new Complex(wMag*Math.cos(lon), wMag*Math.sin(lon));
			final Complex k = new Complex(Math.sqrt(0.5));
			Complex z = F(w.acos(),k);
					
			output[0] = -z.getImaginary()*1000;
			output[1] = z.getReal()*1000;
			return output;
		}
		else {
			final double wMag = Math.tan(lat/2+Math.PI/4);
			final Complex w = new Complex(wMag*Math.cos(lon), wMag*Math.sin(lon));
			final Complex k = new Complex(Math.sqrt(0.5));
			Complex z = F(w.acos(),k);
					
			output[0] = -z.getReal()*1000;
			output[1] = z.getImaginary()*1000;
			return output;
		}
	}
	/*END PROJECTION METHODS*/
	
	
	public static void parse(BufferedReader input, BufferedWriter output, int projection, double latD, double lonD, double thtD) {
		final double lat0 = Math.toRadians(latD);
		final double lon0 = Math.toRadians(lonD);
		final double tht0 = Math.toRadians(thtD+180);
		
		try {
			
			int character;
			do { // skips to the juicy bits
				character = input.read(); // reads some boring stuff
				output.write(character); // writes some boring stuff
			} while (character != 'M'); // stops after writing M
			output.write(input.read()); // reads and writes a space
			
			String number;
			character = input.read();
			do {
				number = "";
				do { // reads the longitude
					number = number + (char)character; // saves the number character
					character = input.read(); // reads the next character
				} while (character != ',');
				double xCoord = Double.parseDouble(number);
				number = "";
				character = input.read(); // reads a hyphen
				do { // gets the latitude in the same way
					number = number + (char)character; // reads the latitude
					character = input.read();
				} while (character != ' ');
				double yCoord = Double.parseDouble(number);
				
				double[] latLon = convertCoords(xCoord, yCoord, lat0,lon0,tht0);
				double[] mapCoords = getNewCoords(latLon[0], latLon[1], projection, lat0,lon0,tht0);
				output.write(mapCoords[0] + "," + mapCoords[1]); // writes the new coordinates in
				
				do {
					output.write(character); // writes
					character = input.read(); // reads
				} while (character == ' ' || character == 'z' || character == 'M' || character == 'L'); // until it hits a hyphen or a quotation mark
			} while (character != '\"'); // the quotation mark means it is done.
			
			do {
				output.write(character);
				character = input.read(); // finishes up
			} while (character != -1);
			
		} catch (IOException e) {
			System.err.println("You suck.");
		}
	}
	
	
	public static double[] getNewCoords(double lat, double lon, int projection, double lat0, double lon0, double tht0) {
		switch (projection) {
			case EQUIRECTANGULAR:
				return equirectangular(lat0,lon0,tht0,lat,lon);
			case POLAR:
				return polar(lat0,lon0,tht0,lat,lon);
			case QUINCUNCIAL:
				return quincuncial(lat0,lon0,tht0,lat,lon);
			case QUINSHIFT:
				return shiftquin(lat0,lon0,tht0,lat,lon);
			default:
				System.err.println("Justin, you forgot to add a projection to the switch case! (or you forgot a break;)");
				return null;
		}
	}

	
	public static double[] convertCoords(double x, double y, double lat0,double lon0,double tht0) { // converts svg coordinates to lat and lon
		final double width = 4378.125;
		final double height = 2534.9375;//2619.25547895;//2434.9375;
		double latF = y*Math.PI/height + 2.533457922281304 + 0.045 - Math.PI; // final latitude from y
		double lonF = x*2*Math.PI/width + 2.5849410045340258 + 0.19198621771937625346160598453375; // final longitude from x
		double[] latLon = new double[2];
			
		Vector r0 = new Vector (1, lat0, lon0);
		Vector rF = new Vector (1, latF, lonF);
		Vector r0XrF = r0.cross(rF);
		Vector r0Xk = r0.cross(Vector.K);
		
		latLon[0] = Math.asin(r0.dot(rF)); // relative latitude
		if (lat0 == Math.PI/2 || lat0 == -Math.PI/2) // accounts for all the 0/0 errors at the poles
			latLon[1] = lonF;
		else
			latLon[1] = Math.acos(r0XrF.dot(r0Xk)/(r0XrF.abs()*r0Xk.abs())); // relative longitude
		 
		if (Double.isNaN(latLon[1]))
			latLon[1] = 0;
		else if (r0XrF.cross(r0Xk).dot(r0)/(r0XrF.abs()*r0Xk.abs()) > 0) // it's a plus-or-minus arccos.
			latLon[1] = 2*Math.PI-latLon[1];
		latLon[1] = (latLon[1]+tht0+Math.PI) % (2*Math.PI);
		
		return latLon;
	}
	
	
	public static final Complex F(Complex phi, final Complex k) { // series solution to incomplete elliptic integral of the first kind
		Complex sum = new Complex(0);
	    Complex i_n = phi;
	    
	    for (int n = 0; n < 100; n ++) {
	        if (n > 0)
	            i_n = i_n.multiply((2.0*n-1)/(2.0*n)).subtract(phi.cos().multiply(phi.sin().pow(2.0*n-1)).divide(2.0*n));
	        sum = sum.add(i_n.multiply(Math.abs(combine(-.5,n))).multiply(k.pow(2.0*n)));
	    }
	    
	    return sum;
	}
	
	
//	public static final Complex F(Complex phi, final Complex k) { // riemann solution to incomplete elliptic integral of the first kind
//		Complex theta = new Complex(0);
//		Complex dTheta = phi.divide(500.0);
//		Complex area = new Complex(0);
//		final Complex one = new Complex(1);
//		
//		while (theta.abs() < phi.abs()) {
//			theta = theta.add(dTheta);
//			area = area.add(one.subtract(k.pow(2).multiply(theta.sin().pow(2))).sqrt().multiply(dTheta));
//		}
//		
//		return area;
//	}
	
	
	public static final double combine(double n, int k) {
		double output = 1;
		for (int i = k; i > 0; i --) {
			output *= (n+i-k)/i;
		}
		return output;
	}
	
	
	public static final Complex abs(Complex z) {
		return new Complex(z.abs());
	}
}