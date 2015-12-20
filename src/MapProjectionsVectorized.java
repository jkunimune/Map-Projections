import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

import ellipticFunctions.Jacobi;
import mfc.field.Complex;

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
	private static final double[] thts = {0,180,-32,-35,-45,71.5,137,145,150,50};
	
	
	
	
	public static void main(String[] args) {
		Scanner in = new Scanner(System.in);
		String response;
		System.out.println("Welcome to the map configurer. You will be asked for a seiries of values. Leaving the field blank at any time will set values to default.");
		
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
		output[0] = (Math.PI/2-lat)*1000/Math.PI * Math.cos(lon);
		output[1] = (Math.PI/2-lat)*1000/Math.PI * Math.sin(lon);
		return output;
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
			default:
				System.err.println("Justin, you forgot to add a projection to the switch case! (or you forgot a break;)");
				return null;
		}
	}

	
	public static double[] convertCoords(double x, double y, double lat0,double lon0,double tht0) { // converts svg coordinates to lat and lon
		final double width = 4378.125;
		final double height = 2619.25547895;//2434.9375;
		double latF = y*Math.PI/height + 2.533457922281304 - Math.PI; // final latitude from y
		double lonF = x*2*Math.PI/width + 2.5849410045340258; // final longitude from x
		Vector r0 = new Vector (1, lat0, lon0);
		Vector rF = new Vector (1, latF, lonF);
		Vector r0XrF = r0.cross(rF);
		Vector r0Xk = r0.cross(Vector.K);
		
		double[] latLon = new double[2];
		latLon[0] = Math.asin(r0.dot(rF)); // relative latitude
		latLon[1] = Math.acos(r0XrF.dot(r0Xk)/(r0XrF.abs()*r0Xk.abs())); // relative longitude
		
		if (Double.isNaN(latLon[1]))
			latLon[1] = 0;
		else if (r0XrF.cross(r0Xk).dot(r0)/(r0XrF.abs()*r0Xk.abs()) > 0) // it's a plus-or-minus arccos.
			latLon[1] = 2*Math.PI-latLon[1];
		latLon[1] = (latLon[1]+tht0) % (2*Math.PI);
		
		return latLon;
	}
}