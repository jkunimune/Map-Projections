import java.awt.Desktop;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

import ellipticFunctions.Jacobi;
import mfc.field.Complex;

/**
 * 
 */

/**
 * @author Justin Kunimune
 *
 */
public class MapProjections implements ActionListener {
	private static final String[] PROJ = {"Equirectangular","Mercator","Gall Stereographic",
			"Cylindrical Equal-Area","Polar","Stereographic","Azimuthal Equal-Area","Orthogonal","Gnomic",
			"Lambert Conical","Winkel Tripel","Van der Grinten","Mollweide","Hammer","Sinusoidal","Lemons",
			"Pierce Quincuncial","Magnifier","Guyou Hemisphere-in-a-Square","Rectus Aequilibrium" };
	private static final int[] DEFW = {1400,1000,1200,1800,1100,1100,1100,1100,1100,1600,1500,1100,1560,1560,1400,1400,1000,1100,1400,1000};
	private static final int[] DEFH = {700, 1000,900, 570, 1100,1100,1100,1100,1100,800, 1500,1100,780 ,780, 500, 700, 1000,1100,700 ,1000};
	private static final String[] TIP3 = {"An equidistant cylindrical map",
											"A conformal cylindrical map",
											"A compromising cylindrical map",
											"An equal-area cylindrical map",
											"An equidistant azimuthal map",
											"A conformal azimuthal map",
											"An equal-area azimuthal map",
											"Represents earth viewed from an infinite distance",
											"Every straight line on the map is a straight line on the sphere",
											"A conical map (conical maps suck; don't use this one)",
											"The compromise map used by National Geographic",
											"A circular compromise map",
											"An equal-area map shaped like an elipse",
											"An equal-area map shaped like an elipse",
											"An equal-area map shaped like a sinusoid",
											"BURN LIFE'S HOUSE DOWN!",
											"A conformal square map that uses complex math",
											"A novelty map that swells the center to disproportionate scale",
											"A reorganized version of Pierce Quincuncial and actually the best map ever",
											"A compromise map shaped like a square" };
	
	private static final String[] FILE = {"Realistic","Altitude","Sillouette","Rivers","HiContrast","Terrain",
			"No_Ice","Biomes","Satellite","Political","Timezones","Historic","Population","Antipode","Flights",
			"Empire","Mars","Stars","Color_Wheel","Grid","Circles","Soccer"};
	private static final String[] TIP1 = {"A realistic rendering of the Earth",
											"Color-coded based on altitude",
											"Land is black; oceans are white.",
											"Land is black; rivers and oceans are white (not recomended for low-resolution maps).",
											"Lots of snow and ice; oceans are black.",
											"Color-coded based on terrain",
											"Color-coded based on terrain, without ice",
											"Color-coded based on biome",
											"A composite satellite image of the earth",
											"Political map with country names removed",
											"A map of different time-zones",
											"An old map by European explorers",
											"Color-coded by population",
											"If you dug straight down, where would you end up?",
											"The most frequent airplane routes",
											"The British Empire at its height",
											"A realistic rendering of Mars",
											"The cosmos, as seen from Earth",
											"Color-coded by latitude and longitude",
											"Each square represents 30 degrees.",
											"10-degree Tissot's indatrix.",
											"A realistic rendering of a football" };
	
	private static final String[] AXES = {"Standard","Transverse","Center of Mass","Jerusalem",
			"Point Nemo","Longest Line","Longest Line Transverse","Cylindrical","Conical","Quincuncial" };
	private static final String[] TIP2 = {"The north pole (standard for most maps)",
											"Offset from standard by 90 degrees",
											"The center of landmass on Earth (Giza)",
											"The city of Jerusalem",
											"Antipode of the point farthest from land",
											"Sets the longest sailable line as the equator",
											"Sets the longest sailable line as the meridian",
											"Perfect for cylindrical maps",
											"Perfect for conical maps",
											"Perfect for the Pierce Quincuncial projection" };
	private static final double[] lats = {90,0,29.9792,31.7833,48.8767,-28.5217,-46.4883,-35,10,60};
	private static final double[] lons = {0,0,31.1344,35.216,56.6067,141.451,16.5305,-13.6064,-115,-6};
	private static final double[] thts = {0,180,-32,-35,-45,71.5,137,145,150,-10};
	
	
	public String command = "";
	
	
	
	
	public static void main(String[] args) {
		BufferedImage input, output;
		int w,h;
		double latD, lonD, thtD;
		int projection;
		
		MapProjections listener = new MapProjections(); // initialization
		JFrame frame = new JFrame("Map Configurer");
	    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    frame.setSize(400,300);
	    
	    while (true) { // make as many maps as you want
	
			JPanel panel = new JPanel();
			JLabel label = new JLabel("Please select a map theme."); // select map theme
			label.setToolTipText("This is the equirectangular image the program will reference.");
			panel.add(label);
			JButton buttn;
			for (int i = 0; i < FILE.length; i ++) {
			    buttn = new JButton(FILE[i]);
			    buttn.setToolTipText(TIP1[i]);
			    buttn.setActionCommand(FILE[i]);
			    buttn.addActionListener(listener);
			    panel.add(buttn);
			}
		    frame.add(panel);
		    frame.setVisible(true);
		    while (listener.isWaiting()) {} // waits for a button to be pressed
		    try {
		    	input = ImageIO.read(new File("input/"+listener.command+".jpg"));
		    } catch (IOException e) {
		    	System.err.println("Where the heck is the image?!");
		    	return;
		    }
		    listener.reset();
		    frame.remove(panel);
			
		    panel = new JPanel();
			label = new JLabel("Pick a projection algorithm."); // select projection
			label.setToolTipText("How will the Earth be mapped onto a plane?");
			panel.add(label);
			for (int i = 0; i < PROJ.length; i ++) {
			    buttn = new JButton(PROJ[i]);
			    buttn.setToolTipText(TIP3[i]);
			    buttn.setActionCommand(String.valueOf(i));
			    buttn.addActionListener(listener);
			    panel.add(buttn);
			}
		    frame.add(panel);
		    frame.setVisible(true);
		    while (listener.isWaiting()) {} // wait for a button to be pressed
		    projection = Integer.parseInt(listener.command);
		    listener.reset();
		    frame.remove(panel);
		    
		    panel = new JPanel();
			label = new JLabel("Choose an axis preset, or make a custom one."); // select axis
			label.setToolTipText("Changing the axis effectively rotates the earth, which can produce some very unconventional maps.");
			panel.add(label);
		    buttn = new JButton("Custom"); // custom button
		    buttn.setToolTipText("Enter coordinates to create a custom axis.");
		    buttn.setActionCommand("-2");
		    buttn.addActionListener(listener);
		    panel.add(buttn);
		    for (int i = 0; i < AXES.length; i ++) { // all the other buttons
		    	buttn = new JButton(AXES[i]);
		    	buttn.setToolTipText(TIP2[i]);
		    	buttn.setActionCommand(String.valueOf(i));
		    	buttn.addActionListener(listener);
		    	panel.add(buttn);
		    }
		    buttn = new JButton("Random"); // random button
		    buttn.setToolTipText("An axis will be chosen at random.");
		    buttn.setActionCommand("-1");
		    buttn.addActionListener(listener);
		    panel.add(buttn);
		    frame.add(panel);
		    frame.setVisible(true);
		    while (listener.isWaiting()) {} // wait for a button to be pressed
		    int n = Integer.parseInt(listener.command);
		    listener.reset();
		    frame.remove(panel);
		    
			if (n >= 0) { // if it is a preset
				latD = lats[n];
				lonD = lons[n];
				thtD = thts[n];
			}
			else if (n == -1) {
				latD = Math.toDegrees(Math.asin(Math.random()*2-1));
				lonD = Math.random()*360;
				thtD = Math.random()*360;
			}
			else { // if it is custom
				panel = new JPanel(); // lets you pick coordinates
				label = new JLabel("<html>Enter coordinates for your axis (lattitude, longitude, orientation).</html>");
				label.setToolTipText("The coordinates you specify will probably appear in the center or at the top of the map.");
				panel.add(label);
				SpinnerModel latModel = new SpinnerNumberModel(0, -90, 90, .1);
				JSpinner lat = new JSpinner(latModel);
				lat.setToolTipText("The lattitude of your desired axis (between -90 and 90)");
				SpinnerModel lonModel = new SpinnerNumberModel(0, -180, 180, .1);
				JSpinner lon = new JSpinner(lonModel);
				lon.setToolTipText("The longitude of your desired axis (between -180 and 180)");
				SpinnerModel thtModel = new SpinnerNumberModel(0, -180, 180, .1);
				JSpinner tht = new JSpinner(thtModel);
				tht.setToolTipText("The orientation of your desired axis (between -180 and 180)");
				panel.add(lat);
				panel.add(lon);
				panel.add(tht);
				buttn = new JButton("Okay");
				buttn.setToolTipText("Press when you are satisfied with your coordinates.");
				buttn.setActionCommand("OK");
				buttn.addActionListener(listener);
				panel.add(buttn);
				frame.add(panel);
				frame.setVisible(true);
				while (listener.isWaiting()) {} // wait for a button to be pressed
				latD = (double)lat.getValue();
				lonD = (double)lon.getValue();
				thtD = (double)tht.getValue();
				listener.reset();
				frame.remove(panel);
			}
		    
		    panel = new JPanel();
			label = new JLabel("Finally, set the dimensions for your map (width, height)."); // select map dimensions
			label.setToolTipText("These will be the dimensions of the JPG file in pixels.");
			panel.add(label);
			SpinnerModel widthModel = new SpinnerNumberModel(DEFW[projection], 400, 20000, 10);
			JSpinner width = new JSpinner(widthModel);
			width.setToolTipText("The width of your map in pixels");
			SpinnerModel heightModel = new SpinnerNumberModel(DEFH[projection], 400, 20000, 10);
			JSpinner height = new JSpinner(heightModel);
			height.setToolTipText("The height of your map in pixels");
			panel.add(width);
			panel.add(height);
		    buttn = new JButton("OK");
		    buttn.setToolTipText("Press when you are satisfied with your dimensions.");
		    buttn.setActionCommand("OK");
		    buttn.addActionListener(listener);
		    panel.add(buttn);
		    frame.add(panel);
		    frame.setVisible(true);
		    while (listener.isWaiting()) {} // wait for a button to be pressed
		    w = (int)(width.getValue());
		    h = (int)(height.getValue());
			output = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
		    listener.reset();
		    frame.remove(panel);
		    
		    panel = new JPanel();
			label = new JLabel("Wait...");
			panel.add(label);
			frame.add(panel);
			frame.setVisible(true);
			map(input,output,projection,latD,lonD,thtD); // This is where all the math happens
			
			saveImage(output);
			
			frame.remove(panel);
			panel = new JPanel();
			label = new JLabel("Done!"); // finished!
			panel.add(label);
			buttn = new JButton("Make Another");
			buttn.setToolTipText("You know you want to..."); // lets you start over
			buttn.setActionCommand("Go!");
			buttn.addActionListener(listener);
			panel.add(buttn);
			buttn = new JButton("Exit");
			buttn.setToolTipText("Please don't deactivate me!"); // lets you close the window
			buttn.setActionCommand("X");
			buttn.addActionListener(listener);
			panel.add(buttn);
			frame.add(panel);
			frame.setVisible(true);
			while (listener.isWaiting()) {}
			if (listener.command.equals("X"))	break; // exits if necessary
			frame.remove(panel);
			listener.reset();
		}
	    
	    frame.setVisible(false);
	    frame.dispose();
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
			return getColor(lat0,lon0,orientation, phi, Math.atan2(width/2.0-x, height/2.0-y), ref);
		else
			return 0;
	}
	
	
	public static int gall(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // a compromise map, similar to mercator
		return getColor(lat0, lon0, orientation,
			      	    2*Math.atan((y-height/2.0) / (height/2.0)), x*2*Math.PI/width, ref);
	}
	
	
	public static int sinusoidal(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // a map shaped like a sinusoid
		return getColor(lat0, lon0, orientation, y*Math.PI/height - Math.PI/2,
				        Math.PI * (x-width/2.0) / (Math.sin(Math.PI*y/height)*width/2.0)+Math.PI, ref);
	}
	
	
	public static int stereographic(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // a shape-preserving infinite map
		double radius = Math.pow(Math.pow(width, -2)+Math.pow(height, -2), -.5) / Math.PI;
		return getColor(lat0, lon0, orientation, 2*Math.atan(Math.hypot(x-width/2, y-height/2) / radius)-Math.PI/2,
                Math.atan2(width/2.0-x, height/2.0-y), ref);
	}
	
	
	public static int gnomic(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // a shape-preserving infinite map
		double radius = Math.pow(Math.pow(width, -2)+Math.pow(height, -2), -.5) / Math.PI;
		return getColor(lat0, lon0, orientation, Math.atan(Math.hypot(x-width/2, y-height/2) / radius)-Math.PI/2,
                Math.atan2(width/2.0-x, height/2.0-y), ref);
	}
	
	
	public static int orthogonal(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // a map that mimics the view from space
		double R = 2*Math.hypot((double)x/width-.5, (double)y/height-.5);
		if (R <= 1)
			return getColor(lat0, lon0, orientation, -Math.acos(R), Math.atan2(x-width/2.0,y-height/2.0), ref);
		else
			return 0;
	}
	
	
	public static int eaCylindrical(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // an equal area cylindrical map
		return getColor(lat0, lon0, orientation, Math.asin((y-height/2.0) / (height/2.0)),
                        x*2*Math.PI / width, ref);
	}
	
	
	public static int lambert(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // a conical projection
		double radius = Math.pow(Math.pow(width, -2)+Math.pow(height, -2), -.5) / Math.PI;
		return getColor(lat0, lon0, orientation, 4.0/3.0*(Math.atan(Math.hypot(x-width/2,y)/(radius)-1)+Math.PI/4) - Math.PI/2,
				        2*Math.atan2(width/2.0-x, -y)-Math.PI, ref);
	}
	
	
	public static int lemons(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // a simple map that is shaped like lemons
		int lemWdt;
		if (width > 12)  lemWdt= width/12; // the pixel width of each lemon
		else             lemWdt = width;
		
		if (Math.abs(x%lemWdt-lemWdt/2.0) < Math.sin(Math.PI*y/height)*lemWdt/2.0) // if it is inside a sin curve
		      return getColor(lat0,lon0,orientation, y*Math.PI/height - Math.PI/2,
		    		          Math.PI * (x%lemWdt-lemWdt/2.0) / (Math.sin(Math.PI*y/height)*lemWdt*6.0) + x/lemWdt*Math.PI/6, ref);
		else
			return 0;
	}
	
	
	public static int eaAzimuth(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // the lambert azimuthal equal area projection
		double R = 4*Math.hypot((double)x/width-.5, (double)y/height-.5);
		if (R <= 2)
			return getColor(lat0, lon0, orientation, Math.asin(R*R/2-1), Math.atan2(x-width/2.0, y-height/2.0), ref);
		else
			return 0;
	}
	
	
	public static int quinshift(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) { // a tessalatable rectangle map
		Complex u = new Complex(3.7116*(0.5*y/height + 1.0*x/width),
				                3.7116*(0.5*y/height - 1.0*x/width)); // don't ask me where 3.7116 come from because I have no idea
		Complex k = new Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus stuff
		Complex ans = Jacobi.cn(u, k);
		double p = 2*Math.atan(ans.abs());
		double theta = Math.atan2(ans.getRe(), ans.getIm());
		double lambda = p-Math.PI/2;
		return getColor(lat0, lon0, orientation, lambda, theta, ref);
	}
	
	
	public static int mollweide(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) {
		double tht = Math.asin(((double)y/height-.5)*2);
		return getColor(lat0,lon0,orientation, Math.asin((2*tht+Math.sin(2*tht))/Math.PI),
				Math.PI+2*Math.PI*((double)x/width-.5)/Math.cos(tht), ref);
	}
	
	
	public static int winkel_tripel(final double lat0, final double lon0, final double orientation,
            final int width, final int height, int x, int y, BufferedImage ref) {
		double phi = 0;
		double lam = 0;
		for (int i = 0; i < 1; i ++) {
			phi = phi - (wtX(phi,lam)-x)/wtXbyPhi(phi,lam);
			lam = lam - (wtX(phi,lam)-x)/wtXbyLam(phi,lam);
			phi = phi - (wtY(phi,lam)-y)/wtYbyPhi(phi,lam);
			lam = lam - (wtY(phi,lam)-y)/wtYbyLam(phi,lam);
		}
		return getColor(lat0,lon0,orientation, phi, lam, ref);
	}
	
	
	public static int grinten(final double lat0, final double lon0, final double orientation,
			final int width, final int height, int x, int y, BufferedImage ref) {
		if (y == height/2) // special case 1: equator
			return getColor(lat0,lon0,orientation, 0, 2*Math.PI*x/width, ref);
		
		double X = 2.0*x/width - 1;
		double Y = 2.0*y/height - 1;
		
		if (x == width/2) // special case 3: meridian
			return getColor(lat0,lon0,orientation, Math.PI/2*Math.sin(2*Math.atan(Y)), Math.PI, ref);
		
		double c1 = -Math.abs(Y)*(1 + X*X + Y*Y);
		double c2 = c1 - 2*Y*Y + X*X;
		double c3 = -2*c1 + 1 + 2*Y*Y + Math.pow(X*X+Y*Y, 2);
		double d = Y*Y/c3 + 1/27.0*(2*Math.pow(c2/c3, 3) - 9*c1*c2/(c3*c3));
		double a1 = 1/c3*(c1-c2*c2/(3*c3));
		double m1 = 2*Math.sqrt(-a1/3);
		double tht1 = Math.acos(3*d/(a1*m1))/3;
		
		return getColor(lat0,lon0,orientation, Math.signum(Y)*Math.PI*(-m1*Math.cos(tht1+Math.PI/3)-c2/(3*c3)),
				Math.PI*(X*X + Y*Y - 1 + Math.sqrt(1 + 2*(X*X-Y*Y)+Math.pow(X*X+Y*Y, 2)))/(2*X) + Math.PI, ref);
	}
	
	
	public static int custom1(final double lat0, final double lon0, final double orientation,
			final int width, final int height, int x, int y, BufferedImage ref) { // a tessalatable square map
		final double xp = 2.0*x/width;
		final double yp = 2.0*y/height-1;
		final double a = .1; // this is arbitrary, but should be between 0 and .2
		double X = xp;
		double Y = yp; // these must be set with Newton's method (I think)
		
		for (int i = 0; i < 1; i ++) { // honestly I can't tell the difference between 1 and 10 iterations, so my equation might be somehow linear
			X = X - (a*Math.sin(Math.PI*X)*Math.pow(Math.cos(Math.PI*yp),5)-X+xp)/(a*Math.PI*Math.cos(Math.PI*X)*Math.pow(Math.cos(Math.PI*yp),5)-1);
			Y = Y - (a*Math.sin(Math.PI*Y)*Math.pow(Math.cos(Math.PI*xp),5)-Y+yp)/(a*Math.PI*Math.cos(Math.PI*Y)*Math.pow(Math.cos(Math.PI*xp),5)-1);
		}
		
		Complex u = new Complex(1.8558*X, 1.8558*Y); // don't ask me where 3.7116 come from because I have no idea
		Complex k = new Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus stuff
		Complex ans = Jacobi.cn(u, k);
		double p = 2*Math.atan(ans.abs());
		double theta = Math.atan2(ans.getRe(), ans.getIm());
		double lambda = p-Math.PI/2;
		return getColor(lat0, lon0, orientation, lambda, theta, ref);
	}
	
	
	public static int custom2(final double lat0, final double lon0, final double orientation,
			final int width, final int height, int x, int y, BufferedImage ref) { // a tessalatable square map
		final double X = Math.PI*x/width-Math.PI/2;
		final double Y = Math.PI*y/height-Math.PI/2;
		final double a = .1; // this is arbitrary, but should be between 0 and .2
		double mx = 1 + a*Math.sqrt(1-Y*Y);
		double my = 1 + a*Math.sqrt(1-X*X);
		if (Double.isNaN(mx))	mx = 1;
		if (Double.isNaN(my))	my = 1;
		double xp = Math.asin(Math.sin(X*mx)/mx);
		double yp = Math.asin(Math.sin(Y*my)/my);
		
		Complex u = new Complex(1.18*(xp+Math.PI/2), 1.18*(yp));
		Complex k = new Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus stuff
		Complex ans = Jacobi.cn(u, k);
		double p = 2*Math.atan(ans.abs());
		double theta = Math.atan2(ans.getRe(), ans.getIm());
		double lambda = p-Math.PI/2;
		return getColor(lat0, lon0, orientation, lambda, theta, ref);
	}
	
	
	public static int magnus(final double lat0, final double lon0, final double orientation,
			final int width, final int height, int x, int y, BufferedImage ref) { // a novelty map that magnifies the center profusely
		double r = 2*Math.hypot((double)x/width-.5, (double)y/height-.5);
		if (r <= 1)
			return getColor(lat0, lon0, orientation, Math.PI/2*(r*r*r*1.8+r*.2-1), Math.atan2(x-width/2.0, y-height/2.0), ref);
		else
			return 0;
	}
	
	
	public static int hammer(final double lat0, final double lon0, final double orientation,
			final int width, final int height, int x, int y, BufferedImage ref) { // similar to Mollweide, but moves distortion from the poles to the edges
		final double X = 4*Math.sqrt(2)*x/width-2*Math.sqrt(2);
		final double Y = 2*Math.sqrt(2)*y/height-Math.sqrt(2);
		final double z = Math.sqrt(1 - Math.pow(X/4,2) - Math.pow(Y/2,2));
		return getColor(lat0, lon0, orientation, Math.asin(z*Y), Math.PI+2*Math.atan(.5*z*X/(2*z*z-1)), ref);
	}
	/*END PROJECTION METHODS*/
	
	
	public static int getColor(final double lat0, final double lon0, final double orientation,
			                   double lat1, double lon1, BufferedImage ref) { // returns the color of any coordinate on earth		
		lon1 += orientation;
		double latitude = Math.asin(Math.sin(lat0)*Math.sin(lat1) + Math.cos(lat0)*Math.cos(lon1)*Math.cos(lat1));
		double longitude;
		double innerFunc = Math.sin(lat1)/Math.cos(lat0)/Math.cos(latitude)-Math.tan(lat0)*Math.tan(latitude); // used for calculating lon
		if (lat0  == Math.PI/2)				// accounts for special case when lat0 = pi/2
			longitude = lon1+Math.PI;
		else if (lat0 == -Math.PI/2)		// accounts for special case when lat0 = -pi/2
			longitude = -lon1;
		else if (Math.abs(innerFunc) > 1) {	// accounts for special case when cos(lat) = --> 0
			if ((lon1 == Math.PI && lat1 < -lat0) || (lon1 != Math.PI && lat1 < lat0))
				longitude = Math.PI+lon0;
			else
				longitude = lon0;
		}
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
				switch (PROJ[projection]) { // methods are selected by the name of the projection
				case "Pierce Quincuncial":
					output.setRGB(x, y, quincuncial(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Equirectangular":
					output.setRGB(x, y, equirectangular(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Mercator":
					output.setRGB(x, y, mercator(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Polar":
					output.setRGB(x, y, polar(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Gall Stereographic":
					output.setRGB(x, y, gall(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Sinusoidal":
					output.setRGB(x, y, sinusoidal(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Stereographic":
					output.setRGB(x, y, stereographic(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Orthogonal":
					output.setRGB(x, y, orthogonal(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Lemons":
					output.setRGB(x, y, lemons(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Azimuthal Equal-Area":
					output.setRGB(x, y, eaAzimuth(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Cylindrical Equal-Area":
					output.setRGB(x, y, eaCylindrical(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Lambert Conical":
					output.setRGB(x, y, lambert(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Gnomic":
					output.setRGB(x, y, gnomic(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Guyou Hemisphere-in-a-Square":
					output.setRGB(x, y, quinshift(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Mollweide":
					output.setRGB(x, y, mollweide(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Winkel Tripel":
					output.setRGB(x, y, winkel_tripel(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Van der Grinten":
					output.setRGB(x, y, grinten(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Rectus Aequilibrium":
					output.setRGB(x, y, custom2(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Hammer":
					output.setRGB(x, y, hammer(lat0,lon0,tht0,width,height,x,y,input));
					break;
				case "Magnifier":
					output.setRGB(x, y, magnus(lat0,lon0,tht0,width,height,x,y,input));
					break;
				default:
					System.err.println("Justin, you forgot to add a projection to the switch case! (or you forgot a break;)");
				}
			}
		}
	}
	
	
	private static void saveImage(BufferedImage img) {
		try {
			File outputFile = new File("output/myMap.jpg");
			ImageIO.write(img, "jpg", outputFile);
			Desktop.getDesktop().open(outputFile);
		} catch (IOException e) {}
	}
	
	
	public static final double wtX(double phi, double lam) {
		return 2/Math.PI*lam + 2*Math.cos(phi)*Math.sin(lam/2)*csccalpha(phi,lam);
	}
	
	public static final double wtXbyPhi(double phi, double lam) {
		return -2*Math.sin(phi)*Math.sin(lam/2)*csccalpha(phi,lam) + 2*Math.cos(phi)*Math.sin(lam/2)*csccalphaByPhi(phi,lam);
	}
	
	public static final double wtXbyLam(double phi, double lam) {
		return 2/Math.PI + Math.cos(phi)*Math.cos(lam/2)*csccalpha(phi,lam) + 2*Math.cos(phi)*Math.sin(lam/2)*csccalphaByLam(phi,lam);
	}
	
	public static final double wtY(double phi, double lam) {
		return phi + Math.sin(phi)*csccalpha(phi,lam);
	}
	
	public static final double wtYbyPhi(double phi, double lam) {
		return 1 + Math.cos(phi)*csccalpha(phi,lam) + Math.sin(phi)*csccalphaByPhi(phi,lam);
	}
	
	public static final double wtYbyLam(double phi, double lam) {
		return Math.sin(phi)*csccalphaByLam(phi,lam);
	}
	
	public static final double csccalpha(double phi, double lam) {
		return Math.acos(Math.cos(phi)*Math.cos(lam/2))/Math.sqrt(1-Math.pow(Math.cos(phi),2)*Math.pow(Math.cos(lam/2),2));
	}
	
	public static final double csccalphaByPhi(double phi, double lam) {
		return (Math.sin(phi)*Math.cos(lam/2) - Math.acos(Math.cos(phi)*Math.cos(lam/2))*Math.pow(Math.cos(lam/2),2)*Math.sin(phi)*Math.cos(phi)*Math.pow(1-Math.pow(Math.cos(phi),2)*Math.pow(Math.cos(lam/2),2),-1.5))
				/ (1 - Math.pow(Math.cos(phi),2)*Math.pow(Math.cos(lam/2),2));
	}
	
	public static final double csccalphaByLam(double phi, double lam) {
		return (Math.sin(lam/2)*Math.cos(phi) - Math.acos(Math.cos(lam/2)*Math.cos(phi))*Math.pow(Math.cos(phi),2)*Math.sin(lam/2)*Math.cos(lam/2)*Math.pow(1-Math.pow(Math.cos(lam/2),2)*Math.pow(Math.cos(phi),2),-1.5))
				/ (2 - 2*Math.pow(Math.cos(phi),2)*Math.pow(Math.cos(lam/2),2));
	}
	
	
	
	
	
	
	public void actionPerformed(ActionEvent e) { // the non-static part of the program acts as a button-listener
		command = e.getActionCommand();
	}
	
	
	public void reset() {
		command = "";
	}
	
	
	public boolean isWaiting() {
		if (command.isEmpty()) {
			System.out.print(""); // this line makes the code work. I've no idea why.
			return true;
		}
		else
			return false;
	}
}