import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class WormholeSimulator {

	public static void main(String[] args) {
		final double r = 1;
		final double l = 4;
		final double rs = 0.1;
		final double speed = 0.03;
		final double d0 = 8;
		final double df = 2;
		final int w = 1920;
		final int h = 1080;
		
		int i = 0;
		double dist = (d0-r)/speed + l/speed + (df-r)/speed;
		
		BufferedImage[] input = new BufferedImage[3];
		try {
	    	input[0] = ImageIO.read(new File("wormhole/input0.jpg"));
	    	input[1] = ImageIO.read(new File("wormhole/input1.jpg"));
	    	input[2] = ImageIO.read(new File("wormhole/input2.jpg"));
	    } catch (IOException e) {
	    	System.err.println("Where the heck is the image?!");
	    	return;
	    }
		
		for (double d = d0; d > r; d -= speed) {
			createFrame(r,l,rs,d,w,h,-1,i,input);
			System.out.println((int)(100*i/dist)+"%");
			i ++;
		}
		for (double d = 0; d < l; d += speed) {
			createFrame(r,l,rs,d,w,h,0,i,input);
			System.out.println((int)(100*i/dist)+"%");
			i ++;
		}
		for (double d = r; d < df; d += speed) {
			createFrame(r,l,rs,d,w,h,1,i,input);
			System.out.println((int)(100*i/dist)+"%");
			i ++;
		}
	}
	
	
	public static void createFrame(double r, double l, double rs, double d, int w, int h, int pos, int i, BufferedImage[] input) {
		BufferedImage frame = new BufferedImage(w,h,BufferedImage.TYPE_INT_RGB);
		for (int x = 0; x < w; x ++)
			for (int y = 0; y < h; y ++)
				setPixel(frame, x, y, r, l, rs, d, w, h, pos, input);
		
		try {
			File outputFile = new File("wormhole/frame"+String.format("%04d", i)+".jpg");
			ImageIO.write(frame, "jpg", outputFile);
		} catch (IOException e) {}
	}

	private static void setPixel(BufferedImage output, int x, int y, double r, double l, double rs, double d,
								 int w, int h, int pos, BufferedImage[] input) {
		double phi = 2*Math.PI*x/w;				// directional angle
		double thtI = Math.PI*y/h;	// azimuthal angle of incidence
		double thtR;					// azimuthal angle of... refraction? reflection? rewormholation?
		int ndx;
		
		if (pos == -1) {	// outside of the worm-hole
			if (Math.abs(thtI) < Math.asin(r/d)) {	// looking into the worm-hole
				thtR = -2.0*Math.asin(d/r*Math.sin(thtI)) - l*Math.tan(Math.asin(d/r*Math.sin(thtI))) - thtI;
				ndx = 1;
			}
			else {	// looking past the worm-hole
				thtR = thtI;
				ndx = 0;
			}
		}
		else if (pos == 0) {	// inside of the worm-hole
			if (Math.abs(thtI) < Math.PI/2-Math.asin(rs/(2*Math.PI*r))) {
				thtR = (l-d)/r*Math.tan(thtI) + thtI;
				ndx = 1;
			}
			else if (Math.abs(thtI) > Math.PI/2+Math.asin(rs/(2*Math.PI*r))) {
				thtR = -d/r*Math.tan(thtI)+(Math.PI-thtI);
				ndx = 0;
			}
			else {
				thtR = 0;
				ndx = 2;
			}
		}
		else { // past the worm-hole
			thtR = thtI;
			ndx = 1;
		}
		if (((int)((thtR/(Math.PI))%2)+2)%2 == 1) // phi is inverted for certain thtRs
			phi = phi+Math.PI;
		phi = (phi+2*Math.PI)%(2*Math.PI);
		thtR = Math.acos(Math.cos(thtR)); // puts it in the range 0 to pi
		
		int x0 = (int)(phi/(2*Math.PI)*input[ndx].getWidth());	// reference coordinates on the image
		int y0 = (int)(thtR*input[ndx].getHeight()/(Math.PI));	// proportional to latitude and longitude
		if (x0 >= input[ndx].getWidth()) // coerces x0 and y0 into bounds
			x0 = input[ndx].getWidth()-1;
		if (x0 < 0)
			x0 = 0;
		if (y0 >= input[ndx].getHeight())
			y0 = input[ndx].getHeight()-1;
		if (y0 < 0)
			y0 = 0;
		
		output.setRGB(x, y, input[ndx].getRGB(x0,y0));
	}

}
