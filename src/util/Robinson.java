package util;

import java.util.Arrays;

public class Robinson {
	private static final double[][] TABLE = {
			{ 00, 05, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90,
				-05,-10,-15,-20,-25,-30,-35,-40,-45,-50,-55,-60,-65,-70,-75,-80,-85,-90 },
			{ 1.000, 0.9986, 0.9954, 0.990, 0.9822, 0.9730, 0.960, 0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322,
					0.9986, 0.9954, 0.990, 0.9822, 0.9730, 0.960, 0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322 },
			{ 0.000, 0.0620, 0.1240, 0.186, 0.2480, 0.3100, 0.372, 0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000,
					-0.0620,-0.1240,-0.186,-0.2480,-0.3100,-0.372,-0.4340,-0.4958,-0.5571,-0.6176,-0.6769,-0.7346,-0.7903,-0.8435,-0.8936,-0.9394,-0.9761,-1.0000 }
	};
	
	private static final int LEN = 19;
	
	
	public static final double plenFromLat(double lat) {
		if (Math.abs(lat) < Math.PI/3)
			return aitkenInterpolate(Math.toDegrees(lat), TABLE[0], TABLE[1]);
		else
			return aitkenInterpolate(Math.toDegrees(lat),
					Arrays.copyOfRange(TABLE[0], LEN-10,LEN-1),
					Arrays.copyOfRange(TABLE[1], LEN-10, LEN-1));
	}
	
	
	public static final double pdfeFromLat(double lat) {
		if (Math.abs(lat) < Math.PI/3)
			return aitkenInterpolate(Math.toDegrees(lat), TABLE[0], TABLE[2]);
		else
			return aitkenInterpolate(Math.toDegrees(lat),
					Arrays.copyOfRange(TABLE[0], LEN-10,LEN-1),
					Arrays.copyOfRange(TABLE[2], LEN-10, LEN-1));
	}
	
	
	public static final double latFromPdfe(double pdfe) {
		if (pdfe < TABLE[2][12])
			return Math.toRadians(aitkenInterpolate(pdfe, TABLE[2], TABLE[0]));
		else
			return Math.toRadians(aitkenInterpolate(pdfe,
					Arrays.copyOfRange(TABLE[2], LEN-10,LEN-1),
					Arrays.copyOfRange(TABLE[0], LEN-10, LEN-1)));
	}
	
	
	public static final double plenFromPdfe(double pdfe) {
		if (pdfe < TABLE[2][12])
			return aitkenInterpolate(pdfe, TABLE[2], TABLE[1]);
		else
			return aitkenInterpolate(pdfe,
					Arrays.copyOfRange(TABLE[2], LEN-15,LEN-1),
					Arrays.copyOfRange(TABLE[1], LEN-15, LEN-1));
	}
	
	
	public static final double aitkenInterpolate(double x,
			double[] X, double[] f) { //map from ith column to jth using aitken interpolation
		final int N = X.length;
		final double[][] fx = new double[N][]; // the table of successive approximations
		
		fx[0] = new double[N]; //fill in the zeroth row
		for (int j = 0; j < N; j ++) {
			fx[0][j] = f[j]; //the initial approximations are horizontal lines
		}
		for (int i = 1; i < N; i ++) { //i+1 is the number of points interpolated on
			fx[i] = new double[N];
			for (int j = i; j < N; j ++) { //the points will be 0, ..., i-1, j
				fx[i][j] = 1/(X[j] - X[i-1])*determ(fx[i-1][i-1], fx[i-1][j],
														X[i-1] - x, X[j] - x); //I won't bother to explain this; go look up Aitken interpolation
			}
		}
		
		return fx[N-1][N-1];
	}
	
	
	private static final double determ(double a, double b, double c, double d) {
		return a*d - b*c;
	}
	
	
	public static final void main(String[] args) {
		System.out.println("Testing Aitken Interpolation:");
		System.out.println("Input points are (-1,1), (-.5,-1), (0,0), (.5,1), and (1,-1)");
		//double[] X = {-1, -.5, 0, .5, 1};
		//double[] Y = {1, -1, 0, 1, -1};
		double[] X = TABLE[2];
		double[] Y = TABLE[1];
		System.out.println("Interpolating points from -1 to 1:");
		for (double x = -1; x <= 1; x += 1/1024.)
			System.out.println(x+", "+aitkenInterpolate(x, X, Y)+";");
	}
}
