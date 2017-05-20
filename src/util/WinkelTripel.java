package util;


public final class WinkelTripel {

	private static final double COS_PHI0 = 2/Math.PI;
	
	
	public static final double f1pX(double phi, double lam) {
		final double d = D(phi,lam);
		final double c = C(phi,lam);
		return (2*d/Math.sqrt(c)*Math.cos(phi)*Math.sin(lam/2) + lam*COS_PHI0)/2.0;
	}
	
	public static final double f2pY(double phi, double lam) {
		final double d = D(phi,lam);
		final double c = C(phi,lam);
		return (d/Math.sqrt(c)*Math.sin(phi) + phi)/2.0;
	}
	
	public static final double df1dphi(double phi, double lam) {
		final double d = D(phi,lam);
		final double c = C(phi,lam);
		return Math.sin(lam)*Math.sin(2*phi)/(4*c) - d/Math.pow(c,1.5)*Math.sin(phi)*Math.sin(lam/2);
	}
	
	public static final double df1dlam(double phi, double lam) {
		final double d = D(phi,lam);
		final double c = C(phi,lam);
		return (Math.pow(Math.cos(phi)*Math.sin(lam/2), 2)/c + d/Math.pow(c,1.5)*Math.cos(phi)*Math.cos(lam/2)*Math.pow(Math.sin(phi),2) + COS_PHI0)/2.0;
	}
	
	public static final double df2dphi(double phi, double lam) {
		final double d = D(phi,lam);
		final double c = C(phi,lam);
		return (Math.pow(Math.sin(phi),2)*Math.cos(lam/2)/c + d/Math.pow(c,1.5)*(1-Math.pow(Math.cos(lam/2),2))*Math.cos(phi) + 1)/2.0;
	}
	
	public static final double df2dlam(double phi, double lam) {
		final double d = D(phi,lam);
		final double c = C(phi,lam);
		return (Math.sin(2*phi)*Math.sin(lam/2)/c - d/Math.pow(c,1.5)*Math.sin(phi)*Math.pow(Math.cos(phi),2)*Math.sin(lam))/8.0;
	}
	
	private static final double D(double phi, double lam) {
		return Math.acos(Math.cos(phi)*Math.cos(lam/2));
	}
	
	private static final double C(double phi, double lam) {
		return 1 - Math.pow(Math.cos(phi)*Math.cos(lam/2), 2);
	}

}
