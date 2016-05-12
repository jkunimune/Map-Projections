
public final class WinkelTripel {
	public static final double X(double phi, double lam) {
		return 2/Math.PI*lam + 2*Math.cos(phi)*Math.sin(lam/2)/sinca(phi,lam);
	}
	
	public static final double Y(double phi, double lam) {
		return phi + Math.sin(phi)/sinca(phi,lam);
	}
	
	public static final double dXdphi(double phi, double lam) {
		final double alpha = a(phi,lam);
		final double beta = b(phi,lam);
		return 2*Math.sin(phi)*Math.sin(lam/2)*alpha/beta - 
				2*Math.sin(phi)*Math.cos(phi)*Math.sin(lam/2)*Math.cos(lam/2)/(beta*beta) + 
				2*Math.sin(phi)*Math.pow(Math.cos(phi), 2)*Math.sin(lam/2)*Math.pow(Math.cos(lam/2), 2)/(beta*beta*beta);
	}
	
	public static final double dXdlam(double phi, double lam) {
		final double alpha = a(phi,lam);
		final double beta = b(phi,lam);
		return 2/Math.PI + Math.cos(phi)*Math.cos(lam/2)*alpha/beta - 
				Math.pow(Math.cos(phi)*Math.sin(lam/2)/beta, 2) + 
				Math.pow(Math.cos(phi), 3)*Math.pow(Math.sin(lam/2), 2)*Math.cos(lam/2)*alpha/(beta*beta*beta);
	}
	
	public static final double dYdphi(double phi, double lam) {
		final double alpha = a(phi,lam);
		final double beta = b(phi,lam);
		return 1 + Math.cos(phi)*alpha/beta - 
				Math.pow(Math.sin(phi), 2)*Math.cos(lam/2)/(beta*beta) + 
				Math.pow(Math.sin(phi), 2)*Math.cos(phi)*Math.pow(Math.cos(lam/2), 2)*alpha/(beta*beta*beta);
	}
	
	public static final double dYdlam(double phi, double lam) {
		final double alpha = a(phi,lam);
		final double beta = b(phi,lam);
		return Math.sin(phi)*Math.pow(Math.cos(phi), 2)*Math.sin(lam/2)*Math.cos(lam/2)*alpha/(2*beta*beta*beta) - 
				Math.sin(phi)*Math.cos(phi)*Math.sin(lam/2)/2*(beta*beta);
	}
	
	public static final double a(double phi, double lam) {
		return Math.acos(Math.cos(phi)*Math.cos(lam/2));
	}
	
	public static final double b(double phi, double lam) {
		return Math.sqrt(1-Math.pow(Math.cos(phi)*Math.cos(lam/2), 2));
	}
	
	public static final double sinca(double phi, double lam) {
		return b(phi,lam)/a(phi,lam);
	}
}
