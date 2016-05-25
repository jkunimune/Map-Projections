
public final class WinkelTripel {
	public static final double X(double phi, double lam) {
		return 2/Math.PI*lam + 2*Math.cos(phi)*Math.sin(lam/2)/sinca(phi,lam);
	}
	
	public static final double Y(double phi, double lam) {
		return phi + Math.sin(phi)/sinca(phi,lam);
	}
	
	public static final double dXdphi(double phi, double lam) {
		final double s = sinca(phi,lam);
		final double dsdP = dsincadphi(phi,lam);
		return (2*Math.sin(phi)*Math.sin(lam/2)*s - 2*Math.cos(phi)*Math.sin(lam/2)*dsdP)/(s*s);
	}
	
	public static final double dXdlam(double phi, double lam) {
		final double s = sinca(phi,lam);
		final double dsdL = dsincadlam(phi,lam);
		return 2/Math.PI + (Math.cos(phi)*Math.cos(lam/2)*s - 2*Math.cos(phi)*Math.sin(lam/2)*dsdL)/(s*s);
	}
	
	public static final double dYdphi(double phi, double lam) {
		final double s = sinca(phi,lam);
		final double dsdP = dsincadphi(phi,lam);
		return 1 + (Math.cos(phi)*s-Math.sin(phi)*dsdP)/(s*s);
	}
	
	public static final double dYdlam(double phi, double lam) {
		final double s = sinca(phi,lam);
		final double dsdL = dsincadlam(phi,lam);
		return -Math.sin(phi)*dsdL/(s*s);
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
	
	public static final double dsincadphi(double phi, double lam) {
		final double alpha = a(phi,lam);
		return (Math.cos(alpha)*alpha-Math.sin(alpha))/(alpha*alpha)*dadphi(phi,lam);
	}
	
	public static final double dsincadlam(double phi, double lam) {
		final double alpha = a(phi,lam);
		return (Math.cos(alpha)*alpha-Math.sin(alpha))/(alpha*alpha)*dadlam(phi,lam);
	}
	
	public static final double dadphi(double phi, double lam) {
		return -Math.sin(phi)*Math.cos(lam/2)/b(phi,lam);
	}
	
	public static final double dadlam(double phi, double lam) {
		return -Math.cos(phi)*Math.sin(lam/2)/(2*b(phi,lam));
	}
}
