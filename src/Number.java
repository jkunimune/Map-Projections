public class Number {
	private double r;
	private double theta;
	
	public static final Number ONE = new Number(1,0);

	
	
	
	public Number(double mag, double ang) {
		r = mag;
		theta = ang;
	}

	
	public Number(double x, double y, boolean cartesian) {
		r = Math.hypot(x, y);
		theta = Math.atan2(y, x);
	}
	
	
	
	
	public Number plus(Number that) {
		return new Number (this.x()+that.x(), this.y()+that.y(), true);
	}
	
	
	public Number plus(double x) {
		return new Number (this.x()+x, this.y(), true);
	}
	
	
	public Number minus(Number that) {
		return new Number (this.x()-that.x(), this.y()-that.y(), true);
	}
	
	
	public Number minus(double x) {
		return new Number (this.x()-x, this.y(), true);
	}
	
	
	public Number times(Number that) {
		return new Number (this.r()*that.r(), this.theta()+that.theta());
	}
	
	
	public Number times(double x) {
		return new Number (r*x, theta);
	}
	
	
	public Number over(double x) {
		return new Number (r/x, theta);
	}
	
	
	public Number sqrd() { // squares it
		return new Number (r*r, theta*2);
	}
	
	
	public Number timesI() { // multiplies it by i
		return new Number (-y(), x());
	}
	
	
	public Number overI() { // multiplies by negative I
		return new Number (y(), -x());
	}
	
	
	public static Number sqrt(Number z) { // sqrts it
		return new Number (Math.sqrt(z.r()), z.theta()/2);
	}
	
	
	public static Number ln(Number z) {
		return new Number (Math.log(z.r()), z.theta(), true);
	}
	
	
	public static Number asin(Number z) {
		return Number.ln(z.timesI().plus(Number.sqrt(Number.ONE.minus(z.sqrd())))).overI();
	}
	
	
	public static double abs(Number z) {
		return z.r();
	}
	
	
	public double r() {
		return r;
	}
	
	
	public double theta() {
		return theta;
	}
	
	
	public double x() {
		return r*Math.cos(theta);
	}
	
	
	public double y() {
		return r*Math.sin(theta);
	}
}