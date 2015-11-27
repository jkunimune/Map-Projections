public final class Vector {
 public static final Vector I = new Vector(1,0,0,true);
 public static final Vector J = new Vector(0,1,0,true);
 public static final Vector K = new Vector(0,0,1,true);
 
 
 private double r; // magnitude
 private double a; // altitude from horizontal
 private double b; // bearing
  
  
  
  public Vector(double newX, double newY, double newZ, boolean cartesian) { // constructs a new vector given horizontal, vertical, and depthual lengths
    r = Math.sqrt(newX*newX + newY*newY + newZ*newZ);
    a = Math.asin(newZ/r); // Z is the positive direction
    b = Math.atan2(newY,newX);
  }
  
  
  public Vector(double newR, double newAlpha, double newBeta) { // constructs a new vector given magnitude, altitude, and bearing
    r = newR;
    a = newAlpha;
    b = newBeta;
  }
  
  
  public Vector(int latD, int latM, int lonD, int lonM) { // constructs a vector of magnitude 1 given degree-minute coordinates
    r = 1;
    a = Math.toRadians(latD+latM/60.0);
    b = Math.toRadians(lonD+lonM/60.0);
  }  
  
  
  
  public final void setR(double newR) {
    r = newR;
  }
  
  
  public final void setA(double newAlpha) {
    a = newAlpha;
  }
  
  
  public final void setB(double newBeta) {
    b = newBeta;
  }
  
  
  public final double getX() { // magnitude of the width component
    return r*Math.cos(a)*Math.cos(b);
  }
  
  
  public final double getY() { // magnitude of the depth component
   return r*Math.cos(a)*Math.sin(b);
  }
  
  
  public final double getZ() { // magnitude of the height component
   return r*Math.sin(a);
  }
  
  
  public final double getR() { // magnitude
    return r;
  }
  
  
  public final double getA() { // altitude
    return a;
  }
  
  
  public final double getB() { // bearing
    return b;
  }
  
  
  public final Vector negative() { // computes the opposite
    return new Vector(r, -a, (2*Math.PI+b)%(2*Math.PI));
  }
  
  
  public final Vector plus(Vector that) { // computes sum with that
    return new Vector(this.getX()+that.getX(), this.getY()+that.getY(), this.getZ()+that.getZ(), true);
  }
  
  
  public final Vector minus(Vector that) { // computes difference with that
    return new Vector(this.getX()-that.getX(), this.getY()-that.getY(), this.getZ()-that.getZ(), true);
  }
  
  
  public final Vector times(double c) { // computes product with c
    return new Vector(c*r, a, b);
  }
  
  
  public final double dot(Vector that) { // computes dot product with that
    return this.getX()*that.getX() + this.getY()*that.getY() + this.getZ()*that.getZ();
  }
  
  
  public final Vector cross(Vector that) { // computes cross product with that
    return new Vector(this.getY()*that.getZ() - this.getZ()*that.getY(),
                      this.getZ()*that.getX() - this.getX()*that.getZ(),
                      this.getX()*that.getY() - this.getY()*that.getX(), true);
  }
  
  
  public final Vector hat() { // makes the magnitude 1
   return new Vector(1, getA(), getB());
  }
  
  
  public final void negate() { // negates
    a = -a;
    b = (2*Math.PI+b)%(2*Math.PI);
  }
  
  
  public final void plusEquals(Vector that) { // adds that
    Vector sum = this.plus(that);
    r = sum.getR();
    a = sum.getA();
    b = sum.getB();
  }
  
  
  public final void minusEquals(Vector that) { // subtracts that
    Vector dif = this.minus(that);
    r = dif.getR();
    a = dif.getA();
    b = dif.getB();
  }
  
  
  public final void timesEquals(double c) { // multiplies by c
    r *= c;
  }
  
  
  public final void crossEquals(Vector that) { // becomes cross product with that
    Vector txt = this.cross(that);
    r = txt.getR();
    a = txt.getA();
    b = txt.getB();
  }
  
  
  public final String toString() {
    return "<"+getX()+", "+getY()+", "+getZ()+">";
  }
  
  
  public final String toStringPolar() {
    return "("+getR()+", "+getA()+", "+getB()+")";
  }
}