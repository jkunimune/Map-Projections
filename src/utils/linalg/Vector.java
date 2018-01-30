/**
 * MIT License
 * 
 * Copyright (c) 2017 Justin Kunimune
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package utils.linalg;


/**
 * A column vector to which you can do linear algebra.
 * This class isn't very efficient, but I'm using it on vectors of max length 3, so
 * 
 * @author jkunimune
 */
public class Vector extends Matrix {
	
	public Vector(int n) {
		super(n, 1);
	}
	
	public Vector(double... values) {
		super(values.length, 1);
		for (int i = 0; i < values.length; i ++)
			this.setElement(i, 0, values[i]);
	}
	
	private Vector(double[][] values) {
		super(values);
		if (values[0].length != 1)
			throw new IllegalArgumentException("Matrix has width "+values[0].length+" and can therefore not be converted to a column vector.");
	}
	
	public static Vector fromMatrix(Matrix mat) {
		return new Vector(mat.getArray());
	}
	
	public static Vector unit(int i, int n) {
		Vector iHat = new Vector(n);
		iHat.setElement(i, 1);
		return iHat;
	}
	
	
	public int getLength() {
		return this.getHeight();
	}
	
	public double getElement(int i) {
		return this.getElement(i, 0);
	}
	
	public void setElement(int i, double val) {
		this.setElement(i, 0, val);
	}
	
	public double[] asArray() {
		double[] arr = new double[this.getLength()];
		for (int i = 0; i < arr.length; i ++)
			arr[i] = this.getElement(i, 0);
		return arr;
	}
	
	public Vector plus(Vector that) {
		return new Vector(super.plus(that).getArray());
	}
	
	public Vector minus(Vector that) {
		return new Vector(super.minus(that).getArray());
	}
	
	public Vector times(double c) {
		return new Vector(super.times(c).getArray());
	}
	
	public double dot(Vector that) {
		return this.transpose().times(that).getElement(0, 0);
	}
	
	public double norm() {
		double s = 0;
		for (int i = 0; i < this.getLength(); i ++)
			s += Math.pow(this.getElement(i), 2);
		return Math.sqrt(s);
	}
	
	public Vector hat() {
		return this.times(1/this.norm());
	}
	
}
