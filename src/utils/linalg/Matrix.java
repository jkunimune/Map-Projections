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
 * A two-dimensional array of numbers to which you can do linear algebra
 * FYI these Matrices index from zero because iatshuip MatLab
 * This class isn't very efficient, but I'm using it on matrices of max length 3, so
 * 
 * @author jkunimune
 */
public class Matrix {
	
	private final double[][] values;
	
	
	public Matrix(int n, int m) {
		this.values = new double[n][m];
	}
	
	public Matrix(double[]... values) {
		this.values = values;
	}
	
	public static Matrix identity(int n) {
		Matrix identity = new Matrix(n, n);
		for (int i = 0; i < n; i ++)
			identity.setElement(i, i, 1);
		return identity;
	}
	
	
	public int getHeight() {
		return this.values.length;
	}
	
	public int getWidth() {
		if (this.getHeight() > 0)
			return this.values[0].length;
		else
			return 0;
	}
	
	public double getElement(int i, int j) {
		return this.values[i][j];
	}
	
	public void setElement(int i, int j, double val) {
		this.values[i][j] = val;
	}
	
	protected double[][] getArray() {
		return values;
	}
	
	public Matrix T() {
		return this.transpose();
	}
	
	public Matrix transpose() {
		Matrix thisT = new Matrix(this.getWidth(), this.getHeight());
		for (int i = 0; i < this.getWidth(); i ++)
			for (int j = 0; j < this.getHeight(); j ++)
				thisT.setElement(i, j, this.getElement(j, i));
		return thisT;
	}
	
	public Matrix inverse() {
		if (this.getWidth() != this.getHeight())
			throw new IllegalArgumentException("Only square matrices, have inverses. "+this+" is not square.");
		
		double det = this.determinant();
		if (det == 0)
			System.err.println(this+" is singular and thus has an infinite determinant.");
		
		Matrix inv = new Matrix(this.getHeight(), this.getWidth());
		for (int i = 0; i < this.getHeight(); i ++) {
			for (int j = 0; j < this.getWidth(); j ++) {
				inv.setElement(i, j, this.cofactor(j, i)/det);
			}
		}
		return inv;
	}
	
	public double determinant() {
		if (this.getWidth() != this.getHeight())
			throw new IllegalArgumentException("Only square matrices have determinants. "+this+" is not square.");
		if (this.getHeight() == 0)
			return 1;
		if (this.getHeight() == 1)
			return this.getElement(0, 0);
		if (this.getHeight() == 2)
			return this.getElement(0, 0) * this.getElement(1, 1) -
					this.getElement(0, 1) * this.getElement(1, 0);
		
		double det = 0;
		for (int j = 0; j < this.getWidth(); j ++) {
			if (j%2 == 0)
				det += this.getElement(0, j)*this.submatrix(0, j).determinant();
			else
				det -= this.getElement(0, j)*this.submatrix(0, j).determinant();
		}
		return det;
	}
	
	public Matrix submatrix(double I, double J) {
		Matrix sub = new Matrix(this.getHeight()-1, this.getWidth()-1);
		for (int i = 0; i < this.getHeight()-1; i ++)
			for (int j = 0; j < this.getWidth()-1; j ++)
				sub.setElement(i, j, this.getElement(i<I ? i : i+1, j<J ? j : j+1));
		return sub;
	}
	
	public double cofactor(double i, double j) {
		if ((i + j)%2 == 0)
			return this.submatrix(i, j).determinant();
		else
			return -this.submatrix(i, j).determinant();
	}
	
	public Matrix plus(Matrix that) {
		if (this.getHeight() != that.getHeight() || this.getWidth() != that.getWidth())
			throw new IllegalArgumentException("Matrix dimensions must match. Cannot multiply\n"+this+" by\n"+that);
		Matrix sum = new Matrix(this.getHeight(), this.getWidth());
		for (int i = 0; i < this.getHeight(); i ++)
			for (int j = 0; j < this.getWidth(); j ++)
				sum.setElement(i, j, this.getElement(i, j) + that.getElement(i, j));
		return sum;
	}
	
	public Matrix minus(Matrix that) {
		return this.plus(that.times(-1));
	}
	
	public Matrix times(double c) {
		Matrix out = new Matrix(this.getHeight(), this.getWidth());
		for (int i = 0; i < this.getHeight(); i ++)
			for (int j = 0; j < this.getWidth(); j ++)
				out.setElement(i, j, this.getElement(i, j)*c);
		return out;
	}
	
	public Matrix times(Matrix that) {
		if (this.getWidth() != that.getHeight())
			throw new IllegalArgumentException("Matrix dimensions must match. Cannot multiply\n"+this+" by\n"+that);
		Matrix product = new Matrix(this.getHeight(), that.getWidth());
		for (int i = 0; i < this.getHeight(); i ++)
			for (int j = 0; j < that.getWidth(); j ++)
				for (int k = 0; k < this.getWidth(); k ++)
					product.setElement(i, j, product.getElement(i, j)
							+ this.getElement(i, k) * that.getElement(k, j));
		return product;
	}
	
	public String toString() {
		StringBuilder str = new StringBuilder("[ ");
		for (int i = 0; i < this.getHeight(); i ++) {
			for (int j = 0; j < this.getWidth(); j ++)
				str.append(this.getElement(i, j)).append(", ");
			str = new StringBuilder(str.substring(0, str.length() - 2) + ";\n  ");
		}
		return str.substring(0, str.length()-4) + " ]";
	}
	
	
	public static void main(String[] args) {
		Matrix a = new Matrix(new double[][] {{1,4,7},{3,0,5},{-1,9,11}});
		System.out.println(a); //TODO: delete this later
		System.out.println(a.determinant());
		System.out.println(a.submatrix(1, 2));
		System.out.println(a.cofactor(1, 2));
		System.out.println(a.inverse());
		System.out.println(a.inverse().transpose());
		System.out.println(a.transpose().times(a.inverse().transpose()));
	}
	
}
