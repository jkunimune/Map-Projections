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
package maps;

/**
 * A projection that uses the same equations as another, but shifts the aspect
 * 
 * @author jkunimune
 */
public class Oblique extends Projection {
	
	private final Projection base;
	private final double[] axis;
	
	
	
	public Oblique(Projection base, String name, double... axis) {
		super(name, base);
		this.base = base;
		this.axis = axis;
	}
	
	
	
	@Override
	public double[] project(double lat, double lon) {
		return base.project(transformFromOblique(lat, lon, axis));
	}
	
	
	@Override
	public double[] inverse(double x, double y) {
		double[] coords = base.inverse(x, y);
		if (coords == null) 	return null;
		return transformToOblique(coords, axis);
	}
	
	
	@Override
	public void initialize(double... params) {
		base.initialize(params);
		this.bounds = base.bounds;
	}
	
}
