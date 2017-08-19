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
public class Transverse extends Projection {
	
	public static final double[] TRANSVERSE_AXIS = {0, 0, 0};
	
	private final Projection base;
	
	
	
	public Transverse(Projection base) {
		super(base);
		this.base = base;
	}
	
	
	
	@Override
	public double[] project(double lat, double lon) {
		return base.project(obliquifySphc(lat, lon, TRANSVERSE_AXIS));
	}
	
	
	@Override
	public double[] inverse(double x, double y) {
		return obliquifyPlnr(base.inverse(x, y), TRANSVERSE_AXIS);
	}
	
	
	@Override
	public void setParameters(double... params) {
		base.setParameters(params);
		this.aspectRatio = base.aspectRatio;
	}
	
}
