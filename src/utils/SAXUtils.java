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
package utils;

import org.xml.sax.Attributes;
import org.xml.sax.helpers.AttributesImpl;

/**
 * A class for working with and manipulating XML stuff.
 * 
 * @author jkunimune
 */
public class SAXUtils {
	
	public static String encode(String s0) { //encode with the ampersand notation
		String s1 = "";
		for (int i = 0; i < s0.length(); i ++) {
			if (s0.charAt(i) >= 128 || s0.charAt(i) == '&')
				s1 += "&#" + (int)s0.charAt(i) + ";";
			else
				s1 += s0.charAt(i);
		}
		return s1;
	}
	
	
	public static Attributes removeAttribute(Attributes attrs, String... ss) {
		AttributesImpl modAttrs = new AttributesImpl(attrs);
		for (String s: ss)
			modAttrs.removeAttribute(modAttrs.getIndex(s));
		return modAttrs;
	}
}
