package utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.Double.parseDouble;
import static java.lang.String.format;

public class Quantity {
	public final double value;
	public final String units;
	
	public Quantity(double value, String units) {
		this.value = value;
		this.units = units;
	}
	
	public static Quantity parseQuantity(String s) {
		Pattern format = Pattern.compile("([0-9.e-]+)\\s*([a-z]*)", Pattern.CASE_INSENSITIVE);
		Matcher match = format.matcher(s);
		if (match.matches())
			return new Quantity(parseDouble(match.group(1)), match.group(2));
		else
			throw new IllegalArgumentException(format(
					"I could not interpret '%s' as a number or as a number with a unit.", s));
	}
}
