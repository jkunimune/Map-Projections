package utils;

public class BoundingBox {
	public final double xMin;
	public final double xMax;
	public final double yMin;
	public final double yMax;
	public final double width;
	public final double height;

	public BoundingBox(double width, double height) {
		this(-width/2, width/2, -height/2, height/2);
	}

	public BoundingBox(double xMin, double xMax, double yMin, double yMax) {
		this.xMin = xMin;
		this.xMax = xMax;
		this.yMin = yMin;
		this.yMax = yMax;
		this.width = xMax - xMin;
		this.height = yMax - yMin;
	}

}