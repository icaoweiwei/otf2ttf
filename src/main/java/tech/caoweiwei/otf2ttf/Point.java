package tech.caoweiwei.otf2ttf;

public class Point {
	public double x;
	public double y;

	public Point() {
	}

	public Point(double x, double y) {
		this.x = x;
		this.y = y;
	}

	@Override
	public String toString() {
		return "Point [x=" + x + ", y=" + y + "]";
	}

	public Point add(Point point) {
		return new Point(this.x + point.x, this.y + point.y);
	}

	public Point sub(Point point) {
		return new Point(this.x - point.x, this.y - point.y);
	}

	public Point mul(double value) {
		return new Point(this.x * value, this.y * value);
	}

	public Point div(double value) {
		return new Point(this.x / value, this.y / value);
	}

	public double dist() {
		return Math.sqrt(this.x * this.x + this.y * this.y);
	}

	public double sqr() {
		return this.x * this.x + this.y * this.y;
	}

	public double dot(Point point) {
		return this.x * point.x + this.y * point.y;
	}
}
