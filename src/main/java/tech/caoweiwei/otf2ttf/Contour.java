package tech.caoweiwei.otf2ttf;

import com.alibaba.fastjson.annotation.JSONField;

public class Contour {
	public double x;
	public double y;
	public boolean on;
	@JSONField(serialize = false)
	public boolean rem;

	public Contour() {
	}

	public Contour(double x, double y, boolean on) {
		this.x = x;
		this.y = y;
		this.on = on;
	}

	@Override
	public String toString() {
		return "Contour [x=" + x + ", y=" + y + ", on=" + on + ", rem=" + rem + "]";
	}
}
