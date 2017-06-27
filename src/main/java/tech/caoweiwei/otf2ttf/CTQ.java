package tech.caoweiwei.otf2ttf;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.Stream;

import org.apache.commons.lang3.ArrayUtils;

import com.alibaba.fastjson.JSONObject;

public class CTQ {
	public double[] cubicToQuad(double p1x, double p1y, double c1x, double c1y, double c2x, double c2y, double p2x,
			double p2y, double errorBound) {
		return new CubicToQuad().cubicToQuad(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y, errorBound);
	}

	public Contour[] removeMids(Contour[] contour, double err) {
		int last = contour.length - 1;
		for (int j = 0; j < contour.length - 1; j++) {
			if (Math.abs(contour[j].x - contour[j + 1].x) < 1 && Math.abs(contour[j].y - contour[j + 1].y) < 1) {
				contour[j + 1].rem = true;
				contour[j].on = true;
			}
		}
		while (last > 0 && Math.abs(contour[0].x - contour[last].x) < 1
				&& Math.abs(contour[0].y - contour[last].y) < 1) {
			contour[last].rem = true;
			contour[0].on = true;
			last -= 1;
		}
		contour = Stream.of(contour).filter(x -> !x.rem).toArray(Contour[]::new);

		last = contour.length - 1;
		for (int j = 1; j < contour.length - 1; j++) {
			if (!contour[j - 1].on && contour[j].on && !contour[j + 1].on) {
				double mx = contour[j - 1].x + contour[j + 1].x;
				double my = contour[j - 1].y + contour[j + 1].y;
				double dy = contour[j - 1].y - contour[j + 1].y;
				if (Math.abs(dy) >= 1 && Math.abs(contour[j].x * 2 - mx) < err
						&& Math.abs(contour[j].y * 2 - my) < err) {
					contour[j].rem = true;
				}
			}
		}
		if (!contour[last].rem && !contour[last].on && contour[0].on && !contour[1].on) {
			double mx = contour[last].x + contour[1].x;
			double my = contour[last].y + contour[1].y;
			if (Math.abs(contour[0].x * 2 - mx) < err && Math.abs(contour[0].y * 2 - my) < err) {
				contour[0].rem = true;
			}
		}
		return Stream.of(contour).filter(x -> !x.rem).toArray(Contour[]::new);
	}

	public Contour[] canonicalStart(Contour[] points) {
		int jm = 0;
		for (int j = 1; j < points.length; j++) {
			if (points[j].x < points[jm].x || (points[j].x == points[jm].x && points[j].y < points[jm].y)) {
				jm = j;
			}
		}
		Contour[] result = ArrayUtils.addAll(Arrays.copyOfRange(points, jm, points.length),
				Arrays.copyOfRange(points, 0, jm));
		ArrayUtils.reverse(result);
		return result;
	}

	public double[] quadSolve(double a, double b, double c) {
		// a*x^2 + b*x + c = 0
		if (a == 0) {
			return b == 0 ? new double[] {} : new double[] { -c / b };
		}
		double D = b * b - 4 * a * c;
		if (D < 0) {
			return new double[] {};
		} else if (D == 0) {
			return new double[] { -b / (2 * a) };
		}
		double DSqrt = Math.sqrt(D);
		return new double[] { (-b - DSqrt) / (2 * a), (-b + DSqrt) / (2 * a) };
	}

	public double[][] splitAt(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4,
			double t) {
		double u = 1 - t, v = t;

		double bx = x1 * u + x2 * v;
		double sx = x2 * u + x3 * v;
		double fx = x3 * u + x4 * v;
		double cx = bx * u + sx * v;
		double ex = sx * u + fx * v;
		double dx = cx * u + ex * v;

		double by = y1 * u + y2 * v;
		double sy = y2 * u + y3 * v;
		double fy = y3 * u + y4 * v;
		double cy = by * u + sy * v;
		double ey = sy * u + fy * v;
		double dy = cy * u + ey * v;

		return new double[][] { new double[] { x1, y1, bx, by, cx, cy, dx, dy },
				new double[] { dx, dy, ex, ey, fx, fy, x4, y4 } };
	}

	public double[][] splitAtTs(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4,
			double[] ts) {
		if (ts.length == 0)
			return new double[][] { new double[] { x1, y1, x2, y2, x3, y3, x4, y4 } };
		double[][] s = splitAt(x1, y1, x2, y2, x3, y3, x4, y4, ts[0]);
		if (ts.length == 1) {
			return s;
		} else {
			return ArrayUtils.addAll(new double[][] { s[0] }, splitAtTs(s[1][0], s[1][1], s[1][2], s[1][3], s[1][4],
					s[1][5], s[1][6], s[1][7], ArrayUtils.subarray(ts, 1, ts.length)));
		}
	}

	public boolean fin(double t) {
		return t > 0.001 && t < 0.999;
	}

	public double substraction(double x, double y) {
		return x - y;
	}

	public double[][] getSplitAtXY(double x1, double y1, double x2, double y2, double x3, double y3, double x4,
			double y4, Double splitAtX, Double splitAtY) {
		double ax = 3 * (-x1 + 3 * x2 - 3 * x3 + x4);
		double bx = 6 * (x1 - 2 * x2 + x3);
		double cx = 3 * (x2 - x1);
		double ay = 3 * (-y1 + 3 * y2 - 3 * y3 + y4);
		double by = 6 * (y1 - 2 * y2 + y3);
		double cy = 3 * (y2 - y1);
		double[] ts = new double[] {};
		if (splitAtX != null) {
			ts = ArrayUtils.addAll(ts, quadSolve(ax, bx, cx));
		}
		if (splitAtY != null) {
			ts = ArrayUtils.addAll(ts, quadSolve(ay, by, cy));
		}
		Double[] temp = Stream.of(ArrayUtils.toObject(ts)).filter(t -> (t > 0.001 && t < 0.999)).toArray(Double[]::new);
		Arrays.sort(temp, new Comparator<Double>() {
			@Override
			public int compare(Double o1, Double o2) {
				return (int) substraction(o1, o2);
			}
		});
		return splitAtTs(x1, y1, x2, y2, x3, y3, x4, y4, ArrayUtils.toPrimitive(temp));
	}

	public double[] handle(Contour z1, Contour z2, Contour z3, Contour z4, Double splitAtX, Double splitAtY,
			double err) {
		double[][] segs = getSplitAtXY(z1.x, z1.y, z2.x, z2.y, z3.x, z3.y, z4.x, z4.y, splitAtX, splitAtY);
		double[] ss = new double[] {};
		for (double[] s : segs) {
			double[] a = cubicToQuad(s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], err);
			for (int j = ss.length != 0 ? 2 : 0; j < a.length; j++) {
				ss = ArrayUtils.add(ss, a[j]);
			}
		}
		return ss;
	}

	public Contour[] toquad(Contour[] contour, Double splitAtX, Double splitAtY, Double err) {
		if (contour.length == 0)
			return new Contour[] {};
		if (contour.length == 1)
			return new Contour[] { contour[0] };
		Contour[] newcontour = new Contour[] {};
		contour = ArrayUtils.add(contour, new Contour(contour[0].x, contour[0].y, true));
		for (int j = 0; j < contour.length; j++) {
			if (contour[j].on) {
				newcontour = ArrayUtils.add(newcontour, contour[j]);
			} else {
				Contour z1 = newcontour[newcontour.length - 1];
				Contour z2 = contour[j];
				Contour z3 = contour[j + 1];
				Contour z4 = contour[j + 2];
				double[] quadzs = handle(z1, z2, z3, z4, splitAtX, splitAtY, err);
				boolean on = false;

				double mx = (z1.x + z4.x) / 2;
				double my = (z1.y + z4.y) / 2;
				double bw = Math.abs(z4.x - z1.x);
				double bh = Math.abs(z4.y - z1.y);
				for (int k = 2; k < quadzs.length - 2; k += 2) {
					double cx = quadzs[k];
					double cy = quadzs[k + 1];
					newcontour = ArrayUtils.add(newcontour, new Contour(cx, cy, on));
					on = !on;
				}
				newcontour = ArrayUtils.add(newcontour, new Contour(z4.x, z4.y, true));
				j += 2;
			}
		}
		return canonicalStart(removeMids(newcontour, err));
	}

	public boolean haspt(Contour[] c) {
		return c != null && c.length > 1;
	}

	public Contour[][] c2qContours(Contour[][] contours, Double splitAtX, Double splitAtY, Double err) {
		Contour[][] ans = new Contour[][] {};
		for (Contour[] c : contours) {
			Contour[] c1 = toquad(c, splitAtX, splitAtY, err == null ? 1 : err);
			if (haspt(c1)) {
				ans = ArrayUtils.add(ans, c1);
			}
		}
		return ans;
	}

	public void exports(JSONObject font, Double splitAtX, Double splitAtY, Double err) {
		font.put("CFF_", null);
		JSONObject glyf = font.getJSONObject("glyf");
		for (String k : glyf.keySet()) {
			JSONObject g = glyf.getJSONObject(k);
			Contour[][] contours = g.getObject("contours", Contour[][].class);
			if (contours != null) {
				contours = c2qContours(contours, splitAtX, splitAtY, err);
			}
			g.put("contours", contours);
			g.put("stemH", null);
			g.put("stemV", null);
			g.put("hintMasks", null);
			g.put("contourMasks", null);
		}
		JSONObject maxp = font.getJSONObject("maxp");
		maxp.put("version", 1.0);
	}
}
