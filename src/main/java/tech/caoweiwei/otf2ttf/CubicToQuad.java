package tech.caoweiwei.otf2ttf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

import org.apache.commons.lang3.ArrayUtils;

public class CubicToQuad {

	public Point[] calcPowerCoefficients(Point p1, Point c1, Point c2, Point p2) {
		// point(t) = p1*(1-t)^3 + c1*t*(1-t)^2 + c2*t^2*(1-t) + p2*t^3 = a*t^3 + b*t^2 + c*t + d
		// for each t value, so
		// a = (p2 - p1) + 3 * (c1 - c2)
		// b = 3 * (p1 + c2) - 6 * c1
		// c = 3 * (c1 - p1)
		// d = p1
		Point a = p2.sub(p1).add(c1.sub(c2).mul(3));
		Point b = p1.add(c2).mul(3).sub(c1.mul(6));
		Point c = c1.sub(p1).mul(3);
		Point d = p1;
		return new Point[] { a, b, c, d };
	}

	public Point calcPoint(Point a, Point b, Point c, Point d, double t) {
		// a*t^3 + b*t^2 + c*t + d = ((a*t + b)*t + c)*t + d
		return a.mul(t).add(b).mul(t).add(c).mul(t).add(d);
	}

	public Point calcPointQuad(Point a, Point b, Point c, double t) {
		// a*t^2 + b*t + c = (a*t + b)*t + c
		return a.mul(t).add(b).mul(t).add(c);
	}

	public Point calcPointDerivative(Point a, Point b, Point c, Point d, double t) {
		// d/dt[a*t^3 + b*t^2 + c*t + d] = 3*a*t^2 + 2*b*t + c = (3*a*t + 2*b)*t + c
		return a.mul(3 * t).add(b.mul(2)).mul(t).add(c);
	}

	public double[] quadSolve(double a, double b, double c) {
		// a*x^2 + b*x + c = 0
		if (a == 0) {
			return (b == 0) ? new double[] {} : new double[] { -c / b };
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

	public double cubicRoot(double x) {
		return (x < 0) ? -Math.pow(-x, 1 / 3) : Math.pow(x, 1 / 3);
	}

	public double[] cubicSolve(double a, double b, double c, double d) {
		// a*x^3 + b*x^2 + c*x + d = 0
		if (a == 0) {
			return quadSolve(b, c, d);
		}
		// solve using Cardan's method, which is described in paper of R.W.D. Nickals
		// http://www.nickalls.org/dick/papers/maths/cubic1993.pdf (doi:10.2307/3619777)
		double xn = -b / (3 * a); // point of symmetry x coordinate
		double yn = ((a * xn + b) * xn + c) * xn + d; // point of symmetry y coordinate
		double deltaSq = (b * b - 3 * a * c) / (9 * a * a); // delta^2
		double hSq = 4 * a * a * Math.pow(deltaSq, 3); // h^2
		double D3 = yn * yn - hSq;
		if (D3 > 0) { // 1 real root
			double D3Sqrt = Math.sqrt(D3);
			return new double[] { xn + cubicRoot((-yn + D3Sqrt) / (2 * a)) + cubicRoot((-yn - D3Sqrt) / (2 * a)) };
		} else if (D3 == 0) { // 2 real roots
			double delta1 = cubicRoot(yn / (2 * a));
			return new double[] { xn - 2 * delta1, xn + delta1 };
		}
		// 3 real roots
		double theta = Math.acos(-yn / Math.sqrt(hSq)) / 3;
		double delta = Math.sqrt(deltaSq);
		return new double[] { xn + 2 * delta * Math.cos(theta), xn + 2 * delta * Math.cos(theta + Math.PI * 2 / 3),
				xn + 2 * delta * Math.cos(theta + Math.PI * 4 / 3) };
	}

	public double minDistanceToQuad(Point point, Point p1, Point c1, Point p2) {
		// f(t) = (1-t)^2 * p1 + 2*t*(1 - t) * c1 + t^2 * p2 = a*t^2 + b*t + c, t in [0, 1],
		// a = p1 + p2 - 2 * c1
		// b = 2 * (c1 - p1)
		// c = p1; a, b, c are vectors because p1, c1, p2 are vectors too
		// The distance between given point and quadratic curve is equal to
		// sqrt((f(t) - point)^2), so these expression has zero derivative by t at points where
		// (f'(t), (f(t) - point)) = 0.
		// Substituting quadratic curve as f(t) one could obtain a cubic equation
		// e3*t^3 + e2*t^2 + e1*t + e0 = 0 with following coefficients:
		// e3 = 2 * a^2
		// e2 = 3 * a*b
		// e1 = (b^2 + 2 * a*(c - point))
		// e0 = (c - point)*b
		// One of the roots of the equation from [0, 1], or t = 0 or t = 1 is a value of t
		// at which the distance between given point and quadratic Bezier curve has minimum.
		// So to find the minimal distance one have to just pick the minimum value of
		// the distance on set {t = 0 | t = 1 | t is root of the equation from [0, 1] }.

		Point a = p1.add(p2).sub(c1.mul(2));
		Point b = c1.sub(p1).mul(2);
		Point c = p1;
		double e3 = 2 * a.sqr();
		double e2 = 3 * a.dot(b);
		double e1 = (b.sqr() + 2 * a.dot(c.sub(point)));
		double e0 = c.sub(point).dot(b);
		Double[] candidates = ArrayUtils.addAll(Stream.of(ArrayUtils.toObject(cubicSolve(e3, e2, e1, e0)))
				.filter(t -> t > 0 && t < 1).toArray(Double[]::new), new Double[] { new Double(0L), new Double(1L) });

		double minDistance = 1e9;
		for (int i = 0; i < candidates.length; i++) {
			double distance = calcPointQuad(a, b, c, candidates[i]).sub(point).dist();
			if (distance < minDistance) {
				minDistance = distance;
			}
		}
		return minDistance;
	}

	public Point[] processSegment(Point a, Point b, Point c, Point d, double t1, double t2) {
		// Find a single control point for given segment of cubic Bezier curve
		// These control point is an interception of tangent lines to the boundary points
		// Let's denote that f(t) is a vector function of parameter t that defines the cubic Bezier
		// curve,
		// f(t1) + f'(t1)*z1 is a parametric equation of tangent line to f(t1) with parameter z1
		// f(t2) + f'(t2)*z2 is the same for point f(t2) and the vector equation
		// f(t1) + f'(t1)*z1 = f(t2) + f'(t2)*z2 defines the values of parameters z1 and z2.
		// Defining fx(t) and fy(t) as the x and y components of vector function f(t) respectively
		// and solving the given system for z1 one could obtain that
		//
		// -(fx(t2) - fx(t1))*fy'(t2) + (fy(t2) - fy(t1))*fx'(t2)
		// z1 = ------------------------------------------------------.
		// -fx'(t1)*fy'(t2) + fx'(t2)*fy'(t1)
		//
		// Let's assign letter D to the denominator and note that if D = 0 it means that the curve
		// actually
		// is a line. Substituting z1 to the equation of tangent line to the point f(t1), one could
		// obtain that
		// cx = [fx'(t1)*(fy(t2)*fx'(t2) - fx(t2)*fy'(t2)) + fx'(t2)*(fx(t1)*fy'(t1) -
		// fy(t1)*fx'(t1))]/D
		// cy = [fy'(t1)*(fy(t2)*fx'(t2) - fx(t2)*fy'(t2)) + fy'(t2)*(fx(t1)*fy'(t1) -
		// fy(t1)*fx'(t1))]/D
		// where c = (cx, cy) is the control point of quadratic Bezier curve.

		Point f1 = calcPoint(a, b, c, d, t1);
		Point f2 = calcPoint(a, b, c, d, t2);
		Point f1_ = calcPointDerivative(a, b, c, d, t1);
		Point f2_ = calcPointDerivative(a, b, c, d, t2);

		double D = -f1_.x * f2_.y + f2_.x * f1_.y;
		if (Math.abs(D) < 1e-8) {
			return new Point[] { f1, f1.add(f2).div(2), f2 }; // straight line segment
		}
		double cx = (f1_.x * (f2.y * f2_.x - f2.x * f2_.y) + f2_.x * (f1.x * f1_.y - f1.y * f1_.x)) / D;
		double cy = (f1_.y * (f2.y * f2_.x - f2.x * f2_.y) + f2_.y * (f1.x * f1_.y - f1.y * f1_.x)) / D;
		return new Point[] { f1, new Point(cx, cy), f2 };
	}

	public boolean isSegmentApproximationClose(Point a, Point b, Point c, Point d, double tmin, double tmax, Point p1,
			Point c1, Point p2, double errorBound) {
		// a,b,c,d define cubic curve
		// tmin, tmax are boundary points on cubic curve
		// p1, c1, p2 define quadratic curve
		// errorBound is maximum allowed distance
		// Try to find maximum distance between one of N points segment of given cubic
		// and corresponding quadratic curve that estimates the cubic one, assuming
		// that the boundary points of cubic and quadratic points are equal.
		//
		// The distance calculation method comes from Hausdorff distance defenition
		// (https://en.wikipedia.org/wiki/Hausdorff_distance), but with following simplifications
		// * it looks for maximum distance only for finite number of points of cubic curve
		// * it doesn't perform reverse check that means selecting set of fixed points on
		// the quadratic curve and looking for the closest points on the cubic curve
		// But this method allows easy estimation of approximation error, so it is enough
		// for practical purposes.

		double n = 10; // number of points + 1
		double dt = (tmax - tmin) / n;
		for (double t = tmin + dt; t < tmax - dt; t += dt) {
			// don't check distance on boundary points because they should be the same
			Point point = calcPoint(a, b, c, d, t);
			if (minDistanceToQuad(point, p1, c1, p2) > errorBound) {
				return false;
			}
		}
		return true;
	}

	public boolean _isApproximationClose(Point a, Point b, Point c, Point d, Point[][] quadCurves, double errorBound) {
		double dt = 1.0 / quadCurves.length;
		for (int i = 0; i < quadCurves.length; i++) {
			Point p1 = quadCurves[i][0];
			Point c1 = quadCurves[i][1];
			Point p2 = quadCurves[i][2];
			if (!isSegmentApproximationClose(a, b, c, d, i * dt, (i + 1) * dt, p1, c1, p2, errorBound)) {
				return false;
			}
		}
		return true;
	}

	public Point[][] fromFlatArray(double[] points) {
		Point[][] result = new Point[][] {};
		double segmentsNumber = (points.length - 2) / 4;
		for (int i = 0; i < segmentsNumber; i++) {
			ArrayUtils.addAll(result,
					new Point[][] { new Point[] { new Point(points[4 * i], points[4 * i + 1]),
							new Point(points[4 * i + 2], points[4 * i + 3]),
							new Point(points[4 * i + 4], points[4 * i + 5]) } });
		}
		return result;
	}

	public double[] toFlatArray(Point[][] quadsList) {
		List<Double> result = new ArrayList<>();
		result.add(quadsList[0][0].x);
		result.add(quadsList[0][0].y);
		for (int i = 0; i < quadsList.length; i++) {
			result.add(quadsList[i][1].x);
			result.add(quadsList[i][1].y);
			result.add(quadsList[i][2].x);
			result.add(quadsList[i][2].y);
		}
		return ArrayUtils.toPrimitive(result.toArray(new Double[result.size()]));
	}

	public boolean isApproximationClose(double p1x, double p1y, double c1x, double c1y, double c2x, double c2y,
			double p2x, double p2y, double[] quads, double errorBound) {
		// TODO: rewrite it in C-style and remove _isApproximationClose
		Point[] pc = calcPowerCoefficients(new Point(p1x, p1y), new Point(c1x, c1y), new Point(c2x, c2y),
				new Point(p2x, p2y));
		return _isApproximationClose(pc[0], pc[1], pc[2], pc[3], fromFlatArray(quads), errorBound);
	}

	/*
	 * Split cubic b¨¦zier curve into two cubic curves, see details here:
	 * https://math.stackexchange.com/questions/877725
	 */
	public double[][] subdivideCubic(double x1, double y1, double x2, double y2, double x3, double y3, double x4,
			double y4, double t) {
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

	public double byNumber(double x, double y) {
		return x - y;
	}

	/*
	 * Find inflection points on a cubic curve, algorithm is similar to this one:
	 * http://www.caffeineowl.com/graphics/2d/vectorial/cubic-inflexion.html
	 */
	public double[] solveInflections(double x1, double y1, double x2, double y2, double x3, double y3, double x4,
			double y4) {
		double p = -(x4 * (y1 - 2 * y2 + y3)) + x3 * (2 * y1 - 3 * y2 + y4) + x1 * (y2 - 2 * y3 + y4)
				- x2 * (y1 - 3 * y3 + 2 * y4);
		double q = x4 * (y1 - y2) + 3 * x3 * (-y1 + y2) + x2 * (2 * y1 - 3 * y3 + y4) - x1 * (2 * y2 - 3 * y3 + y4);
		double r = x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3);

		Double[] temp = Stream.of(ArrayUtils.toObject(quadSolve(p, q, r))).filter(t -> t > 1e-8 && t < (1 - 1e-8))
				.toArray(Double[]::new);
		Arrays.sort(temp, new Comparator<Double>() {
			@Override
			public int compare(Double o1, Double o2) {
				return (int) byNumber(o1, o2);
			}
		});
		return ArrayUtils.toPrimitive(temp);
	}

	/*
	 * Approximate cubic Bezier curve defined with base points p1, p2 and control points c1, c2 with
	 * with a few quadratic Bezier curves. The function uses tangent method to find quadratic
	 * approximation of cubic curve segment and simplified Hausdorff distance to determine number of
	 * segments that is enough to make error small. In general the method is the same as described
	 * here: https://fontforge.github.io/bezier.html.
	 */
	public double[] _cubicToQuad(double p1x, double p1y, double c1x, double c1y, double c2x, double c2y, double p2x,
			double p2y, double errorBound) {
		Point p1 = new Point(p1x, p1y);
		Point c1 = new Point(c1x, c1y);
		Point c2 = new Point(c2x, c2y);
		Point p2 = new Point(p2x, p2y);
		Point[] pc = calcPowerCoefficients(p1, c1, c2, p2);
		Point a = pc[0], b = pc[1], c = pc[2], d = pc[3];

		Point[][] approximation = new Point[][] {};
		for (int segmentsCount = 1; segmentsCount <= 8; segmentsCount++) {
			approximation = new Point[][] {};
			for (double t = 0.0; t < 1; t += 1.0 / segmentsCount) {
				approximation = ArrayUtils.add(approximation, processSegment(a, b, c, d, t, t + 1 / segmentsCount));
			}
			if (segmentsCount == 1 && (approximation[0][1].sub(p1).dot(c1.sub(p1)) < 0
					|| approximation[0][1].sub(p2).dot(c2.sub(p2)) < 0)) {
				// approximation concave, while the curve is convex (or vice versa)
				continue;
			}
			if (_isApproximationClose(a, b, c, d, approximation, errorBound)) {
				break;
			}
		}
		return toFlatArray(approximation);
	}

	/*
	 * If this curve has any inflection points, split the curve and call _cubicToQuad function on
	 * each resulting curve.
	 */
	public double[] cubicToQuad(double p1x, double p1y, double c1x, double c1y, double c2x, double c2y, double p2x,
			double p2y, double errorBound) {
		double[] inflections = solveInflections(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y);

		if (inflections.length == 0) {
			return _cubicToQuad(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y, errorBound);
		}

		double[] result = new double[] {};
		double[] curve = new double[] { p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y };
		double prevPoint = 0;
		double[] quad;
		double[][] split;

		for (int inflectionIdx = 0; inflectionIdx < inflections.length; inflectionIdx++) {
			split = subdivideCubic(curve[0], curve[1], curve[2], curve[3], curve[4], curve[5], curve[6], curve[7],
					// we make a new curve, so adjust inflection point accordingly
					1.0 - (1.0 - inflections[inflectionIdx]) / (1.0 - prevPoint));

			quad = _cubicToQuad(split[0][0], split[0][1], split[0][2], split[0][3], split[0][4], split[0][5],
					split[0][6], split[0][7], errorBound);

			result = ArrayUtils.addAll(result, ArrayUtils.subarray(quad, 0, quad.length - 2));
			curve = split[1];
			prevPoint = inflections[inflectionIdx];
		}

		quad = _cubicToQuad(curve[0], curve[1], curve[2], curve[3], curve[4], curve[5], curve[6], curve[7], errorBound);

		return ArrayUtils.addAll(result, quad);
	}

}
