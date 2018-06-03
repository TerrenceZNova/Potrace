package potrace;

class Math {

    /* integer arithmetic */
    static int sign(int x) {
        return Integer.compare(x, 0);
    }

    static int sign(double x) {
        return ((x) > 0 ? 1 : (x) < 0 ? -1 : 0);
    }

    static int abs(int a) {
        return ((a) > 0 ? (a) : -(a));
    }

    static int min(int a, int b) {
        return ((a) < (b) ? (a) : (b));
    }

    static int max(int a, int b) {
        return ((a) > (b) ? (a) : (b));
    }

    static int sq(int a) {
        return ((a) * (a));
    }

    static int cu(int a) {
        return ((a) * (a) * (a));
    }

    static int mod(int a, int n) {
        return a >= n ? a % n : a >= 0 ? a : n - 1 - (-1 - a) % n;
    }


    /* calculate p1 x p2 */
    static int xprod(IntPoint p1, IntPoint p2) {
        return p1.X * p2.Y - p1.Y * p2.X;
    }

    /* calculate p1 x p2 */
    static double xprod(DoublePoint p1, DoublePoint p2) {
        return p1.X * p2.Y - p1.Y * p2.X;
    }

    /* return 1 if a <= b < c < a, in a cyclic sense (mod n) */
    static boolean cyclic(int a, int b, int c) {
        if (a <= c) {
            return (a <= b && b < c);
        } else {
            return (a <= b || b < c);
        }
    }

    static int floordiv(int a, int n) {
        return a >= 0 ? a / n : -1 - (-1 - a) / n;
    }

    static void pointslope(Path pp, int i, int j, DoublePoint ctr, DoublePoint dir) {
        assert (i < j);

        int n = pp.pt.length;
        SumStruct[] sums = pp.Sums;

        double x, y, x2, xy, y2;
        double k;
        double a, b, c, lambda2, l;
        int r = 0; /* rotations from i to j */

        while (j >= n) {
            j -= n;
            r += 1;
        }
        while (i >= n) {
            i -= n;
            r -= 1;
        }
        while (j < 0) {
            j += n;
            r -= 1;
        }
        while (i < 0) {
            i += n;
            r += 1;
        }

        x = sums[j + 1].X - sums[i].X + r * sums[n].X;
        y = sums[j + 1].Y - sums[i].Y + r * sums[n].Y;
        x2 = sums[j + 1].X2 - sums[i].X2 + r * sums[n].X2;
        xy = sums[j + 1].XY - sums[i].XY + r * sums[n].XY;
        y2 = sums[j + 1].Y2 - sums[i].Y2 + r * sums[n].Y2;
        k = j + 1 - i + r * n;

        //ctr = new DoublePoint(x / k, y / k);
        ctr.X = x / k;
        ctr.Y = y / k;

        a = (x2 - (double) x * x / k) / k;
        b = (xy - (double) x * y / k) / k;
        c = (y2 - (double) y * y / k) / k;

        lambda2 = (a + c + java.lang.Math.sqrt((a - c) * (a - c) + 4 * b * b)) / 2; // larger e.value

        /* now find e.vector for lambda2 */
        a -= lambda2;
        c -= lambda2;

        if (java.lang.Math.abs(a) >= java.lang.Math.abs(c)) {
            l = java.lang.Math.sqrt(a * a + b * b);
            if (l != 0) {
                //dir = new DoublePoint(-b / l, a / l);
                dir.X = -b / l;
                dir.Y = a / l;
            }
        } else {
            l = java.lang.Math.sqrt(c * c + b * b);
            if (l != 0) {
                //dir = new DoublePoint(-c / l, b / l);
                dir.X = -c / l;
                dir.Y = b / l;
            }
        }
        if (l == 0) {
            // Sometimes this can happen when k=4: the two eigenvalues coincide.

            //dir = new DoublePoint();
            dir = null;     //这里这么改可能会出问题
        }
    }

    /* Apply quadratic form Q to vector w = (w.x,w.y) */
    static double quadform(double[][] Q, DoublePoint w) {
        double[] v = {w.X, w.Y, 1};
        int i, j;
        double sum = 0;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                sum += v[i] * Q[i][j] * v[j];
            }
        }

        return sum;
    }


    /**
     * Range over the straight line segment [a,b] when lambda ranges over [0,1].
     *
     * @param lambda Scale
     * @param a      Start Point
     * @param b      Stop Point
     * @return Point on the segment
     */
    static DoublePoint interval(double lambda, DoublePoint a, DoublePoint b) {
        return new DoublePoint(a.X + lambda * (b.X - a.X), a.Y + lambda * (b.Y - a.Y));
    }


    /* ddenom/dpara have the property that the square of radius 1 centered
           at p1 intersects the line p0p2 iff |dpara(p0,p1,p2)| <= ddenom(p0,p2) */
    static double ddenom(DoublePoint p0, DoublePoint p2) {
        IntPoint r = dorth_infty(p0, p2);

        return r.Y * (p2.X - p0.X) - r.X * (p2.Y - p0.Y);
    }


    /**
     * return a direction that is 90 degrees counterclockwise from p2-p0,
     * but then restricted to one of the major wind directions (n, nw, w, etc).
     *
     * @param p0
     * @param p2
     * @return
     */
    static IntPoint dorth_infty(DoublePoint p0, DoublePoint p2) {

        return new IntPoint(-sign(p2.Y - p0.Y), sign(p2.X - p0.X));
    }


    /* return (p1-p0)x(p2-p0), the area of the parallelogram */
    static double dpara(DoublePoint p0, DoublePoint p1, DoublePoint p2) {

        double x1, y1, x2, y2;

        x1 = p1.X - p0.X;
        y1 = p1.Y - p0.Y;
        x2 = p2.X - p0.X;
        y2 = p2.Y - p0.Y;

        return x1 * y2 - x2 * y1;
    }


    /**
     * Calculate distance between two points.
     *
     * @param p
     * @param q
     * @return
     */
    static double ddist(DoublePoint p, DoublePoint q) {
        return java.lang.Math.sqrt((p.X - q.X) * (p.X - q.X) + (p.Y - q.Y) * (p.Y - q.Y));
    }


    /* calculate (p1-p0)x(p3-p2) */
    static double cprod(DoublePoint p0, DoublePoint p1, DoublePoint p2, DoublePoint p3) {
        double x1, y1, x2, y2;

        x1 = p1.X - p0.X;
        y1 = p1.Y - p0.Y;
        x2 = p3.X - p2.X;
        y2 = p3.Y - p2.Y;

        return x1 * y2 - x2 * y1;
    }


    /* calculate (p1-p0)*(p3-p2) */
    static double iprod1(DoublePoint p0, DoublePoint p1, DoublePoint p2, DoublePoint p3) {
        double x1, y1, x2, y2;

        x1 = p1.X - p0.X;
        y1 = p1.Y - p0.Y;
        x2 = p3.X - p2.X;
        y2 = p3.Y - p2.Y;

        return x1 * x2 + y1 * y2;
    }


    /**
     * Calculates the point t in [0..1] on the (convex) bezier curve
     * (p0,p1,p2,p3) which is tangent to q1-q0. Return -1.0 if there is no
     * solution in [0..1].
     *
     * @param p0
     * @param p1
     * @param p2
     * @param p3
     * @param q0
     * @param q1
     * @return
     */
    static double Tangent(
            DoublePoint p0, DoublePoint p1, DoublePoint p2, DoublePoint p3,
            DoublePoint q0, DoublePoint q1) {

        double A, B, C;   /* (1-t)^2 A + 2(1-t)t B + t^2 C = 0 */
        double a, b, c;   /* a t^2 + b t + c = 0 */
        double d, s, r1, r2;

        A = cprod(p0, p1, q0, q1);
        B = cprod(p1, p2, q0, q1);
        C = cprod(p2, p3, q0, q1);

        a = A - 2 * B + C;
        b = -2 * A + 2 * B;
        c = A;

        d = b * b - 4 * a * c;

        if (a == 0 || d < 0) {
            return -1.0;
        }

        s = java.lang.Math.sqrt(d);

        r1 = (-b + s) / (2 * a);
        r2 = (-b - s) / (2 * a);

        if (r1 >= 0 && r1 <= 1) {
            return r1;
        } else if (r2 >= 0 && r2 <= 1) {
            return r2;
        } else {
            return -1.0;
        }
    }


    /* calculate (p1-p0)*(p2-p0) */
    static double iprod(DoublePoint p0, DoublePoint p1, DoublePoint p2) {
        double x1, y1, x2, y2;

        x1 = p1.X - p0.X;
        y1 = p1.Y - p0.Y;
        x2 = p2.X - p0.X;
        y2 = p2.Y - p0.Y;

        return x1 * x2 + y1 * y2;
    }


    /**
     * Calculates a point of a bezier curve.
     *
     * @param t
     * @param p0
     * @param p1
     * @param p2
     * @param p3
     * @return
     */
    static DoublePoint Bezier(double t, DoublePoint p0, DoublePoint p1, DoublePoint p2, DoublePoint p3) {
        double s = 1 - t;
        double x = s * s * s * p0.X + 3 * (s * s * t) * p1.X + 3 * (t * t * s) * p2.X + t * t * t * p3.X;
        double y = s * s * s * p0.Y + 3 * (s * s * t) * p1.Y + 3 * (t * t * s) * p2.Y + t * t * t * p3.Y;
        return new DoublePoint(x, y);
    }
}
