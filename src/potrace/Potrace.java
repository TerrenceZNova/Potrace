package potrace;


import org.opencv.core.Mat;
import org.opencv.core.MatOfPoint;
import org.opencv.imgcodecs.Imgcodecs;
import org.opencv.imgproc.Imgproc;
import org.opencv.core.MatOfPoint;


import java.awt.*;
import java.util.List;
import java.util.ArrayList;
import java.util.LinkedList;


public class Potrace {


    /**
     * This function converts source image into a boolean array
     *
     * @param srcImage
     * @return the boolean array where true stands for background and false
     * stands for foreground
     * @apiNote the output array has larger size than source image cause
     * background border is added around the image by 1 pixel
     */
    public static boolean[][] bitmapToBinary(Mat srcImage) {

        assert (srcImage != null && !srcImage.empty());

        Mat binaryImage = new Mat();
        if (srcImage.channels() == 3) {
            Mat grayImage = new Mat();
            Imgproc.cvtColor(srcImage, grayImage, Imgproc.COLOR_BGR2GRAY);
            Imgproc.threshold(grayImage, binaryImage,
                    207, 255, Imgproc.THRESH_BINARY);
        } else {
            Imgproc.threshold(srcImage, binaryImage,
                    207, 255, Imgproc.THRESH_BINARY);
        }

        boolean[][] result = new boolean[srcImage.rows() + 2][srcImage.cols() + 2];
        for (int row = 0; row < binaryImage.rows(); row++) {
            for (int col = 0; col < binaryImage.cols(); col++) {
                if ((int) binaryImage.get(row, col)[0] == 255) {
                    result[row + 1][col + 1] = true;
                } else {
                    result[row + 1][col + 1] = false;
                }
            }
        }

        // add white borders
        for (int col = 0; col < result[0].length; col++) {
            result[0][col] = true;
            result[result.length - 1][col] = true;
        }
        for (int row = 1; row < result.length - 1; row++) {
            result[row][0] = true;
            result[row][result[0].length - 1] = true;
        }

        return result;
    }


    /**
     * This function attempts to trace the given bitmap using the
     * given tracing parameters.
     *
     * @param bm                the given bitmap in binary array form
     * @param ListOfCurveArrays the curve result of Potrace
     */
    public static void potrace_trace(boolean[][] bm, ArrayList ListOfCurveArrays) {

        //Step 1
        ArrayList plistp = new ArrayList();
        bm_to_pathlist(bm, plistp);

        //Step 2
        process_path(plistp);

        //Step 3
        PathList_to_ListOfCurveArrays(plistp, ListOfCurveArrays);
    }


    /**
     * Decompose the given bitmap into paths. Returns a list of Path
     * objects with the fields len, pt, area filled
     *
     * @param bm     A binary bitmap which holds the imageinformations.
     * @param plistp List of Path objects
     */
    private static void bm_to_pathlist(boolean[][] bm, ArrayList<ArrayList<Path>> plistp) {

        Integer x = 0, y = 0;
        Point point;

        while (true) {
            point = findNext(bm, x, y);
            if (point == null) {
                break;
            } else {
                getContour(bm, point.x, point.y, plistp);
            }
        }

    }


    /**
     * Get contour from binary matrix
     *
     * @param bm
     * @param x
     * @param y
     * @param plistp
     */
    private static void getContour(boolean[][] bm, Integer x, Integer y,
                                   ArrayList<ArrayList<Path>> plistp) {

        //Step 1: Find Path
        Path Contur = findpath(bm, new IntPoint(x, y));

        //Step 2:
        Xor_Path(bm, Contur);

        outputMatrix(bm);

        ArrayList<Path> PolyPath = new ArrayList<Path>();

        // only area > turdsize is taken
        if (Contur.area > Constant.turdsize) {
            PolyPath.add(Contur); // Path with index 0 is a conture
            plistp.add(PolyPath);
        }


        Point point;
        while (true) {
            point = findNext(bm, x, y, Contur);
            if (point == null) {
                break;
            }
            x = point.x;
            y = point.y;


            Path Hole = findpath(bm, new IntPoint(x, y));
            Xor_Path(bm, Hole);

            if (Hole.area > Constant.turdsize) {
                PolyPath.add(Hole); // Path with index > 0 is a hole,
            }

            while (true) {    // 13.07.12 von if auf while
                point = findNext(bm, x, y, Hole);
                if (point == null) {
                    break;
                }
                x = point.x;
                y = point.y;
                getContour(bm, x, y, plistp);
            }
        }
    }


    /**
     * Searches a x and a y such that source[x,y] = true and source[x+1,y] false.
     * If this not exists, false will be returned else the result is true.
     *
     * @param bm a binary matrix
     * @param x  x index in the source Matrix
     * @param y  y index in the source Matrix
     * @return
     */
    private static Point findNext(boolean[][] bm, Integer x, Integer y) {

        for (int col = 1; col < bm[0].length - 1; col++) {
            for (int row = 0; row < bm.length - 1; row++) {
                // black found
                if (!bm[row + 1][col]) {
                    x = row;
                    y = col;
                    return new Point(row, col);
                }
            }
        }

        return null;

    }


    /**
     * Searches a x and a y inside the Path P such that source[x,y] = true and source[x+1,y] false.
     * If this not exists, false will be returned else the result is true.
     *
     * @param Matrix a Binary Matrix
     * @param x      x index in the source Matrix
     * @param y      y index in the source Matrix
     * @param P
     * @return
     */
    private static Point findNext(boolean[][] Matrix, int x, int y, Path P) {

        int i = 0;
        int n = P.pt.length;
        ArrayList<MonotonInterval> MonotonIntervals = P.MonotonIntervals;
        if (MonotonIntervals.size() == 0) {
            return null;
        }

        MonotonInterval MI = MonotonIntervals.get(0);
        MI.ResetCurrentID(n);
        y = P.pt[MI.CurrentID].Y;
        ArrayList<MonotonInterval> CurrentIntervals = new ArrayList<MonotonInterval>();
        CurrentIntervals.add(MI);
        MI.CurrentID = MI.Min();

        while (i + 1 < MonotonIntervals.size()
                && (MonotonIntervals.get(i + 1)).MinY(P.pt) == y) {

            MI = MonotonIntervals.get(i + 1);
            MI.ResetCurrentID(n);
            CurrentIntervals.add(MI);
            i++;
        }

        while (CurrentIntervals.size() > 0) {

            for (int k = 0; k < CurrentIntervals.size() - 1; k++) {
                int x1 = P.pt[(CurrentIntervals.get(k)).CurrentID].X + 1;
                int x2 = P.pt[(CurrentIntervals.get(k + 1)).CurrentID].X;

                for (x = x1; x <= x2; x++)
                    if (!Matrix[x][y]) {
                        x--;
                        return new Point(x, y);
                    }
                k++;
            }

            y++;
            for (int j = CurrentIntervals.size() - 1; j >= 0; j--) {

                MonotonInterval M = CurrentIntervals.get(j);

                if (y > M.MaxY(P.pt)) {
                    CurrentIntervals.remove(j);
                    continue;
                }
                int CID = M.CurrentID;
                do {
                    if (M.Increasing) {
                        CID = Math.mod(CID + 1, n);
                    } else {
                        CID = Math.mod(CID - 1, n);
                    }
                } while (P.pt[CID].Y < y);

                M.CurrentID = CID;
            }

            // Add Items of MonotonIntervals with Miny==y
            while (i + 1 < MonotonIntervals.size() && (MonotonIntervals.get(i + 1)).MinY(P.pt) == y) {

                MonotonInterval NewInt = MonotonIntervals.get(i + 1);
                int j = 0;
                // search the correct x-Position
                int _x = P.pt[NewInt.Min()].X;

                while (j < CurrentIntervals.size()
                        && (_x > P.pt[(CurrentIntervals.get(j)).CurrentID].X)) {
                    j++;
                }
                CurrentIntervals.add(j, NewInt);
                NewInt.ResetCurrentID(n);
                i++;
            }
        }

        return null;
    }


    private static void process_path(ArrayList<ArrayList<Path>> plistp) {

        Path p;
        ArrayList<Path> plist;

        /* call downstream function with each path */
        for (int j = 0; j < plistp.size(); j++) {
            plist = plistp.get(j);

            for (int i = 0; i < plist.size(); i++) {
                p = plist.get(i);
                calc_sums(p);
                calc_lon(p);
                BestPolygon(p);
                AdjustVertices(p);
                smooth(p.Curves, 1, Constant.alphamax);
                if (Constant.curveoptimizing) {
                    opticurve(p, Constant.opttolerance);
                    p.FCurves = p.OptimizedCurves;
                } else {
                    p.FCurves = p.Curves;
                }
                p.Curves = p.FCurves;
            }
        }
    }


    /**
     * Preparation: fill in the sum* fields of a path (used for later
     * rapid summing).
     *
     * @param pp Path for which the preparation will be done
     */
    private static void calc_sums(Path pp) {

        int i, x, y;
        int n = pp.pt.length;
        pp.Sums = new SumStruct[n + 1];

        // origin
        int x0 = pp.pt[0].X;
        int y0 = pp.pt[0].Y;

        // preparatory computation for later fast summing
        //pp->sums[0].x2 = pp->sums[0].xy = pp->sums[0].y2 = pp->sums[0].x = pp->sums[0].y = 0;
        pp.Sums[0].X2 = pp.Sums[0].XY = pp.Sums[0].Y2 = pp.Sums[0].X = pp.Sums[0].Y = 0;

        for (i = 0; i < n; i++) {
            x = pp.pt[i].X - x0;
            y = pp.pt[i].Y - y0;
            pp.Sums[i + 1].X = pp.Sums[i].X + x;
            pp.Sums[i + 1].Y = pp.Sums[i].Y + y;
            pp.Sums[i + 1].X2 = pp.Sums[i].X2 + x * x;
            pp.Sums[i + 1].XY = pp.Sums[i].XY + x * y;
            pp.Sums[i + 1].Y2 = pp.Sums[i].Y2 + y * y;
        }
    }

    private static void calc_lon(Path pp) {

        int i, j, k, k1;
        int a, b, c, d;
        int[] ct = {0, 0, 0, 0};
        int dir;
        IntPoint[] constraint = new IntPoint[2];
        IntPoint cur;
        IntPoint off;
        IntPoint dk;  /* direction of k-k1 */
        IntPoint[] pt = pp.pt;


        int n = pt.length;
        int[] Pivot = new int[n];
        int[] nc = new int[n];

        /* initialize the nc data structure. Point from each point to the
           furthest future point to which it is connected by a vertical or
           horizontal segment. We take advantage of the fact that there is
           always a direction change at 0 (due to the path decomposition
           algorithm). But even if this were not so, there is no harm, as
           in practice, correctness does not depend on the word "furthest"
           above.  */

        k = 0;
        for (i = n - 1; i >= 0; i--) {
            if (pt[i].X != pt[k].X && pt[i].Y != pt[k].Y) {
                k = i + 1;  /* necessarily i<n-1 in this case */
            }
            nc[i] = k;
        }

        pp.Lon = new int[n];

        // determine pivot points: for each i, let pivk[i] be the furthest k
        // such that all j with i<j<k lie on a line connecting i,k.
        for (i = n - 1; i >= 0; i--) {

            boolean foundkFlag = false;

            ct[0] = ct[1] = ct[2] = ct[3] = 0;

            /* keep track of "directions" that have occurred */
            dir = (3 + 3 * (pt[Math.mod(i + 1, n)].X - pt[i].X)
                    + (pt[Math.mod(i + 1, n)].Y - pt[i].Y)) / 2;
            ct[dir]++;


            /* find the next k such that no straight line from i to k */
            k = nc[i];
            k1 = i;
            while (true) {
                dir = (3 + 3 * Math.sign(pt[k].X - pt[k1].X) + Math.sign(pt[k].Y - pt[k1].Y)) / 2;
                ct[dir]++;

                /* if all four "directions" have occurred, cut this path */
                if ((ct[0] == 1) && (ct[1] == 1) && (ct[2] == 1) && (ct[3] == 1)) {
                    Pivot[i] = k1;
                    //goto foundk;
                    foundkFlag = true;
                    break;
                }

                cur = new IntPoint(pt[k].X - pt[i].X, pt[k].Y - pt[i].Y);

                /* see if current constraint is violated */
                if (Math.xprod(constraint[0], cur) < 0
                        || Math.xprod(constraint[1], cur) > 0) {
                    break;
                }

                /* else, update constraint */
                if (Math.abs(cur.X) <= 1 && Math.abs(cur.Y) <= 1) {
                    /* no constraint */
                } else {
                    off = new IntPoint(cur.X + ((cur.Y >= 0 && (cur.Y > 0 || cur.X < 0)) ? 1 : -1), cur.Y + ((cur.X <= 0 && (cur.X < 0 || cur.Y < 0)) ? 1 : -1));
                    if (Math.xprod(constraint[0], off) >= 0) {
                        constraint[0] = off;
                    }
                    off = new IntPoint(cur.X + ((cur.Y <= 0 && (cur.Y < 0 || cur.X < 0)) ? 1 : -1), cur.Y + ((cur.X >= 0 && (cur.X > 0 || cur.Y < 0)) ? 1 : -1));
                    if (Math.xprod(constraint[1], off) <= 0) {
                        constraint[1] = off;
                    }
                }
                k1 = k;
                k = nc[k1];
                if (!Math.cyclic(k, i, k1)) {
                    break;
                }
            }
            if (foundkFlag) {
                continue;
            }

            //constraint_viol:

            /* k1 was the last "corner" satisfying the current constraint, and
               k is the first one violating it. We now need to find the last
               point along k1..k which satisfied the constraint. */

            dk = new IntPoint(Math.sign(pt[k].X - pt[k1].X), Math.sign(pt[k].Y - pt[k1].Y));
            cur = new IntPoint(pt[k1].X - pt[i].X, pt[k1].Y - pt[i].Y);

            /* find largest integer j such that xprod(constraint[0], cur+j*dk)
               >= 0 and xprod(constraint[1], cur+j*dk) <= 0. Use bilinearity
               of xprod. */
            a = Math.xprod(constraint[0], cur);
            b = Math.xprod(constraint[0], dk);
            c = Math.xprod(constraint[1], cur);
            d = Math.xprod(constraint[1], dk);

            /* find largest integer j such that a+j*b>=0 and c+j*d<=0. This
               can be solved with integer arithmetic. */
            j = Integer.MAX_VALUE;
            if (b < 0) {
                j = Math.floordiv(a, -b);
            }
            if (d > 0) {
                j = Math.min(j, Math.floordiv(-c, d));
            }
            Pivot[i] = Math.mod(k1 + j, n);
            //foundk:
        }

        /* for i */
        /* clean up: for each i, let lon[i] be the largest k such that for
           all i' with i<=i'<k, i'<k<=pivk[i']. */

        j = Pivot[n - 1];
        pp.Lon[n - 1] = j;

        for (i = n - 2; i >= 0; i--) {
            if (Math.cyclic(i + 1, Pivot[i], j)) {
                j = Pivot[i];
            }
            pp.Lon[i] = j;

        }

        for (i = n - 1; Math.cyclic(Math.mod(i + 1, n), j, pp.Lon[i]); i--) {
            pp.Lon[i] = j;
        }
    }

    /**
     * find the optimal polygon. Fill in the m and po components. Return 1
     * on failure with errno set, else 0. Non-cyclic version: assumes i=0
     * is in the polygon.
     * Fixme implement cyclic version.
     *
     * @param pp
     */
    private static void BestPolygon(Path pp) {

        int i, j, m, k;
        int n = pp.pt.length;

        double[] pen = new double[n + 1]; /* pen[n+1]: penalty vector */
        int[] prev = new int[n + 1];   /* prev[n+1]: best path pointer vector */
        int[] clip0 = new int[n];  /* clip0[n]: longest segment pointer, non-cyclic */
        int[] clip1 = new int[n + 1];  /* clip1[n+1]: backwards segment pointer, non-cyclic */
        int[] seg0 = new int[n + 1];    /* seg0[m+1]: forward segment bounds, m<=n */
        int[] seg1 = new int[n + 1];   /* seg1[m+1]: backward segment bounds, m<=n */

        double thispen;
        double best;
        int c;

        /* calculate clipped paths */
        for (i = 0; i < n; i++) {
            c = Math.mod(pp.Lon[Math.mod(i - 1, n)] - 1, n);

            if (c == i) {
                c = Math.mod(i + 1, n);
            }
            if (c < i) {
                clip0[i] = n;
            } else {
                clip0[i] = c;
            }
        }

        /* calculate backwards path clipping, non-cyclic. j <= clip0[i] if
           clip1[j] <= i, for i,j=0..n. */
        j = 1;
        for (i = 0; i < n; i++) {
            while (j <= clip0[i]) {
                clip1[j] = i;
                j++;
            }
        }

        /* calculate seg0[j] = longest path from 0 with j segments */
        i = 0;
        for (j = 0; i < n; j++) {
            seg0[j] = i;
            i = clip0[i];
        }
        seg0[j] = n;
        m = j;

        /* calculate seg1[j] = longest path to n with m-j segments */
        i = n;
        for (j = m; j > 0; j--) {
            seg1[j] = i;
            i = clip1[i];
        }
        seg1[0] = 0;

        /* now find the shortest path with m segments, based on penalty3 */
        /* note: the outer 2 loops jointly have at most n interations, thus
           the worst-case behavior here is quadratic. In practice, it is
           close to linear since the inner loop tends to be short. */
        pen[0] = 0;
        for (j = 1; j <= m; j++) {
            for (i = seg1[j]; i <= seg0[j]; i++) {
                best = -1;
                for (k = seg0[j - 1]; k >= clip1[i]; k--) {
                    thispen = penalty3(pp, k, i) + pen[k];
                    if (best < 0 || thispen < best) {
                        prev[i] = k;
                        best = thispen;
                    }
                }
                pen[i] = best;
            }
        }


        /* read off shortest path */

        int[] B = new int[m];

        pp.po = new int[m];
        for (i = n, j = m - 1; i > 0; j--) {
            i = prev[i];
            B[j] = i;
        }

        /*if ((m > 0) && (mod(pp.Lon[m - 1] - 1, n) <= B[1])) {// reduce
            B[0] = B[m - 1];
            pp.po = new int[m - 1];
            for (i = 0; i < m - 1; i++)
                pp.po[i] = B[i];

        } else {
        }*/

        pp.po = B;
    }

    /**
     * Stage 3: vertex adjustment (Sec. 2.3.1).
     * <p>
     * Calculate "optimal" point-slope representation for each line segment.
     * Adjust vertices of optimal polygon: calculate the intersection of
     * the two "optimal" line segments, then move it into the unit square
     * if it lies outside. Return 1 with errno set on error; 0 on success.
     *
     * @param pp
     */
    private static void AdjustVertices(Path pp) {

        int m = pp.po.length;
        int[] po = pp.po;
        IntPoint[] pt = pp.pt;
        int n = pt.length;


        int x0 = pt[0].X;
        int y0 = pt[0].Y;

        DoublePoint[] ctr = new DoublePoint[m];      /* ctr[m] */
        DoublePoint[] dir = new DoublePoint[m];      /* dir[m] */

        double[][][] q = new double[m][3][3];
        double[] v = new double[3];
        double d;
        int i, j, k, l;
        DoublePoint s;
        pp.Curves = new PrivCurve(m);

        // Calculate "optimal" point-slope representation for each line segment.
        for (i = 0; i < m; i++) {
            j = po[Math.mod(i + 1, m)];
            j = Math.mod(j - po[i], n) + po[i];
            Math.pointslope(pp, po[i], j, ctr[i], dir[i]);
        }
        // Represents each line segment as a singular quadratic form; the
        // distance of a point (x,y) from the line segment will be
        // (x,y,1)Q(x,y,1)^t, where Q=q[i].
        for (i = 0; i < m; i++) {
            d = dir[i].X * dir[i].X + dir[i].Y * dir[i].Y;

            if (d == 0.0) {
                for (j = 0; j < 3; j++) {
                    for (k = 0; k < 3; k++) {
                        q[i][j][k] = 0;
                    }
                }
            } else {
                v[0] = dir[i].Y;
                v[1] = -dir[i].X;
                v[2] = -v[1] * ctr[i].Y - v[0] * ctr[i].X;
                for (l = 0; l < 3; l++) {
                    for (k = 0; k < 3; k++) {
                        q[i][l][k] = v[l] * v[k] / d;
                    }
                }
            }
        }

        /* now calculate the "intersections" of consecutive segments.
               Instead of using the actual intersection, we find the point
               within a given unit square which minimizes the square distance to
               the two lines. */
        for (i = 0; i < m; i++) {
            double[][] Q = new double[3][3];
            DoublePoint w;
            double dx, dy;
            double det;
            double min, cand; /* minimum and candidate for minimum of quad. form */
            double xmin, ymin;    /* coordinates of minimum */
            int z;

            /* let s be the vertex, in coordinates relative to x0/y0 */
            s = new DoublePoint(pt[po[i]].X - x0, pt[po[i]].Y - y0);

            /* intersect segments i-1 and i */

            j = Math.mod(i - 1, m);

            /* add quadratic forms */
            for (l = 0; l < 3; l++) {
                for (k = 0; k < 3; k++) {
                    Q[l][k] = q[j][l][k] + q[i][l][k];
                }
            }
            while (true) {
                /* minimize the quadratic form Q on the unit square */
                /* find intersection */
                det = Q[0][0] * Q[1][1] - Q[0][1] * Q[1][0];
                if (det != 0.0) {
                    double wx = (-Q[0][2] * Q[1][1] + Q[1][2] * Q[0][1]) / det;
                    double wy = (Q[0][2] * Q[1][0] - Q[1][2] * Q[0][0]) / det;
                    w = new DoublePoint(wx, wy);
                    break;
                }

                /* matrix is singular - lines are parallel. Add another,
                   orthogonal axis, through the center of the unit square */

                if (Q[0][0] > Q[1][1]) {

                    v[0] = -Q[0][1];
                    v[1] = Q[0][0];

                } else if (Q[1][1] != 0) {
                    // nur if (Q[1,1])
                    v[0] = -Q[1][1];
                    v[1] = Q[1][0];

                } else {

                    v[0] = 1;
                    v[1] = 0;
                }

                d = v[0] * v[0] + v[1] * v[1];
                v[2] = -v[1] * s.Y - v[0] * s.X;
                for (l = 0; l < 3; l++) {
                    for (k = 0; k < 3; k++) {
                        Q[l][k] += v[l] * v[k] / d;
                    }
                }
            }
            dx = java.lang.Math.abs(w.X - s.X);
            dy = java.lang.Math.abs(w.Y - s.Y);

            if (dx <= .5 && dy <= .5) {
                // - 1 because we have a additional border set to the bitmap
                pp.Curves.Vertex[i] = new DoublePoint(w.X + x0, w.Y + y0);
                continue;
            }

            // the minimum was not in the unit square; now minimize quadratic
            // on boundary of square.
            min = Math.quadform(Q, s);
            xmin = s.X;
            ymin = s.Y;

            boolean fixxFlag = false;
            if (Q[0][0] == 0.0) {
                //goto fixx;
                fixxFlag = true;
            }

            if (!fixxFlag) {
                for (z = 0; z < 2; z++) {
                    // value of the y-coordinate.
                    double wY = s.Y - 0.5 + z;
                    double wX = -(Q[0][1] * wY + Q[0][2]) / Q[0][0];
                    w = new DoublePoint(wX, wY);
                    dx = java.lang.Math.abs(wX - s.X);
                    cand = Math.quadform(Q, w);
                    if (dx <= .5 && cand < min) {
                        min = cand;
                        xmin = w.X;
                        ymin = w.Y;
                    }
                }
            }

            //fixx:
            boolean cornersFlag = false;
            if (Q[1][1] == 0.0) {
                //goto corners;
                cornersFlag = true;
            }
            if (!cornersFlag) {
                for (z = 0; z < 2; z++) {   /* value of the x-coordinate */
                    double wX = s.X - 0.5 + z;
                    double wY = -(Q[1][0] * wX + Q[1][2]) / Q[1][1];
                    w = new DoublePoint(wX, wY);
                    dy = java.lang.Math.abs(wY - s.Y);
                    cand = Math.quadform(Q, w);
                    if (dy <= .5 && cand < min) {
                        min = cand;
                        xmin = w.X;
                        ymin = w.Y;
                    }
                }
            }

            //corners:
            /* check four corners */
            for (l = 0; l < 2; l++) {
                for (k = 0; k < 2; k++) {
                    w = new DoublePoint(s.X - 0.5 + l, s.Y - 0.5 + k);
                    cand = Math.quadform(Q, w);
                    if (cand < min) {
                        min = cand;
                        xmin = w.X;
                        ymin = w.Y;
                    }
                }
            }
            // - 1 because we have a additional border set to the bitmap
            pp.Curves.Vertex[i] = new DoublePoint(xmin + x0 - 1, ymin + y0 - 1);

            //continue;
        }
    }

    /**
     * Stage 4: smoothing and corner analysis (Sec. 2.3.3)
     * Always succeeds and returns 0
     *
     * @param curve
     * @param sign
     * @param alphamax
     */
    private static void smooth(PrivCurve curve, int sign, double alphamax) {

        int m = curve.Count;

        int i, j, k;
        double dd, denom, alpha;
        DoublePoint p2, p3, p4;

        if (sign == '-') {
            /* reverse orientation of negative paths */
            for (i = 0, j = m - 1; i < j; i++, j--) {
                DoublePoint tmp;
                tmp = curve.Vertex[i];
                curve.Vertex[i] = curve.Vertex[j];
                curve.Vertex[j] = tmp;
            }
        }

        /* examine each vertex and find its best fit */
        for (i = 0; i < m; i++) {
            j = Math.mod(i + 1, m);
            k = Math.mod(i + 2, m);
            p4 = Math.interval(1 / 2.0, curve.Vertex[k], curve.Vertex[j]);

            denom = Math.ddenom(curve.Vertex[i], curve.Vertex[k]);
            if (denom != 0.0) {
                dd = Math.dpara(curve.Vertex[i], curve.Vertex[j], curve.Vertex[k]) / denom;
                dd = java.lang.Math.abs(dd);
                alpha = dd > 1 ? (1 - 1.0 / dd) : 0;
                alpha = alpha / 0.75;
            } else {
                alpha = 4 / 3.0;
            }
            curve.Alpha0[j] = alpha;     /* remember "original" value of alpha */

            if (alpha > alphamax) {  /* pointed corner */
                curve.Tag[j] = Constant.POTRACE_CORNER;
                //curve.c[j][1] = curve->vertex[j];
                curve.ControlPoints[j][1] = curve.Vertex[j];
                curve.ControlPoints[j][2] = p4;
            } else {
                if (alpha < 0.55) {
                    alpha = 0.55;
                } else if (alpha > 1) {
                    alpha = 1;
                }
                p2 = Math.interval(.5 + .5 * alpha, curve.Vertex[i], curve.Vertex[j]);
                p3 = Math.interval(.5 + .5 * alpha, curve.Vertex[k], curve.Vertex[j]);
                curve.Tag[j] = Constant.POTRACE_CURVETO;
                curve.ControlPoints[j][0] = p2;
                curve.ControlPoints[j][1] = p3;
                curve.ControlPoints[j][2] = p4;
            }
            curve.Alpha[j] = alpha;    /* store the "cropped" value of alpha */
            curve.Beta[j] = 0.5;
        }
    }


    /**
     * optimize the path p, replacing sequences of Bezier segments by a
     * single segment when possible. Return 0 on success, 1 with errno set
     * on failure.
     *
     * @param pp
     * @param opttolerance
     */
    private static void opticurve(Path pp, double opttolerance) {

        int m = pp.Curves.Count;
        int[] pt = new int[m + 1];     /* pt[m+1] */
        double[] pen = new double[m + 1];  /* pen[m+1] */
        int[] len = new int[m + 1];     /* len[m+1] */
        opti[] opt = new opti[m + 1];    /* opt[m+1] */
        int[] convc = new int[m];       /* conv[m]: pre-computed convexities */
        double[] areac = new double[m + 1];  /* cumarea[m+1]: cache for fast area computation */

        int om;
        int i, j;
        boolean r;
        opti o = new opti();
        DoublePoint p0;
        int i1;
        double area;
        double alpha;
        double[] s;
        double[] t;

        /* pre-calculate convexity: +1 = right turn, -1 = left turn, 0 = corner */
        for (i = 0; i < m; i++) {
            if (pp.Curves.Tag[i] == Constant.POTRACE_CURVETO) {
                convc[i]
                        = Math.sign(
                        Math.dpara(
                                pp.Curves.Vertex[Math.mod(i - 1, m)],
                                pp.Curves.Vertex[i],
                                pp.Curves.Vertex[Math.mod(i + 1, m)]));
            } else {
                convc[i] = 0;
            }
        }

        /* pre-calculate areas */
        area = 0.0;
        areac[0] = 0.0;
        p0 = pp.Curves.Vertex[0];
        for (i = 0; i < m; i++) {
            i1 = Math.mod(i + 1, m);
            if (pp.Curves.Tag[i1] == Constant.POTRACE_CURVETO) {
                alpha = pp.Curves.Alpha[i1];
                area +=
                        0.3 * alpha * (4 - alpha) * Math.dpara(
                                pp.Curves.ControlPoints[i][2],
                                pp.Curves.Vertex[i1], pp.Curves.ControlPoints[i1][2]) / 2;
                area +=
                        Math.dpara(p0,
                                pp.Curves.ControlPoints[i][2],
                                pp.Curves.ControlPoints[i1][2]) / 2;
            }
            areac[i + 1] = area;
        }

        pt[0] = -1;
        pen[0] = 0;
        len[0] = 0;

        /* Fixme: we always start from a fixed point -- should find the best
          curve cyclically ### */

        for (j = 1; j <= m; j++) {
            /* calculate best path from 0 to j */
            pt[j] = j - 1;
            pen[j] = pen[j - 1];
            len[j] = len[j - 1] + 1;

            for (i = j - 2; i >= 0; i--) {
                r = opti_penalty(pp, i, Math.mod(j, m), o, opttolerance, convc, areac);
                if (r) {
                    break;
                }
                if (len[j] > len[i] + 1 || (len[j] == len[i] + 1 && pen[j] > pen[i] + o.pen)) {
                    pt[j] = i;
                    pen[j] = pen[i] + o.pen;
                    len[j] = len[i] + 1;
                    opt[j] = o;
                }
            }
        }
        om = len[m];
        pp.OptimizedCurves = new PrivCurve(om);

        s = new double[om];
        t = new double[om];


        j = m;
        for (i = om - 1; i >= 0; i--) {
            if (pt[j] == j - 1) {
                pp.OptimizedCurves.Tag[i] = pp.Curves.Tag[Math.mod(j, m)];

                pp.OptimizedCurves.ControlPoints[i][0]
                        = pp.Curves.ControlPoints[Math.mod(j, m)][0];
                pp.OptimizedCurves.ControlPoints[i][1]
                        = pp.Curves.ControlPoints[Math.mod(j, m)][1];
                pp.OptimizedCurves.ControlPoints[i][2]
                        = pp.Curves.ControlPoints[Math.mod(j, m)][2];
                pp.OptimizedCurves.Vertex[i]
                        = pp.Curves.Vertex[Math.mod(j, m)];
                pp.OptimizedCurves.Alpha[i]
                        = pp.Curves.Alpha[Math.mod(j, m)];
                pp.OptimizedCurves.Alpha0[i]
                        = pp.Curves.Alpha0[Math.mod(j, m)];
                pp.OptimizedCurves.Beta[i]
                        = pp.Curves.Beta[Math.mod(j, m)];
                s[i] = t[i] = 1.0;
            } else {
                pp.OptimizedCurves.Tag[i]
                        = Constant.POTRACE_CURVETO;
                pp.OptimizedCurves.ControlPoints[i][0]
                        = opt[j].c[0];
                pp.OptimizedCurves.ControlPoints[i][1]
                        = opt[j].c[1];
                pp.OptimizedCurves.ControlPoints[i][2]
                        = pp.Curves.ControlPoints[Math.mod(j, m)][2];
                pp.OptimizedCurves.Vertex[i]
                        = Math.interval(
                        opt[j].s, pp.Curves.ControlPoints[Math.mod(j, m)][2],
                        pp.Curves.Vertex[Math.mod(j, m)]);
                pp.OptimizedCurves.Alpha[i]
                        = opt[j].alpha;
                pp.OptimizedCurves.Alpha0[i]
                        = opt[j].alpha;

                s[i] = opt[j].s;
                t[i] = opt[j].t;
            }

            j = pt[j];
        }

        /* calculate beta parameters */
        for (i = 0; i < om; i++) {
            i1 = Math.mod(i + 1, om);
            pp.OptimizedCurves.Beta[i] = s[i] / (s[i] + t[i1]);
        }
    }


    /* calculate best fit from i+.5 to j+.5.  Assume i<j (cyclically).
           Return 0 and set badness and parameters (alpha, beta), if
           possible. Return 1 if impossible. */
    private static boolean opti_penalty(
            Path pp, int i, int j, opti res, double opttolerance, int[] convc, double[] areac) {

        int m = pp.Curves.Count;
        int k, k1, k2, conv, i1;
        double area, alpha, d, d1, d2;
        DoublePoint p0, p1, p2, p3, pt;
        double A, R, A1, A2, A3, A4;
        double s, t;

        /* check convexity, corner-freeness, and maximum bend < 179 degrees */

        if (i == j) {  /* sanity - a full loop can never be an opticurve */
            return true;
        }

        k = i;
        i1 = Math.mod(i + 1, m);
        k1 = Math.mod(k + 1, m);
        conv = convc[k1];

        if (conv == 0) {
            return true;
        }

        d = Math.ddist(pp.Curves.Vertex[i], pp.Curves.Vertex[i1]);

        for (k = k1; k != j; k = k1) {
            k1 = Math.mod(k + 1, m);
            k2 = Math.mod(k + 2, m);

            if (convc[k1] != conv) {
                return true;
            }

            if (
                    Math.sign(Math.cprod(
                            pp.Curves.Vertex[i],
                            pp.Curves.Vertex[i1],
                            pp.Curves.Vertex[k1],
                            pp.Curves.Vertex[k2])) != conv) {
                return true;
            }
            if (Math.iprod1(pp.Curves.Vertex[i],
                    pp.Curves.Vertex[i1],
                    pp.Curves.Vertex[k1],
                    pp.Curves.Vertex[k2])
                    < d
                    * Math.ddist(pp.Curves.Vertex[k1], pp.Curves.Vertex[k2])
                    * Constant.COS179) {
                return true;
            }
        }

        /* the curve we're working in: */
        p0 = pp.Curves.ControlPoints[Math.mod(i, m)][2];
        p1 = pp.Curves.Vertex[Math.mod(i + 1, m)];
        p2 = pp.Curves.Vertex[Math.mod(j, m)];
        p3 = pp.Curves.ControlPoints[Math.mod(j, m)][2];

        /* determine its area */
        area = areac[j] - areac[i];
        area -= Math.dpara(
                pp.Curves.Vertex[0],
                pp.Curves.ControlPoints[i][2],
                pp.Curves.ControlPoints[j][2]) / 2;

        if (i >= j) {
            area += areac[m];
        }

        /* find intersection o of p0p1 and p2p3. Let t,s such that o =
           interval(t,p0,p1) = interval(s,p3,p2). Let A be the area of the
           triangle (p0,o,p3). */

        A1 = Math.dpara(p0, p1, p2);
        A2 = Math.dpara(p0, p1, p3);
        A3 = Math.dpara(p0, p2, p3);
        /* A4 = dpara(p1, p2, p3); */
        A4 = A1 + A3 - A2;

        if (A2 == A1) {  /* this should never happen */
            return true;
        }

        t = A3 / (A3 - A4);
        s = A2 / (A2 - A1);
        A = A2 * t / 2.0;

        if (A == 0.0) {  /* this should never happen */
            return true;
        }

        R = area / A;     /* relative area */
        alpha = 2 - java.lang.Math.sqrt(4 - R / 0.3);  /* overall alpha for p0-o-p3 curve */
        res.c = new DoublePoint[2];
        res.c[0] = Math.interval(t * alpha, p0, p1);
        res.c[1] = Math.interval(s * alpha, p3, p2);
        res.alpha = alpha;
        res.t = t;
        res.s = s;

        p1 = res.c[0];
        p2 = res.c[1];  /* the proposed curve is now (p0,p1,p2,p3) */

        res.pen = 0;

        /* calculate penalty */
        /* check tangency with edges */
        for (k = Math.mod(i + 1, m); k != j; k = k1) {
            k1 = Math.mod(k + 1, m);
            t = Math.Tangent(p0, p1, p2, p3, pp.Curves.Vertex[k], pp.Curves.Vertex[k1]);
            if (t < -.5) {
                return true;
            }
            pt = Math.Bezier(t, p0, p1, p2, p3);
            d = Math.ddist(pp.Curves.Vertex[k], pp.Curves.Vertex[k1]);
            if (d == 0.0) {  /* this should never happen */

                return true;
            }
            d1 = Math.dpara(pp.Curves.Vertex[k], pp.Curves.Vertex[k1], pt) / d;
            if (java.lang.Math.abs(d1) > opttolerance) {
                return true;
            }
            if (Math.iprod(pp.Curves.Vertex[k], pp.Curves.Vertex[k1], pt) < 0
                    || Math.iprod(pp.Curves.Vertex[k1], pp.Curves.Vertex[k], pt) < 0) {
                return true;
            }
            res.pen += d1 * d1;
        }

        /* check corners */
        for (k = i; k != j; k = k1) {
            k1 = Math.mod(k + 1, m);
            t = Math.Tangent(p0, p1, p2, p3, pp.Curves.ControlPoints[k][2], pp.Curves.ControlPoints[k1][2]);
            if (t < -.5) {
                return true;
            }
            pt = Math.Bezier(t, p0, p1, p2, p3);
            d = Math.ddist(pp.Curves.ControlPoints[k][2], pp.Curves.ControlPoints[k1][2]);
            if (d == 0.0) {  /* this should never happen */
                return true;
            }
            d1 = Math.dpara(pp.Curves.ControlPoints[k][2],
                    pp.Curves.ControlPoints[k1][2], pt) / d;
            d2 = Math.dpara(pp.Curves.ControlPoints[k][2],
                    pp.Curves.ControlPoints[k1][2],
                    pp.Curves.Vertex[k1]) / d;
            d2 *= 0.75 * pp.Curves.Alpha[k1];
            if (d2 < 0) {
                d1 = -d1;
                d2 = -d2;
            }
            if (d1 < d2 - opttolerance) {
                return true;
            }
            if (d1 < d2) {
                res.pen += (d1 - d2) * (d1 - d2);
            }
        }

        return false;
    }


    /**
     * Stage 2: calculate the optimal polygon (Sec. 2.2.2-2.2.4).
     * Auxiliary function: calculate the penalty of an edge from i to j in
     * the given path. This needs the "lon" and "sum*" data.
     *
     * @param pp
     * @param i
     * @param j
     * @return
     */
    private static double penalty3(Path pp, int i, int j) {

        int n = pp.pt.length;

        /* assume 0<=i<j<=n  */
        double x, y, x2, xy, y2;
        double k;
        double a, b, c, s;
        double px, py, ex, ey;
        SumStruct[] sums = pp.Sums;
        IntPoint[] pt = pp.pt;

        int r = 0; /* rotations from i to j */
        if (j >= n) {
            j -= n;
            r += 1;
        }

        x = sums[j + 1].X - sums[i].X + r * sums[n].X;
        y = sums[j + 1].Y - sums[i].Y + r * sums[n].Y;
        x2 = sums[j + 1].X2 - sums[i].X2 + r * sums[n].X2;
        xy = sums[j + 1].XY - sums[i].XY + r * sums[n].XY;
        y2 = sums[j + 1].Y2 - sums[i].Y2 + r * sums[n].Y2;
        k = j + 1 - i + r * n;

        px = (pt[i].X + pt[j].X) / 2.0 - pt[0].X;
        py = (pt[i].Y + pt[j].Y) / 2.0 - pt[0].Y;
        ey = (pt[j].X - pt[i].X);
        ex = -(pt[j].Y - pt[i].Y);

        a = ((x2 - 2 * x * px) / k + px * px);
        b = ((xy - x * py - y * px) / k + px * py);
        c = ((y2 - 2 * y * py) / k + py * py);

        s = ex * ex * a + 2 * ex * ey * b + ey * ey * c;

        return java.lang.Math.sqrt(s);
    }


    private static void PathList_to_ListOfCurveArrays(ArrayList plistp, ArrayList ListOfCurveArrays) {

    }

    /**
     * Compute a path in the binary matrix.
     * Start path at the point (x0,x1), which must be an upper left corner
     * of the path. Also compute the area enclosed by the path. Return a
     * new path_t object, or NULL on error (note that a legitimate path
     * cannot have length 0).
     * We omit turnpolicies and sign
     *
     * @param Matrix Binary Matrix
     * @param Start  start searching point
     * @return
     */
    private static Path findpath(boolean[][] Matrix, IntPoint Start) {

        ArrayList<IntPoint> L = new ArrayList<IntPoint>();

        Direction Dir = Direction.North;
        int x;
        int y;
        int area = 0;
        int diry = -1;
        x = Start.X;
        y = Start.Y;

        do {
            L.add(new IntPoint(x, y));
            int _y = y;
            int[] result = findNextTrace(Matrix, x, y, Dir);
            x = result[0];
            y = result[1];
            Dir = Direction.valueOf(Direction.getName(result[2]));
            diry = _y - y;
            area += x * diry;
        } while ((x != Start.X) || (y != Start.Y));

        if (L.size() == 0) {
            return null;
        }

        Path result = new Path();
        result.pt = new IntPoint[L.size()];
        result.area = area;

        for (int i = 0; i < L.size(); i++) {
            result.pt[i] = L.get(i);
        }

        // Shift 1 to be compatible with Potrace
        if (result.pt.length > 0) {
            IntPoint P = result.pt[result.pt.length - 1];
            for (int i = result.pt.length - 1; i >= 0; i--) {
                if (i > 0) {
                    result.pt[i] = result.pt[i - 1];
                } else {
                    result.pt[0] = P;
                }
            }
        }

        result.MonotonIntervals = GetMonotonIntervals(result.pt);

        return result;
    }


    /**
     * @param Matrix
     * @param x
     * @param y
     * @param Dir
     * @return
     */
    private static int[] findNextTrace(boolean[][] Matrix, int x, int y, Direction Dir) {

        if (Dir == Direction.West) {
            if (!Matrix[x + 1][y + 1]) {
                y++;
                Dir = Direction.North;
            } else {
                if (!Matrix[x + 1][y]) {
                    x++;
                    Dir = Direction.West;
                } else {
                    y--;
                    Dir = Direction.South;
                }
            }
        } else if (Dir == Direction.South) {
            if (!Matrix[x + 1][y]) {
                x++;
                Dir = Direction.West;
            } else {
                if (!Matrix[x][y]) {
                    y--;
                    Dir = Direction.South;
                } else {
                    x--;
                    Dir = Direction.East;
                }
            }
        } else if (Dir == Direction.East) {
            if (!Matrix[x][y]) {
                y--;
                Dir = Direction.South;
            } else {
                if (!Matrix[x][y + 1]) {
                    x--;
                    Dir = Direction.East;
                } else {
                    y++;
                    Dir = Direction.North;
                }
            }
        } else if (Dir == Direction.North) {
            if (!Matrix[x][y + 1]) {
                x--;
                Dir = Direction.East;
            } else {
                if (!Matrix[x + 1][y + 1]) {
                    y++;
                    Dir = Direction.North;
                } else {
                    x++;
                    Dir = Direction.West;
                }
            }
        }

        int[] result = new int[3];
        result[0] = x;
        result[1] = y;
        result[2] = Dir.getValue();

        return result;
    }


    /**
     * Divide contour into several segments according to segment direction
     *
     * @param Pts
     * @return
     */
    private static ArrayList<MonotonInterval> GetMonotonIntervals(IntPoint[] Pts) {

        ArrayList<MonotonInterval> result = new ArrayList<MonotonInterval>();
        int n = Pts.length;
        if (n == 0) {
            return result;
        }


        //Step 1: Divide contour into several segments
        ArrayList<MonotonInterval> L = new ArrayList<MonotonInterval>();

        //----- Start with Strong Monoton (Pts[i].y < Pts[i+1].y) or (Pts[i].y > Pts[i+1].y)
        int FirstStrongMonoton = 0;
        while (Pts[FirstStrongMonoton].Y == Pts[FirstStrongMonoton + 1].Y) {
            FirstStrongMonoton++;
        }
        boolean Up = (Pts[FirstStrongMonoton].Y < Pts[FirstStrongMonoton + 1].Y);
        MonotonInterval Interval = new MonotonInterval(Up, FirstStrongMonoton, FirstStrongMonoton);
        L.add(Interval);

        int i = FirstStrongMonoton;
        do {
            // Interval.to = i;
            if (Pts[i].Y == Pts[Math.mod(i + 1, n)].Y
                    || Up == (Pts[i].Y < Pts[Math.mod(i + 1, n)].Y)) {
                Interval.to = i;
            } else {
                Up = (Pts[i].Y < Pts[Math.mod(i + 1, n)].Y);
                Interval = new MonotonInterval(Up, i, i);
                L.add(Interval);
            }
            i = Math.mod(i + 1, n);
        } while (i != FirstStrongMonoton);


        //Step 2: Make the number of segments even
        if (L.size() / 2 * 2 != L.size()) {// Connect the Last with first
            MonotonInterval M0 = L.get(0);
            MonotonInterval ML = L.get(L.size() - 1);
            M0.from = ML.from;
            L.remove(L.size() - 1);
        }


        //Step 3: Order the segments by y-value(first-key) and x-value(second-key)
        //where the segment with lowest x-value and y-value ranks first

        //----- order now by the min y - value of interval to result
        // and as second Key by the x-value
        //
        while (L.size() > 0) {
            MonotonInterval M = L.get(0);
            i = 0;
            // order by y-value
            while (i < result.size()
                    && Pts[M.Min()].Y > Pts[(result.get(i)).Min()].Y) {
                i++;
            }
            // order by x- value as second Key
            while (i < result.size()
                    && Pts[M.Min()].Y == Pts[(result.get(i)).Min()].Y
                    && (Pts[M.Min()].X > (Pts[(result.get(i)).Min()].X))) {
                i++;
            }
            result.add(i, M);
            L.remove(0);
        }
        return result;
    }

    /**
     * Invert color with path
     *
     * @param Matrix
     * @param P
     */
    private static void Xor_Path(boolean[][] Matrix, Path P) {

        int i = 0;
        int n = P.pt.length;

        ArrayList<MonotonInterval> MonotonIntervals = P.MonotonIntervals;
        if (MonotonIntervals.size() == 0) {
            return;
        }
        MonotonInterval MI = MonotonIntervals.get(0);
        MI.ResetCurrentID(n);

        int y = P.pt[MI.CurrentID].Y;
        ArrayList<MonotonInterval> CurrentIntervals = new ArrayList<MonotonInterval>();
        CurrentIntervals.add(MI);
        MI.CurrentID = MI.Min();

        while (i + 1 < MonotonIntervals.size()
                && MonotonIntervals.get(i + 1).MinY(P.pt) == y) {
            MI = MonotonIntervals.get(i + 1);
            MI.ResetCurrentID(n);
            CurrentIntervals.add(MI);
            i++;
        }

        while (CurrentIntervals.size() > 0) {   // invertLine

            for (int k = 0; k < CurrentIntervals.size() - 1; k++) {

                int x1 = P.pt[(CurrentIntervals.get(k)).CurrentID].X + 1;
                int x2 = P.pt[(CurrentIntervals.get(k + 1)).CurrentID].X;

                for (int x = x1; x <= x2; x++) {
                    Matrix[x][y] = !Matrix[x][y];
                }
                k++;
            }

            y++;
            for (int j = CurrentIntervals.size() - 1; j >= 0; j--) {

                MonotonInterval M = CurrentIntervals.get(j);

                if (y > M.MaxY(P.pt)) {
                    CurrentIntervals.remove(j);
                    continue;
                }
                int CID = M.CurrentID;
                do {
                    if (M.Increasing)
                        CID = Math.mod(CID + 1, n);
                    else
                        CID = Math.mod(CID - 1, n);
                } while (P.pt[CID].Y < y);

                M.CurrentID = CID;
            }

            // Add Items of MonotonIntervals with Down.y==y
            while (i + 1 < MonotonIntervals.size()
                    && (MonotonIntervals.get(i + 1)).MinY(P.pt) == y) {

                MonotonInterval NewInt = MonotonIntervals.get(i + 1);
                int j = 0;

                // search the correct x-Position
                int _x = P.pt[NewInt.Min()].X;
                while (j < CurrentIntervals.size()
                        && _x > P.pt[(CurrentIntervals.get(j)).CurrentID].X) {
                    j++;
                }
                CurrentIntervals.add(j, NewInt);
                NewInt.ResetCurrentID(n);
                i++;
            }
        }
    }

    //Below is for test
    public static void main(String[] args) {

        String dllPath = "C:\\OpenCV\\opencv\\build\\java\\x64\\opencv_java320.dll";
        System.load(dllPath);

        testBm_to_pathlist();
    }

    private static void testBitmapToBinary() {

        String filePath = "E:\\Java_Projects\\Potrace\\resources\\sourceEntireImages\\4a.png";
        Mat srcImage = Imgcodecs.imread(filePath);
        boolean[][] result = bitmapToBinary(srcImage);

        outputMatrix(result);
    }

    private static void testBm_to_pathlist() {

        String filePath = "E:\\Java_Projects\\Potrace\\resources\\sourceEntireImages\\11a.png";
        Mat srcImage = Imgcodecs.imread(filePath);


        boolean[][] matrix = bitmapToBinary(srcImage);

        outputMatrix(matrix);

        ArrayList<ArrayList<Path>> ListOfCurveArray = new ArrayList<>();
        bm_to_pathlist(matrix, ListOfCurveArray);

        System.out.println();
    }

    /**
     * Output the given binary matrix in 0/1 form
     *
     * @param matrix
     */
    private static void outputMatrix(boolean[][] matrix) {

        for (int row = 0; row < matrix.length; row++) {
            for (int col = 0; col < matrix[row].length; col++) {

                if (matrix[row][col]) {
                    System.out.print(0);
                } else {
                    System.out.print(1);
                }
            }
            System.out.println();
        }

        System.out.println();
        System.out.println();
    }
}
