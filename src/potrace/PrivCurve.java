package potrace;

public class PrivCurve {


    //Number of segments.
    int Count;

    //Tag[n]: POTRACE_CORNER or POTRACE_CURVETO.
    int[] Tag;

    //c[n][i]: control points.
    //c[n][0] is unused for tag[n]=POTRACE_CORNER
    DoublePoint[][] ControlPoints;


    //for POTRACE_CORNER, this equals c[1](*c)[3];
    //c[n][i]: control points.
    //c[n][0] is unused for tag[n]=POTRACE_CORNER
    DoublePoint[] Vertex;


    //only for POTRACE_CURVETO.
    double[] Alpha;


    //"uncropped" alpha parameter - for debug output only.
    double[] Alpha0;
    double[] Beta;

    /**
     * Constructor
     *
     * @param count
     */
    public PrivCurve(int count) {

        this.Count = count;
        this.Tag = new int[Count];
        this.ControlPoints = new DoublePoint[Count][3];
        this.Vertex = new DoublePoint[Count];
        this.Alpha = new double[Count];
        this.Alpha0 = new double[Count];
        this.Beta = new double[Count];
    }
}
