package potrace;

public class Constant {

    //----------------------Potrace Constants and aux functions
    final int POTRACE_CORNER = 1;
    final int POTRACE_CURVETO = 2;
    final static double COS179 = java.lang.Math.cos(179 * java.lang.Math.PI / 180);

    //Area of largest path to be ignored.
    final static int turdsize = 2;

    /// Corner threshold.
    final static double alphamax = 1.0;

    /// Use curve optimization.
    /// optimize the path p, replacing sequences of Bezier segments by a single segment when possible.
    final static boolean curveoptimizing = true;

    //Curve optimization tolerance
    final static double opttolerance = 0.2;
}
