package potrace;


/**
 * a type for the result of opti_penalty
 */
class opti {

    double pen;       /* penalty */
    DoublePoint[] c;   /* curve parameters */
    double t, s;       /* curve parameters */
    double alpha;       /* curve parameter */
}
