package potrace;

class Math {

    /* integer arithmetic */
    public static int sign(int x) {
        return Integer.compare(x, 0);
    }

    public static int sign(double x) {
        return ((x) > 0 ? 1 : (x) < 0 ? -1 : 0);
    }

    public static int abs(int a) {
        return ((a) > 0 ? (a) : -(a));
    }

    public static int min(int a, int b) {
        return ((a) < (b) ? (a) : (b));
    }

    public static int max(int a, int b) {
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

    static int floordiv(int a, int n) {
        return a >= 0 ? a / n : -1 - (-1 - a) / n;
    }
}
