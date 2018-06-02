package potrace;

class MonotonInterval {

    public boolean Increasing;
    public int from;
    public int to;

    public void ResetCurrentID(int modulo) {
        if (!Increasing)
            CurrentID = Math.mod(Min() + 1, modulo);
        else
            CurrentID = Min();
    }

    public int CurrentID; // only used by Invert

    public MonotonInterval(boolean Increasing, int from, int to) {
        this.Increasing = Increasing;
        this.from = from;
        this.to = to;
    }

    public int Min() {
        if (Increasing) return from;
        return to;
    }

    public int MinY(IntPoint[] Pts) {
        return Pts[Min()].Y;
    }

    public int MaxY(IntPoint[] Pts) {
        return Pts[Max()].Y;
    }

    public int Max() {
        if (!Increasing) return from;
        return to;
    }
}
