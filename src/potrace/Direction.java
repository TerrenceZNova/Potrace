package potrace;

enum Direction {

    North(1),
    East(2),
    South(3),
    West(4);

    private int value;

    private Direction(int value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }

    public static String getName(int value) {
        switch (value) {
            case 1: {
                return "North";
            }
            case 2: {
                return "East";
            }
            case 3: {
                return "South";
            }
            case 4: {
                return "West";
            }
            default: {
                return null;
            }
        }
    }
}
