package etmo.util.math;

public class Distance {
    public static double getDistance(double[] p1, double[] p2) {
        return getDistance(p1, p2, 2);
    }

    public static double getDistance(double[] p1, double[] p2, int d) {
        double distance = 0;
        for (int i = 0; i < p1.length; i++) {
            distance += Math.abs(Math.pow(p1[i] - p2[i], d));
        }
        return Math.pow(distance, 1/d);
    }
}
