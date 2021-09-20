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

    public static double getWassersteinDistance(double[][] p1, double[][] p2) {
        double distance = 0;
        int size = p1.length;
        double[][] distanceMat = new double[size][size];
        for (int i = 0; i < size - 1; i++){
            distanceMat[i][i] = Double.MAX_VALUE;
            for (int j = i + 1; j < size; j++){
                distanceMat[i][j] = distanceMat[j][i] = getDistance(p1[i], p2[j]);
            }
        }

        int[] assignment = new HungarianAlgorithm(distanceMat).execute();
        for (int i = 0; i < size; i++){
            distance += distanceMat[i][assignment[i]];
        }

        return distance / size;
    }
}
