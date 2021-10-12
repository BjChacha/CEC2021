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
        return Math.pow(distance, 1.0 / d);
    }

    public static double getCosineSimilarity(double[] v1, double[] v2) {
        double similarity = 0;
        similarity = Vector.vecDot(v1, v2) / (Vector.vecModule(v1) * Vector.vecModule(v2));
        return similarity;
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

    public static double getCoralLoss(double[][] p1, double[][] p2) {
        int d = p1[0].length; 
        double dist = 0;
        double[][] sigma1 = Matrix.getMatSigma(p1);
        double[][] sigma2 = Matrix.getMatSigma(p2);

        for (int i = 0; i < d; i ++) {
            for (int j = 0; j < d; j ++) {
                dist += Math.pow(sigma1[i][j] - sigma2[i][j], 2);
            }
        }
        return Math.sqrt(dist) / (4.0 * Math.pow(d, 2));
    }

    public static double getCoralLossWithSigma(double[][] sigma1, double[][] sigma2) {
        int d = sigma1.length; 
        double dist = 0;
        for (int i = 0; i < d; i ++) {
            for (int j = 0; j < d; j ++) {
                dist += Math.pow(sigma1[i][j] - sigma2[i][j], 2);
            }
        }
        return Math.sqrt(dist) / (4.0 * Math.pow(d, 2));
    }

    public static double getCorrelationMatrixDistance(double[][] sigma1, double[][] sigma2) {
        double dist = 0;
        dist = Matrix.matTrace(Matrix.matMul(sigma1, sigma2)) / Matrix.matNorm(sigma1) / Matrix.matNorm(sigma2);
        return 1 - dist;
    }
}
