package etmo.metaheuristics.matbml2;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.dimensionalityreduction.PCA;
import org.nd4j.linalg.ops.transforms.Transforms;

class WassersteinDistance {
    public static double getWD(double[][] p1, double[][] p2){
        double distance = 0;

        double[][] d1 = getDistribution(p1);
        double[][] d2 = getDistribution(p2);

        int d = Math.min(d1[0].length, d2[0].length);
        double dis1 = 0;
        double dis2 = 0;

        for (int i = 0; i < d; i++){
            dis1 += Math.pow(d1[0][i] - d2[0][i], 2);
            dis2 += Math.pow(d2[1][i] - d2[1][i], 2);
        }

        distance = Math.sqrt(dis1 + dis2);
        return distance;
    }

    private static double[][] getDistribution(double[][] p1) {
        double[][] res = new double[2][];
        int n = p1.length;
        if (n < 1) {
            System.out.println("Get Distribution Error: population length is 0.");
            return res;
        }
        int d = p1[0].length;

        // mean
        res[0] = new double[d];
        for (int i = 0; i < d; i++){
            double sum = 0;
            for (int j = 0; j < n; j++)
                sum += p1[j][i];
            res[0][i] = (sum / n);
        }

        // std
        res[1] = new double[d];
        for (int i = 0; i < d; i++){
            double sum = 0;
            for (int j = 0; j < n; j++)
                sum += Math.pow(p1[j][i] - res[0][i], 2);
            res[1][i] = Math.sqrt(sum / (n - 1));
        }

        return res;
    }

    public static double getWD2(double[][] p1, double[][] p2){
        int size = p1.length;
        double[][] distanceMat = new double[size][size];
        for (int i = 0; i < size - 1; i++){
            distanceMat[i][i] = Double.MAX_VALUE;
            for (int j = i + 1; j < size; j++){
                distanceMat[i][j] = distanceMat[j][i] = Utils.distVector(p1[i], p2[j]);
            }
        }

        int[] assignment = new HungarianAlgorithm(distanceMat).execute();
        double distance = 0;
        for (int i = 0; i < size; i++){
            distance += distanceMat[i][assignment[i]];
        }

        return distance / size;
    }
}
