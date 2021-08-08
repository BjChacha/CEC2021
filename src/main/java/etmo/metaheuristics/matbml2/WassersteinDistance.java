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
        INDArray P1 = new NDArray(p1);
        INDArray P2 = new NDArray(p2);

        INDArray reduced1 = PCA.pca(P1, 5, false);
        INDArray reduced2 = PCA.pca(P2, 5, false);

        INDArray[] d1 = PCA.covarianceMatrix(reduced1);
        INDArray[] d2 = PCA.covarianceMatrix(reduced2);

        INDArray mean1 = d1[1];
        INDArray cov1 = d1[0];
        INDArray mean2 = d2[1];
        INDArray cov2 = d2[0];

        double distance1 = 0;
        double distance2 = 0;

        distance1 = Transforms.pow(mean1.sub(mean2), 2).sumNumber().doubleValue();
        INDArray tmp = Transforms.pow(cov1, 0.5);
        distance2 = (Transforms.pow((Transforms.pow(cov1, 0.5).sub(Transforms.pow(cov2, 0.5))), 2)).sumNumber().doubleValue();

        double distance = distance1 + distance2;
        return distance;
    }
}
