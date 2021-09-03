package etmo.metaheuristics.ematomkt;

import java.util.ArrayList;
import java.util.Arrays;

import etmo.core.SolutionSet;
import etmo.util.JMException;
import etmo.util.math.Distance;

public class Utils {
    public static double MMD(double[][] mat1, double[][] mat2, double sigma) {
        double mmd = 0;
        int d = mat1[0].length;

        double[][] K = rbfDot(mat1, mat1, sigma);
        double[][] L = rbfDot(mat2, mat2, sigma);
        double[][] KL = rbfDot(mat1, mat2, sigma);

        for (int i = 0; i < d; i++) {
            double sumK = 0;
            double sumL = 0;
            double sumKL = 0;
            for (int j = 0; j < d; j++) {
                sumK += K[i][j];
                sumL += L[i][j];
                sumKL += KL[i][j];
            }
            mmd += (sumK / Math.pow(d, 2) + sumL / Math.pow(d, 2) - 2 * sumKL / Math.pow(d, 2));
        }

        mmd = Math.sqrt(mmd);
        return mmd;
    }

    public static double[][] rbfDot(double[][] mat1, double[][] mat2, double sigma) {
        int d = mat1[0].length;
        int n = mat1.length;
        double[] squareSum1 = new double[d];
        double[] squareSum2 = new double[d];
        double[][] crossMat = new double[d][d];
        double[][] resMat = new double[d][d];

        for (int j = 0; j < d; j++){
            for (int i = 0; i < n; i++) {
                squareSum1[j] += Math.pow(mat1[i][j], 2);
                squareSum2[j] += Math.pow(mat2[i][j], 2);
            }

            for (int i = 0; i < d; i++) {
                for (int k = 0; k < n; k++) {
                    crossMat[j][i] += mat1[k][j] * mat2[k][i];
                }
            }
        }

        for (int i = 0; i < d; i++){
            for (int j = 0; j < d; j++){
                resMat[i][j] = squareSum1[i] + squareSum2[j] - 2 * crossMat[i][j];
                resMat[i][j] = Math.exp(-resMat[i][j] / 2 / Math.pow(sigma, 2));
            }
        }

        return resMat;
    }

    // public static double[][][] kmeans2(SolutionSet union, int clusterNum, int populationSize) {

    // }

    public static double[][][] kmeans(SolutionSet union, int clusterNum, int populationSize) throws JMException {
        double[][][] clusters = new double[clusterNum][][];

        // 随机设置质点
        int[] perm = randomPermutation(union.size(), clusterNum);
        double[][] centroids = new double[clusterNum][];
        ArrayList<ArrayList<double[]>> tmpClusters = new ArrayList<>();
        for (int i = 0; i < clusterNum; i++) {
            tmpClusters.add(new ArrayList<double[]>());
            centroids[i] = union.get(perm[i]).getDecisionVariablesInDouble();
        }


        boolean clusterChanged = true;
        while (clusterChanged) {
            // 1. 每个点选择最近的质点
            for (int i = 0; i < union.size(); i++){
                double[] selectedPoint = union.get(i).getDecisionVariablesInDouble();

                int selectedClusterIndex = -1;
                double minDistance = Double.MAX_VALUE;
                for (int j = 0; j < clusterNum; j++) {
                    double distance = Distance.getDistance(selectedPoint, centroids[j], 1);
                    if (distance < minDistance) {
                        minDistance = distance;
                        selectedClusterIndex = j;
                    }
                }
                tmpClusters.get(selectedClusterIndex).add(selectedPoint);
                if (i < populationSize)
                    union.get(i).setFlag(selectedClusterIndex);
                
            }

            // 2. 更新每个簇的质点
            clusterChanged = false;
            for (int i = 0; i < clusterNum; i++) {
                int n = tmpClusters.get(i).size();
                int d = centroids[i].length;
                double[] newCentroid = new double[d];
                for (int j = 0; j < d; j++) {
                    for (int k = 0; k < n; k++) {
                        newCentroid[j] += tmpClusters.get(i).get(k)[j];
                    }
                    newCentroid[j] /= n;

                    if (Math.abs(newCentroid[j] - centroids[i][j]) > 1e-10)
                        clusterChanged = true;
                    centroids[i][j] = newCentroid[j];
                }
            }

            // 3. 旧簇清空
            if (clusterChanged){
                for (int i = 0; i < clusterNum; i++) {
                    tmpClusters.get(i).clear();
                }
            }
            else {
                // 4. 若存在空簇，则重新初始化
                boolean isEmpty = false;
                for (int i = 0; i < clusterNum; i++) {
                    if (tmpClusters.get(i).size() == 0){
                        isEmpty = true;
                        break;
                    }
                }
                if (isEmpty){
                    System.out.println("Empty cluster detected.");
                    clusterChanged = true;
                    int[] tperm = randomPermutation(union.size(), clusterNum);
                    for (int i = 0; i < clusterNum; i++) {
                        centroids[i] = union.get(tperm[i]).getDecisionVariablesInDouble();
                    }
                    for (int i = 0; i < clusterNum; i++) {
                        tmpClusters.get(i).clear();
                    }
                }
            }
        }

        for (int i = 0; i < clusterNum; i++) {
            clusters[i] = tmpClusters.get(i).toArray(new double[tmpClusters.get(i).size()][]);
        }

        return clusters;
    }

    public static int[] randomPermutation(int range, int size) {
        /*
        range: generate array with value in [0,range)
        */
        if (range < size) {
            System.out.println("Error: range cannot less than size. setting size = range...");
            size = range;
        }

        int[] arr = new int[range];
        for (int i = 0; i < range; i++)
            arr[i] = i;

        for (int i = 0; i < range - 1; i++) {
            int j = etmo.util.PseudoRandom.randInt(i, range - 1);
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }

        return Arrays.copyOfRange(arr, 0, size);
	}

    public static int[] randomPermutation(int range) {
        return randomPermutation(range, range);
    }
}