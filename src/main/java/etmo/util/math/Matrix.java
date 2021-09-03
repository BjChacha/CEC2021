package etmo.util.math;

import java.util.ArrayList;

public class Matrix {
    public static double[][] matMul(double[][] mat1, double[][] mat2) {
        // mat1: d1 x d2
        // mat2: d2 x d3
        // output: d1 x d3
        int d1 = mat1.length;
        int d2 = mat1[0].length;
        int d3 = mat2[0].length;
        double[][] output = new double[d1][d3];

        return output;
    }

    public static double[][] matT(double[][] mat) {
        // mat: d1 x d2
        // output: d2 x d1
        int d1 = mat.length;
        int d2 = mat[0].length;
        double[][] output = new double[d1][d2];

        return output;
    }

    public static double[] getMeanOfMat(ArrayList<double[]> mat) {
        double[][] newMat = new double[mat.size()][];
        for (int i = 0; i < newMat.length; i++) {
            newMat[i] = mat.get(i);
        }
        return getMeanOfMat(newMat);
    }

    public static double[] getMeanOfMat(double[][] mat) {
        int n = mat.length;
        int d = mat[0].length;
        double[] mean = new double[d];
        for (int j = 0; j < d; j++) {
            for (int i = 0; i < n; i++) {
                mean[j] += mat[i][j];
            }
            mean[j] /= n;
        }
        return mean;
    }

    public static double[] getStdOfMat(ArrayList<double[]> mat) {
        double[][] newMat = new double[mat.size()][];
        for (int i = 0; i < newMat.length; i++) {
            newMat[i] = mat.get(i);
        }
        return getStdOfMat(newMat);
    }

    public static double[] getStdOfMat(double[][] mat) {
        int n = mat.length;
        int d = mat[0].length;
        double[] std = new double[d];
        double[] mean = getMeanOfMat(mat);
        for (int j = 0; j < d; j++) {
            for (int i = 0; i < n; i++) {
                std[j] += Math.pow(mat[i][j] - mean[j], 2);
            }
            std[j] = Math.sqrt(std[j] / Math.max(n - 1, 1));
        }
        return std;
    }
}
