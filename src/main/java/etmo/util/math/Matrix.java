package etmo.util.math;

import java.util.ArrayList;
import java.util.Random;

public class Matrix {
    public static double[][] matMul(double[][] mat1, double[][] mat2) {
        // mat1: d1 x d2
        // mat2: d2 x d3
        // output: d1 x d3
        assert mat1[0].length == mat2.length;

        int d1 = mat1.length;
        int d2 = mat1[0].length;
        int d3 = mat2[0].length;
        double[][] output = new double[d1][d3];
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d3; j++) {
                for (int k = 0; k < d2; k++) {
                    output[i][j] += mat1[i][k] * mat2[k][j];
                }
            }
        }

        return output;
    }

    public static double[][] matAdd(double[][] mat1, double[][] mat2) {
        assert mat1.length == mat2.length && mat1[0].length == mat2[0].length;

        int d1 = mat1.length;
        int d2 = mat1[0].length;
        double[][] output = new double[d1][d2];

        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                output[i][j] = mat1[i][j] + mat2[i][j];
            }
        }

        return output;
    }

    public static double[][] matTranspose(double[][] mat) {
        // mat: d1 x d2
        // output: d2 x d1
        int d1 = mat.length;
        int d2 = mat[0].length;
        double[][] output = new double[d2][d1];
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                output[j][i] = mat[i][j];
            }
        }
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

    public static double[] randomUnitVector(int d) {
        double[] vector = new double[d];
        double sum = 0;
        Random r = new Random();
        for (int i = 0; i < d; i ++) {
            vector[i] = r.nextGaussian();
            sum += Math.pow(vector[i], 2);
        }
        sum = Math.sqrt(sum);
        for (int i = 0; i < d; i ++) {
            vector[i] /= sum;
        }
        return vector;
    }

    public static double[] unifyVector(double[] vector) {
        double sum = 0;
        double[] output = new double[vector.length];
        for (double e: vector) {
            sum += Math.pow(e, 2);
        }
        sum = Math.sqrt(sum);
        for (int i = 0; i < output.length; i++) {
            output[i] = vector[i] / sum;
        }
        return output;
    }

    public static double[] normailizeVector(double[] vector) {
        double[] output = new double[vector.length];
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        for (int i = 0; i < vector.length; i++) {
            min = Math.min(min, vector[i]);
            max = Math.max(max, vector[i]);
        }
        for (int i = 0; i < vector.length; i++) {
            output[i] = (vector[i] - min) / (max - min);
        }
        return output;
    }

    public static void eAdd_(double[][] mat, double e) {
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[i].length; j++) {
                mat[i][j] += e;
            }
        }
    }

    public static void eMul_(double[][] mat, double e) {
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[i].length; j++) {
                mat[i][j] *= e;
            }
        }
    }
}
