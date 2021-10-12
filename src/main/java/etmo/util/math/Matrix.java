package etmo.util.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

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

    public static double[][] matSub(double[][] mat1, double[][] mat2) {
        assert mat1.length == mat2.length && mat1[0].length == mat2[0].length;
        int d1 = mat1.length;
        int d2 = mat1[0].length;
        double[][] output = new double[d1][d2];

        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                output[i][j] = mat1[i][j] - mat2[i][j];
            }
        }

        return output;
    }


    public static double[][] matElemMul(double[][] mat, double factor) {
        double[][] output = mat.clone();
        Arrays.setAll(output, i -> Vector.vecElemMul(mat[i], factor));
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


    public static double[][] matSlice(double[][] mat, int rowStart, int rowEnd, int colStart, int colEnd) {
        assert rowStart >= 0 && rowEnd < mat.length && rowStart < rowEnd 
            && colStart >= 0 && colEnd < mat.length && colStart < colEnd;
        double[][] newMat = new double[rowEnd-rowStart+1][colEnd-colStart+1];
        for (int i = rowStart; i <= rowEnd; i++) {
            for (int j = colStart; j <= colEnd; j++) {
                newMat[i-rowStart][j-colStart] = mat[i][j];
            }
        }
        return newMat;
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


    public static void matAdd_(double[][] mat, double e) {
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[i].length; j++) {
                mat[i][j] += e;
            }
        }
    }

    public static void matMul_(double[][] mat, double e) {
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[i].length; j++) {
                mat[i][j] *= e;
            }
        }
    }

    public static double[][] getRandomMat(int d1, int d2) {
        assert d1 > 0 && d2 > 0;
        double[][] mat = new double[d1][d2];
        Random rand = new Random();
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                mat[i][j] = rand.nextDouble();
            }
        }
        return mat;
    }

    public static double[][] getMatSigma(double[][] mat) {
        assert mat.length > 1 && mat[0].length > 1;
        int m = mat.length;
        int n = mat[0].length;
        double[][] sigma = new double[n][n];

        double[] mean = new double[n];
        Arrays.setAll(mean, i -> { 
            double sum = 0;
            for (int j = 0; j < m; j ++)
                sum += mat[j][i];
            return sum / m;
        });

        double[][] X = mat.clone();
        Arrays.setAll(X, i -> Vector.vecSub(X[i], mean));

        sigma = matElemMul(matMul(matTranspose(X), X), 1.0 / (m - 1));

        return sigma;
    }

    public static double[][] matSqrt(double[][] mat) {
        double[][] output = null;
        RealMatrix m = new Array2DRowRealMatrix(mat);
        SingularValueDecomposition svd = new SingularValueDecomposition(m);

        double[][] U = svd.getU().getData();
        double[][] S = svd.getS().getData();
        double[][] VT = svd.getVT().getData();
        for (int i = 0; i < S.length; i ++) {
            S[i][i] = Math.sqrt(S[i][i]);
        }
        output = matMul(matMul(U, S), VT);

        return output;
    }

    public static double matTrace(double[][] mat) {
        assert mat.length > 0 && mat.length == mat[0].length;
        double trace = 0;
        for (int i = 0; i < mat.length; i ++) {
            trace += mat[i][i];
        }
        return trace;
    }

    public static double matNorm(double[][] mat) {
        double res = 0;
        for (int i = 0; i < mat.length; i ++) {
            for (int j = 0; j < mat[i].length; j ++) {
                res += Math.pow(mat[i][j], 2);
            }
        }
        return res;
    }

    public static boolean matIsNonSingular(double[][] mat) {
        RealMatrix m = new Array2DRowRealMatrix(mat);
        SingularValueDecomposition svd = new SingularValueDecomposition(m);
        DecompositionSolver solver = svd.getSolver();
        return solver.isNonSingular();
    }
}
