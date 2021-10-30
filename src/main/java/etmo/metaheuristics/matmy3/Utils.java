//  Utils.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package etmo.metaheuristics.matmy3;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.nd4j.linalg.cpu.nativecpu.NDArray;

import etmo.util.PseudoRandom;
import etmo.util.math.Vector;

/**
 * Utilities methods to used by MOEA/D
 */
public class Utils {

	public static double distVector(double[] vector1, double[] vector2) {
		int dim = vector1.length;
		double sum = 0;
		for (int n = 0; n < dim; n++) {
			sum += Math.pow(vector1[n] - vector2[n], 2);
		}
		return Math.sqrt(sum);
	} // distVector

	public static void minFastSort(double x[], int idx[], int n, int m) {
		for (int i = 0; i < m; i++) {
			for (int j = i + 1; j < n; j++) {
				if (x[i] > x[j]) {
					double temp = x[i];
					x[i] = x[j];
					x[j] = temp;
					int id = idx[i];
					idx[i] = idx[j];
					idx[j] = id;
				} // if
			}
		} // for

	} // minFastSort
	
	public static void maxFastSort(double x[], int idx[], int n, int m) {
		for (int i = 0; i < m; i++) {
			for (int j = i + 1; j < n; j++) {
				if (x[i] < x[j]) {
					double temp = x[i];
					x[i] = x[j];
					x[j] = temp;
					int id = idx[i];
					idx[i] = idx[j];
					idx[j] = id;
				} // if
			}
		} // for

	} // minFastSort

	public static void randomPermutation(int[] perm, int size) {
		int[] index = new int[size];
		boolean[] flag = new boolean[size];

		for (int n = 0; n < size; n++) {
			index[n] = n;
			flag[n] = true;
		}

		int num = 0;
		while (num < size) {
			int start = etmo.util.PseudoRandom.randInt(0, size - 1);
			// int start = int(size*nd_uni(&rnd_uni_init));
			while (true) {
				if (flag[start]) {
					perm[num] = index[start];
					flag[start] = false;
					num++;
					break;
				}
				if (start == (size - 1)) {
					start = 0;
				} else {
					start++;
				}
			}
		} // while
	} // randomPermutation

	static void QuickSort(double[] array, int[] idx, int from, int to) {
		if (from < to) {
			double temp = array[to];
			int tempIdx = idx[to];
			int i = from - 1;
			for (int j = from; j < to; j++) {
				if (array[j] <= temp) {
					i++;
					double tempValue = array[j];
					array[j] = array[i];
					array[i] = tempValue;
					int tempIndex = idx[j];
					idx[j] = idx[i];
					idx[i] = tempIndex;
				}
			}
			array[to] = array[i + 1];
			array[i + 1] = temp;
			idx[to] = idx[i + 1];
			idx[i + 1] = tempIdx;
			QuickSort(array, idx, from, i);
			QuickSort(array, idx, i + 1, to);
		}
	}

	public static int roulette(double[] arr){
		double[] sm = softmax(arr);

		double s = 0;
		double p = PseudoRandom.randDouble();
		int idx;

		for (idx = 0; idx < sm.length - 1; idx++) {
			s += sm[idx];
			if (s >= p) break;
		}
		return idx;
	}

	public static int rouletteExceptZero(double[] arr){
		double[] sm = softmaxExceptZero(arr);

		double s = 0;
		double p = PseudoRandom.randDouble();
		int idx;

		for (idx = 0; idx < sm.length - 1; idx++) {
			s += sm[idx];
			if (s >= p) break;
		}
		return idx;
	}

	public static double[] softmax(double[] arr){
		double[] res = new double[arr.length];
		double sum = 0;

		for (int i = 0; i < res.length; i++){
			res[i] = Math.exp(arr[i]);
			sum += res[i];
		}

		for (int i = 0; i < res.length; i++) {
			res[i] /= sum;
		}
		return res;
	}

	public static double[] softmaxExceptZero(double[] arr){
		double[] res = new double[arr.length];
		double sum = 0;

		for (int i = 0; i < res.length; i++){
			if (arr[i] > 0) {
				res[i] = Math.exp(arr[i]);
				sum += res[i];
			}
		}

		for (int i = 0; i < res.length; i++) {
			res[i] /= sum;
		}
		return res;
	}

    public static double[] subspaceAlignment(double[][] mat1, double[][] mat2, double[] input){
        // 迁移方向: mat1->mat2
        RealMatrix realMatrix1 = MatrixUtils.createRealMatrix(new double[][] {input});

        RealMatrix baseVectors1 = PCAGetBaseVectors(mat1);
        RealMatrix baseVectors2 = PCAGetBaseVectors(mat2, baseVectors1.getColumnDimension());

        RealMatrix M = baseVectors1.transpose().multiply(baseVectors2);

        RealMatrix resMatrix = realMatrix1.multiply(baseVectors1);
        resMatrix = resMatrix.multiply(M);
        resMatrix = resMatrix.multiply(baseVectors2.transpose());

        double[] res = resMatrix.getData()[0];
        Vector.vecClip_(res, 0, 1);
        return res;
    }

	// PCA
    public static RealMatrix PCAGetBaseVectors(double[][] mat){
        // 归一化+缩放，使均值为零，且最大最小值之差为1.
        MatNormalize(mat);
        // 获得协方差矩阵
        RealMatrix covarianceMatrix = new Covariance(MatrixUtils.createRealMatrix(mat)).getCovarianceMatrix();
        
        // 奇异值分解
        SingularValueDecomposition svd = new SingularValueDecomposition(covarianceMatrix);
        // 取前95%信息的奇异值数量
        double[] singleValues = svd.getSingularValues();
        double sum = 0;
        for (double value: singleValues)
            sum += value;
        int d = 0;
        double curSum = 0;
        for (double value: singleValues){
            curSum += value;
            d ++;
            if (curSum >= sum * 0.95)
                break;
        }

        RealMatrix U = svd.getU();
        RealMatrix baseVectors = U.getSubMatrix(0,U.getRowDimension() - 1,0,d - 1);

        // 输入：n x D   输出：d x D   降维：D -> d
        return baseVectors;
    }

	public static RealMatrix PCAGetBaseVectors(double[][] mat, int d){
        // 归一化+缩放，使均值为零，且最大最小值之差为1.
        MatNormalize(mat);
        // 获得协方差矩阵
        RealMatrix covarianceMatrix = new Covariance(MatrixUtils.createRealMatrix(mat)).getCovarianceMatrix();

        // 奇异值分解
        SingularValueDecomposition svd = new SingularValueDecomposition(covarianceMatrix);
        RealMatrix U = svd.getU();
        RealMatrix baseVectors = U.getSubMatrix(0,U.getRowDimension() - 1,0,d - 1);

        // 输入：n x D   输出：d x D   降维：D -> d
        return baseVectors;
    }

	static void MatNormalize(double[][] mat){
        int row = mat.length;
        int col = mat[0].length;

        double[] means = new double[col];
        double[] max = new double[col];
        double[] min = new double[col];
        for (int i = 0; i < col; i++) {
            float sum = 0;
            max[i] = Double.MIN_VALUE;
            min[i] = Double.MAX_VALUE;
            // 获得每一列的均值、最大值、最小值
            for (int j = 0; j < row; j++) {
                if (max[i] < mat[j][i])
                    max[i] = mat[j][i];
                if (min[i] > mat[j][i])
                    min[i] = mat[j][i];
                sum += mat[j][i];
            }
            means[i] = sum / row;
            // 数据预处理
            for (int j = 0; j < row; j++) {
                if (max[i] == min[i])
                    mat[j][i] = 0;
                else
                    mat[j][i] = (mat[j][i] - means[i]) / (max[i] - min[i]);
            }
        }
    }
}
