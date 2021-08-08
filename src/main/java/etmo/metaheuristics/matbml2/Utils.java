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

package etmo.metaheuristics.matbml2;

import etmo.util.PseudoRandom;

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
}
