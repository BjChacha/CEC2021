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

package etmo.metaheuristics.matmy2;

import etmo.metaheuristics.matmy2.libs.MaTAlgorithm;
import etmo.util.PseudoRandom;

import java.util.Arrays;

public class Utils {
    public static double distVector(double[] vector1, double[] vector2) {
        int dim = vector1.length;
        double sum = 0;
        for (int n = 0; n < dim; n++) {
            sum += (vector1[n] - vector2[n]) * (vector1[n] - vector2[n]);
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

    public static void randomPermutation(int[] perm, int size, int range){
        if (size > range){
            System.out.println("randomPermutation error: size should be less than range.");
            return;
        }
        int[] index = new int[range];
        randomPermutation(index, range);

        for (int i = 0; i < perm.length; i++){
            perm[i] = index[i];
        }
    }

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

    public static double[] softMaxOnlyPositive(double[] arr){
        double[] res = new double[arr.length];
        double sum = 0;
        double max = 0;
        for (int i = 0; i < arr.length; i++) {
            if (max < arr[i])
                max = arr[i];
        }
        for (int i = 0; i < res.length; i++){
            if (arr[i] > 0){
                res[i] = Math.exp(arr[i] / max);
                sum += res[i];
            }
            else{
                res[i] = 0;
            }
        }

        for (int i = 0; i < res.length; i++) {
            res[i] /= sum;
        }
        return res;
    }

    public static int roulette(double[] arr){
        double s = 0;
        double p = PseudoRandom.randDouble();
        int idx;
        // 轮盘赌算法
        for (idx = 0; idx < arr.length; idx++) {
            s += arr[idx];
            if (s >= p)
                break;
        }
        if (idx >= arr.length)
            idx = arr.length - 1;
        return idx;
    }

    public static double[] vectorMinus(double[] vec1, double[] vec2){
        if (vec1.length != vec2.length){
            System.out.println("Error: only can minus vectors with identical length.");
            return vec1;
        }

        double[] res = new double[vec1.length];
        for (int i = 0; i < res.length; i++)
            res[i] = vec1[i] - vec2[i];
        return res;
    }

    public static double calVectorAngle(double[] vec1, double[] vec2){
        if (vec1.length != vec2.length){
            System.out.println("Error: only can minus vectors with identical length.");
            return 1;
        }

        double v1 = 0;
        double v2 = 0;
        double dotP = 0;
        for (int i = 0; i < vec1.length; i++){
            v1 += vec1[i] * vec1[i];
            v2 += vec2[i] * vec2[i];
            dotP += vec1[i] * vec2[i];
        }
        v1 = Math.sqrt(v1);
        v2 = Math.sqrt(v2);

        double res = Math.abs(dotP / (v1 * v2));
        res = Math.acos(res);

        if (Double.isNaN(res) && v1 == 0)
            res = Double.MAX_VALUE;

        return res;
    }
}
