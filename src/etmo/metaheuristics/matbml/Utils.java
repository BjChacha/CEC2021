package etmo.metaheuristics.matbml;

import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.JMException;
import etmo.util.PseudoRandom;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Utils {
    static public List<List<Integer>> AGNES(SolutionSet[] populations, int clusterNum) throws JMException {
        List<List<Integer>> clusters = new ArrayList<>();
        // 1. each population is one cluster
        for (int k = 0; k < populations.length; k++){
            List<Integer> cluster = new ArrayList<>();
            cluster.add(k);
            clusters.add(cluster);
        }

        int[] idx = new int[2];
        while (clusters.size() > clusterNum) {
            // 2. find the closest clusters pair
            double minDistance = Double.MAX_VALUE;
            for (int i = 0; i < clusters.size() - 1; i++) {
                for (int j = i + 1; j < clusters.size(); j++) {
                    SolutionSet[] c1 = new SolutionSet[clusters.get(i).size()];
                    SolutionSet[] c2 = new SolutionSet[clusters.get(j).size()];
                    for (int ii = 0; ii < clusters.get(i).size(); ii++)
                        c1[ii] = populations[clusters.get(i).get(ii)];
                    for (int jj = 0; jj < clusters.get(j).size(); jj++)
                        c2[jj] = populations[clusters.get(j).get(jj)];

                    double dis = getMaxDistanceBetweenTasks(c1, c2);
                    if (dis < minDistance){
                        idx[0] = i;
                        idx[1] = j;
                        minDistance = dis;
                    }
                }
            }

            // 3. merge two clusters found above.
             for (int i = 0; i < clusters.get(idx[1]).size(); i++)
                clusters.get(idx[0]).add(clusters.get(idx[1]).get(i));
            clusters.remove(idx[1]);
        }

        return clusters;
    }

    static double getMaxDistanceBetweenTasks(SolutionSet[] PA, SolutionSet[] PB) throws JMException {
        double distance = 0;
        for (int k1 = 0; k1 < PA.length; k1++){
            for (int k2 = 0; k2 < PB.length; k2++){
                for (int i = 0; i < PA[k1].size(); i++){
                    for (int j = 0; j < PB[k2].size(); j++){
                        distance = Math.max(distance, getDistanceBetweenIndividuals(PA[k1].get(i), PB[k2].get(j)));
                    }
                }
            }
        }
        return distance;
    }

    static double getMinDistanceBetweenTasks(SolutionSet[] PA, SolutionSet[] PB) throws JMException {
        double distance = Double.MAX_VALUE;
        for (int k1 = 0; k1 < PA.length; k1++){
            for (int k2 = 0; k2 < PB.length; k2++){
                for (int i = 0; i < PA[k1].size(); i++){
                    for (int j = 0; j < PB[k2].size(); j++){
                        distance = Math.min(distance, getDistanceBetweenIndividuals(PA[k1].get(i), PB[k2].get(j)));
                    }
                }
            }
        }
        return distance;
    }

    static double getAvgDistanceBetweenTasks(SolutionSet[] PA, SolutionSet[] PB) throws JMException {
        double distance = 0;
        int CA = 0;
        int CB = 0;
        for (SolutionSet p: PA){ CA += p.size(); }
        for (SolutionSet p: PB){ CB += p.size(); }

        for (int k1 = 0; k1 < PA.length; k1++){
            for (int k2 = 0; k2 < PB.length; k2++){
                for (int i = 0; i < PA[k1].size(); i++){
                    for (int j = 0; j < PB[k2].size(); j++){
                        distance += getDistanceBetweenIndividuals(PA[k1].get(i), PB[k2].get(j));
                    }
                }
            }
        }
        distance /= (CA * CB);
        return distance;
    }

    static double getDistanceBetweenIndividuals(Solution pA, Solution pB) throws JMException {
        double distance = 0;
        int length = pA.getDecisionVariables().length;
        for (int i = 0; i < length; i++){
            distance += Math.pow((pA.getDecisionVariables(i) - pB.getDecisionVariables(i)), 2);
        }
        distance = Math.sqrt(distance);
        return distance;
    }

    static double getDistanceBetweenIndividuals(double[] pA, double[] pB) throws JMException {
        double distance = 0;
        int length = pA.length;
        for (int i = 0; i < length; i++){
            distance += Math.pow((pA[i] - pB[i]), 2);
        }
        distance = Math.sqrt(distance);
        return distance;
    }

    public static int[] permutation(int range, int size) {
        if (range < size){
            System.out.println("Error: range can not less than size.");
            return null;
        }
        int[] perm = new int[range];
        int[] index = new int[range];
        boolean[] flag = new boolean[range];

        for (int n = 0; n < range; n++) {
            index[n] = n;
            flag[n] = true;
        }

        int num = 0;
        while (num < range) {
            int start = etmo.util.PseudoRandom.randInt(0, range - 1);
            while (true) {
                if (flag[start]) {
                    perm[num] = index[start];
                    flag[start] = false;
                    num++;
                    break;
                }
                if (start == (range - 1)) {
                    start = 0;
                } else {
                    start++;
                }
            }
        } // while
        perm = Arrays.copyOfRange(perm, 0, size);
        return perm;
    }

    public static int roulette(double[] arr){
        arr = softmax(arr);

        double s = 0;
        double p = PseudoRandom.randDouble();
        int idx;

        for (idx = 0; idx < arr.length - 1; idx++) {
            s += arr[idx];
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
}
