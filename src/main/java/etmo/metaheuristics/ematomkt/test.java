package etmo.metaheuristics.ematomkt;

import etmo.util.math.Matrix;

public class test {
    public static void main(String[] args) {
        // double[][] a = {{0.7831,0.7365,0.6442,0.5060}, {0.3743,0.1399,0.3644,0.3721}, {0.5057,0.8613,0.4254,0.4585}};
        // double[][] b = {{0.9000,0.9007,0.2387,0.8727}, {0.9372,0.0356,0.3433,0.1190}, {0.0617,0.1486,0.6308,0.7271}};
        // Utils.MMD(a,b, 1.0);
        // for (int i = 0; i < 10; i++)
        //     System.out.println(Arrays.toString(Utils.randomPermutation(10, 3)));

        double[][] union = {
            {0.1,0.2,0.3,0.4,0.5},
            {0.1,0.1,0.1,0.1,0.1},
            {0.1,0.2,0.2,0.3,0.2},
            {0.2,0.2,0.3,0.3,0.3},
            {0.4,0.4,0.4,0.4,0.4},
            {0.4,0.5,0.5,0.6,0.6},
            {0.6,0.5,0.4,0.3,0.4},
            {0.5,0.1,0.2,0.3,0.5},
            {0.9,0.9,0.9,0.8,0.8},
            {0.7,0.7,0.7,0.7,0.7}};
        // double[][][] res = Utils.kmeans(union, 3);
        double[] mean = Matrix.getMeanOfMat(union);
        double[] std = Matrix.getStdOfMat(union);
        System.out.println("pause");
    }
}
