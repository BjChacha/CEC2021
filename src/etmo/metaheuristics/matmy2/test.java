package etmo.metaheuristics.matmy2;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.api.rng.Random;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.IOException;
import java.util.Arrays;

public class test {
    public static void main(String[] args) throws IOException {
//        double[] vector1 = new double[]{1, 0};
//        double[] vector2 = new double[]{1, 1};
//        double[] vector3 = new double[]{0, 1};
//        double[] vector4 = new double[]{-1, 0};
//
//        // 0
//        System.out.println(Utils.calVectorAngle(vector1, vector1));
//        // 0.78
//        System.out.println(Utils.calVectorAngle(vector1, vector2));
//        // 1.57
//        System.out.println(Utils.calVectorAngle(vector1, vector3));
//        // 0
//        System.out.println(Utils.calVectorAngle(vector1, vector4));
//        for (int i = 0; i < 10; i++)
//        {
//            System.out.println(PseudoRandom.randInt(0, 1));
//        }

//        double[][] mat1 = {{-1,1,0},{-4,3,0},{1,0,2}};
////        double[][] mat2 = {{1,3,2},{2,4,2},{3,7,2}};
//        double[][] mat3 = {
//                {1,3,5,2,1},
//                {2,4,7,2,2},
//                {4,9,5,5,8},
//                {9,2,7,5,9},
//                {1,2,3,4,1},
//        };
//        double[][] mat4 = {
//                {9,2,7,5,9},
//                {2,4,7,2,2},
//                {1,3,5,2,1},
//                {1,2,3,4,1},
//                {4,9,5,5,8},
//        };
////        double[][] res1 = Utils.MappingViaPCA(mat1, mat2);
//
////        INDArray in = Nd4j.rand(100, 50);
////        INDArray out = Nd4j.rand(100, 50);
////
//////        AutoEncoder ae = new AutoEncoder(in.toDoubleMatrix(), out.toDoubleMatrix());
//////        ae.train();
//////        ae.predict(in.toDoubleMatrix());
////        double[][] res = Utils.MappingViaPCA1(in.toDoubleMatrix(), out.toDoubleMatrix());
////        double[][] res = Utils.Sample(mat1, 1);
////        LTR.test(mat3);
//        Utils.MappingViaGAN(mat3, mat4, 10);
//        System.out.println("a");

        int times = 5;
        double[] tmp = new double[]{1, 3, 5, 7, 9};
        Arrays.sort(tmp);
        double best, worst, mean, median, std = 0;
        best = tmp[0];
        worst = tmp[times-1];
        mean = Arrays.stream(tmp).sum() / times;
        median = tmp[times/2];
        for (double e: tmp){
            std += Math.pow(e - mean, 2);
        }
        std = Math.sqrt(std / times);
        System.out.println("\tBest: " + best);
        System.out.println("\tWorst: " + worst);
        System.out.println("\tMean: " + mean);
        System.out.println("\tMedian: " + median);
        System.out.println("\tStd: " + std);
    }
}
