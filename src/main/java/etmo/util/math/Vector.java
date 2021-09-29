package etmo.util.math;

import java.util.Arrays;
import java.util.Random;

public class Vector {
    public static double[] vecSub(double[] vec1, double[] vec2) {
        assert vec1.length == vec2.length;
        double[] output = new double[vec1.length];
        Arrays.setAll(output, i -> vec1[i] - vec2[i]);
        return output;
    }

    public static void vecSub_(double[] vec1, double[] vec2) {
        assert vec1.length == vec2.length;
        Arrays.setAll(vec1, i -> vec1[i] - vec2[i]);
    }

    public static double[] vecAdd(double[] vec1, double[] vec2) {
        assert vec1.length == vec2.length;
        double[] output = new double[vec1.length];
        Arrays.setAll(output, i -> vec1[i] + vec2[i]);
        return output;
    }

    public static void vecAdd_(double[] vec1, double[] vec2) {
        assert vec1.length == vec2.length;
        // for (int i = 0; i < vec1.length; i ++) {
        //     vec1[i] = vec1[i] - vec2[i];
        // }
        Arrays.setAll(vec1, i -> vec1[i] + vec2[i]);
    }

    public static double vecDot(double[] vec1, double[] vec2) {
        assert vec1.length == vec2.length;
        double result = 0;
        for (int i = 0; i < vec1.length; i++) {
            result += vec1[i] * vec2[i];
        }
        return result;
    }

    public static double[] vecElemMul(double[] vec, double factor) {
        double[] output = new double[vec.length];
        Arrays.setAll(output, i -> vec[i] * factor);
        return output;
    }

    public static void vecElemMul_(double[] vec, double factor) {
        Arrays.setAll(vec, i -> vec[i] * factor);
    }

    public static double[] vecElemDiv(double[] vec1, double[] vec2) {
        double[] output = new double[vec1.length];
        Arrays.setAll(output, i -> vec1[i] / (vec2[i] + 1e-13));
        return output;
    }

    public static void vecClip_(double[] vec, double lb, double ub) {
        for (int i = 0; i < vec.length; i++) {
            vec[i] = Math.min(vec[i], ub);
            vec[i] = Math.max(vec[i], lb);
        }
    }

    public static double vecModule(double[] vector) {
        double module = 0;
        for (double e: vector) {
            module += Math.pow(e, 2);
        }
        module = Math.sqrt(module);
        return module;
    }

    public static double[] getRandomUnitVector(int d) {
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

    public static double[] vecUnify(double[] vector) {
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

    public static double[] vecNormailize(double[] vector) {
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

    public static double[] copy(double[] vector) {
        double[] vectorCopy = new double[vector.length];
        for(int i = 0; i < vector.length; i++) {
            vectorCopy[i] = vector[i];
        }
        return vectorCopy;
    }
}
