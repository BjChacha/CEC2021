package etmo.util.math;

import java.util.Random;

public class Probability {
    public static double sampleByNorm(double mean, double std) {
        return new Random().nextGaussian() * std + mean;
    }

    public static double[] sampleByNorm(double[] mean, double[] std) {
        double[] sample = new double[mean.length];
        for (int i = 0; i < sample.length; i++) {
            sample[i] = sampleByNorm(mean[i], std[i]);
            sample[i] = Math.max(0, sample[i]);
            sample[i] = Math.min(1, sample[i]);
        }
        return sample;
    }

}
