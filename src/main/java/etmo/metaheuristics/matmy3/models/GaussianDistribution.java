package etmo.metaheuristics.matmy3.models;

import java.util.Arrays;
import java.util.Random;

public class GaussianDistribution extends AbstractDistribution {
    int d;
    double[] means;
    double[] stds;
    
    public GaussianDistribution(double mean, double std) {
        this.d = 1;
        this.means = new double[] {mean};
        this.stds = new double[] {std};
    }

    public GaussianDistribution(double mean, double std, int d) {
        this.d = d;
        this.means = new double[d];
        this.stds = new double[d];
        Arrays.fill(this.means, mean);
        Arrays.fill(this.stds, std);
    }

    public GaussianDistribution(double[] means, double[] stds) {
        assert means.length == stds.length;
        this.d = means.length;
        this.means = means.clone();
        this.stds = stds.clone();
    }

    public double[] sample() {
        double[] sample = new double[d];
        Random randomGenerator = new Random();
        Arrays.parallelSetAll(sample, i -> randomGenerator.nextGaussian() * stds[i] + means[i]);
        return sample;
    }

    public double[][] sample(int count) {
        double[][] samples = new double[count][d];
        Arrays.parallelSetAll(samples, i -> sample());
        return samples;
    }
}
