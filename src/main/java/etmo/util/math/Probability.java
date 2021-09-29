package etmo.util.math;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

public class Probability {
    public static double sampleByNorm(double mean, double std) {
        return new Random().nextGaussian() * std + mean;
    }

    public static double[] sampleByNorm(double[] mean, double[] std) {
        double[] sample = new double[mean.length];
        Arrays.setAll(sample, i -> sampleByNorm(mean[i], std[i]));
        return sample;
    }

    public static double[] sampleByNorm(double[] mean, double[][] sigma) {
        double[][] simgaClone = sigma.clone();      
        for (int i = 0; i < simgaClone.length; i++) {
            simgaClone[i][i] *= 1.001;
        }

        double[] sample = new double[mean.length];
        Random randomer = new Random();
        for (int i = 0; i < sample.length; i++) {
            sample[i] = randomer.nextGaussian();
        }

        RealMatrix mX = new Array2DRowRealMatrix(sample);
        RealMatrix mSigma = new Array2DRowRealMatrix(simgaClone);
        SingularValueDecomposition svd = new SingularValueDecomposition(mSigma);
        RealMatrix S = svd.getS();
        for (int i = 0; i < S.getColumnDimension(); i++) {
            S.setEntry(i, i, Math.sqrt(S.getEntry(i, i)));
        }
        RealMatrix mSigmaSqrt = svd.getVT().multiply(svd.getU().multiply(S));
        RealMatrix mMean =  new Array2DRowRealMatrix(mean);
        RealMatrix result = mSigmaSqrt.multiply(mX).add(mMean);
        return result.transpose().getData()[0];
    }
}
