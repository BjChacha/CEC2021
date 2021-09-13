package etmo.metaheuristics.matmy3.models;

import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.listeners.ScoreIterationListener;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.learning.config.AdaGrad;
import org.nd4j.linalg.lossfunctions.LossFunctions.LossFunction;

public class Regression {
    int epoch;
    double lr;
    int inD;
    int outD;
    int hiddenD;

    MultiLayerNetwork net;

    public Regression(int epoch, double lr, int inD, int outD, int hiddenD) {
        this.epoch = epoch;
        this.lr = lr;
        this.inD = inD;
        this.outD = outD;
        this.hiddenD = hiddenD;

        var config = new NeuralNetConfiguration.Builder()
            .weightInit(WeightInit.NORMAL)
            .updater(new AdaGrad(lr))
            .list()
            .layer(0, new DenseLayer.Builder().nIn(inD).nOut(hiddenD).activation(Activation.TANH).biasInit(0).build())
            .layer(1, new OutputLayer.Builder().nIn(hiddenD).nOut(outD).activation(Activation.HARDSIGMOID).lossFunction(LossFunction.MSE).build())
            .build();
        net = new MultiLayerNetwork(config);
        net.setListeners(new ScoreIterationListener(10));
        net.init();
    }

    public void train(double[][] x, double[][] y) {
        INDArray input = new NDArray(x);
        INDArray output = new NDArray(y);

        net.clear();
        for (int i = 0; i < epoch; i++) {
            net.fit(input, output);
        }
    }    

    public double[][] predict(double[][] x) {
        INDArray input = new NDArray(x);
        INDArray output = net.output(input);
        return output.toDoubleMatrix();
    }

    public void reset() {
        var config = new NeuralNetConfiguration.Builder()
        .weightInit(WeightInit.ZERO)
        .updater(new AdaGrad(lr))
        .list()
        .layer(0, new DenseLayer.Builder().nIn(inD).nOut(hiddenD).activation(Activation.TANH).build())
        .layer(1, new OutputLayer.Builder().nIn(hiddenD).nOut(outD).activation(Activation.HARDSIGMOID).build())
        .build();
        net = new MultiLayerNetwork(config);
        net.setListeners(new ScoreIterationListener(10));
        net.init();
    }

    public void setLr(double lr) {
        net.setLearningRate(lr);
    }
}
