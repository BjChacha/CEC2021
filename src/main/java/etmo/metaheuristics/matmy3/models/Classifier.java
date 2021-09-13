package etmo.metaheuristics.matmy3.models;

import java.util.Arrays;

import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.listeners.ScoreIterationListener;
import org.nd4j.evaluation.classification.Evaluation;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.learning.config.AdaGrad;
import org.nd4j.linalg.lossfunctions.LossFunctions.LossFunction;


public class Classifier {
    int epoch;
    double lr;
    int inD;
    int hiddenD;

    MultiLayerNetwork net;

    public Classifier(int epoch, double lr, int inD, int hiddenD) {
        this.epoch = epoch;
        this.lr = lr;
        this.inD = inD;
        this.hiddenD = hiddenD;

        var config = new NeuralNetConfiguration.Builder()
            .weightInit(WeightInit.NORMAL)
            .updater(new AdaGrad(lr))
            .list()
            .layer(0, new DenseLayer.Builder().nIn(inD).nOut((int)(hiddenD*1.5)).activation(Activation.TANH).biasInit(0).build())
            .layer(1, new DenseLayer.Builder().nIn((int)(hiddenD*1.5)).nOut(hiddenD/2).activation(Activation.TANH).biasInit(0).build())
            // .layer(2, new DenseLayer.Builder().nIn(hiddenD*2).nOut(hiddenD).activation(Activation.TANH).biasInit(0).build())
            // .layer(3, new DenseLayer.Builder().nIn(hiddenD).nOut(hiddenD/2).activation(Activation.TANH).biasInit(0).build())
            .layer(2, new OutputLayer.Builder().nIn(hiddenD/2).nOut(1).activation(Activation.SIGMOID).lossFunction(LossFunction.XENT).build())
            .build();

        net = new MultiLayerNetwork(config);
        // net.setListeners(new ScoreIterationListener(10));
        net.init();
    }

    public void train(double[][] positive, double[][] negative) {
        int pLength = positive.length;
        int nLength = negative.length;
        INDArray pX = new NDArray(positive);
        INDArray nX = new NDArray(negative);
        DataSet data = DataSet.merge(Arrays.asList(
            new DataSet(pX, Nd4j.zeros(pLength, 1)),
            new DataSet(nX, Nd4j.ones(nLength, 1))));

        for (int i = 0; i < epoch; i++) {
            net.fit(data);
        }
    }

    public boolean judge(double[] x, double threshold) {
        return net.output(new NDArray(new double[][] {x})).getDouble(0) < threshold;
    }

    public double predict(double[] x) {
        return net.output(new NDArray(new double[][] {x})).getDouble(0); 
    }

    public void evaluate(double[][] testData) {
        Evaluation evaluator = new Evaluation(1);
        int length = testData.length / 2;
        // INDArray pLabel = Nd4j.zeros(length, 1);
        // INDArray nLabel = Nd4j.ones(length, 1);
        INDArray pLabel = Nd4j.ones(length, 1);
        INDArray nLabel = Nd4j.zeros(length, 1);
        INDArray labels = Nd4j.vstack(pLabel, nLabel);
        INDArray output = net.output(new NDArray(testData));
        evaluator.eval(labels, output);
        System.out.println(evaluator.stats());
    }
}
