package etmo.metaheuristics.matmy3.models;

import java.util.Arrays;

import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.layers.BatchNormalization;
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
import org.nd4j.linalg.dataset.SplitTestAndTrain;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.learning.config.AdaGrad;
import org.nd4j.linalg.lossfunctions.LossFunctions.LossFunction;


public class Classifier {
    int epoch;
    double lr;
    int inD;
    int hiddenD;
    double f1;

    MultiLayerNetwork net;

    public Classifier(int epoch, double lr, int inD, int hiddenD) {
        this.epoch = epoch;
        this.lr = lr;
        this.inD = inD;
        this.hiddenD = hiddenD;

        var config = new NeuralNetConfiguration.Builder()
            .weightInit(WeightInit.NORMAL)
            .updater(new AdaGrad(lr))
            .l2(2e-5)
            .list()
            // .layer(new BatchNormalization.Builder().nIn(inD).nOut(inD).build())
            .layer(new DenseLayer.Builder().nIn(inD).nOut((int)(hiddenD / 2)).activation(Activation.TANH).biasInit(0).build())
            .layer(new DenseLayer.Builder().nIn(hiddenD / 2).nOut(hiddenD / 4).activation(Activation.TANH).biasInit(0).build())
            .layer(new OutputLayer.Builder().nIn(hiddenD / 4).nOut(1).activation(Activation.SIGMOID).biasInit(0).lossFunction(LossFunction.XENT).build())
            .build();

        net = new MultiLayerNetwork(config);
        net.setListeners(new ScoreIterationListener(epoch / 5));
        net.init();
    }

    public void train(double[][] positive, double[][] negative) {
        train(positive, negative, epoch);
    }

    public void train(double[][] positive, double[][] negative, int newEpoch) {
        int pLength = positive.length;
        int nLength = negative.length;
        INDArray pX = new NDArray(positive);
        INDArray nX = new NDArray(negative);
        DataSet data = DataSet.merge(Arrays.asList(
            new DataSet(pX, Nd4j.zeros(pLength, 1)),
            new DataSet(nX, Nd4j.ones(nLength, 1))));

        SplitTestAndTrain split = data.splitTestAndTrain(0.9);
        DataSet train = split.getTrain();
        DataSet test = split.getTest();
        for (int i = 0; i < newEpoch; i++) {
            net.fit(train);
        }
        Evaluation evaluator = new Evaluation(1);
        evaluator.eval(test.getLabels(), net.output(test.getFeatures()));
        this.f1 = evaluator.f1();
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
        INDArray pLabel = Nd4j.zeros(length, 1);
        INDArray nLabel = Nd4j.ones(length, 1);
        // INDArray pLabel = Nd4j.ones(length, 1);
        // INDArray nLabel = Nd4j.zeros(length, 1);
        INDArray labels = Nd4j.vstack(pLabel, nLabel);
        INDArray output = net.output(new NDArray(testData));
        evaluator.eval(labels, output);
        System.out.println(evaluator.stats());
    }

    public double getF1() {
        return this.f1;
    }

}
