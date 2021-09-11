package etmo.metaheuristics.matmy3;

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
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.learning.config.AdaGrad;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.learning.config.Sgd;
import org.nd4j.linalg.lossfunctions.LossFunctions.LossFunction;

public class EmsembleClassifier {
    int epoch;
    double lr;
    int inD;
    int classifierNum;

    MultiLayerNetwork[] nets;

    public EmsembleClassifier(int epoch, double lr, int inD, int classifierNum) {
        this.epoch = epoch;
        this.lr = lr;
        this.inD = inD;
        this.classifierNum = classifierNum;

        int splitedD = inD / classifierNum;
        nets = new MultiLayerNetwork[classifierNum];
        for (int i = 0; i < classifierNum; i ++) {
            var config = new NeuralNetConfiguration.Builder()
                .weightInit(WeightInit.XAVIER)
                .updater(new Adam(lr))
                // .l1(1e-3)
                .list()
                .layer(new BatchNormalization.Builder().nIn(splitedD).nOut(splitedD).build())
                .layer(new DenseLayer.Builder().nIn(splitedD).nOut((int)(splitedD * 1.5)).activation(Activation.TANH).biasInit(0).build())
                // .layer(new BatchNormalization.Builder().nIn((int)(splitedD * 1.5)).nOut((int)(splitedD * 1.5)).build())
                // .layer(new DenseLayer.Builder().nIn((int)(splitedD * 1.5)).nOut((int)(splitedD * 1.5)).activation(Activation.TANH).biasInit(0).build())
                .layer(new OutputLayer.Builder().nIn((int)(splitedD * 1.5)).nOut(1).activation(Activation.SIGMOID).lossFunction(LossFunction.XENT).build())
                .build();

            nets[i] = new MultiLayerNetwork(config);
            nets[i].setListeners(new ScoreIterationListener(10));
            nets[i].init();
        }
    }

    public void train(double[][] positive, double[][] negative) {
        int pLength = positive.length;
        int nLength = negative.length;
        double[][][] splitedPositive = splitData(positive);
        double[][][] splitedNegative = splitData(negative);
        for (int i = 0; i < classifierNum; i++) {
            INDArray pX = new NDArray(splitedPositive[i]);
            INDArray nX = new NDArray(splitedNegative[i]);
            DataSet data = DataSet.merge(Arrays.asList(
                new DataSet(pX, Nd4j.zeros(pLength, 1)),
                new DataSet(nX, Nd4j.ones(nLength, 1))));
            
            for (int e = 0; e < epoch; e++) {
                // for (DataSet batch: data.dataSetBatches(10)) {
                //     nets[i].fit(batch);
                // }
                nets[i].fit(data);
            }
        }
    }

    public boolean judge(double[] x, double classifiedThreshold, double voteThreshold) {
        double[][] splitedX = splitData(x);
        int positiveCount = 0;
        for (int i = 0; i < classifierNum; i++) {
            if (nets[i].output(new NDArray(new double[][] {splitedX[i]})).getDouble(0) < classifiedThreshold) {
                positiveCount ++;
            }
        }
        return ((double)positiveCount / classifierNum) >= voteThreshold;
    }

    private double[][][] splitData(double[][] data) {
        int n = data.length;
        int d = data[0].length;
        int splitedD = d / classifierNum;
        double[][][] splitedData = new double[classifierNum][n][splitedD];

        for (int i = 0; i < d; i++) {
            for (int j = 0; j < n; j ++) {
                splitedData[i/splitedD][j][i%splitedD] = data[j][i];
            }
        }

        return splitedData;
    }

    private double[][] splitData(double[] data) {
        int d = data.length;
        int splitedD = d / classifierNum;
        double[][] splitedData = new double[classifierNum][splitedD];
        for (int i = 0; i < d; i++) {
            splitedData[i/splitedD][i%splitedD] = data[i];
        }
        return splitedData;
    }

    public void evaluate(double[][] testData) {
        // int TP = 0, FP = 0, TN = 0, FN = 0;
        int length = testData.length / 2;
        INDArray pLabel = Nd4j.ones(length, 1);
        INDArray nLabel = Nd4j.zeros(length, 1);
        INDArray labels = Nd4j.vstack(pLabel, nLabel);
        double[][][] splitedTestData = splitData(testData);
        for (int i = 0; i < classifierNum; i++) {  
            Evaluation evaluator = new Evaluation(1);
            INDArray output = nets[i].output(new NDArray(splitedTestData[i]));
            evaluator.eval(labels, output);
            System.out.println("====== Classifier " + i + " ======");
            System.out.println(evaluator.stats());
        }
    }
}
