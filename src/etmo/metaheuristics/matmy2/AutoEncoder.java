package etmo.metaheuristics.matmy2;

import org.deeplearning4j.nn.api.OptimizationAlgorithm;
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

public class AutoEncoder {
    private static final int hiddenNum = 20;

    int inN_;
    int outN_;
    int hN_;

    int dN_;

    NDArray input_;
    NDArray output_;

    MultiLayerNetwork net;

    public AutoEncoder(double[][] features, double[][] labels){
        dN_ = features.length;
        inN_ = features[0].length;
        outN_ = labels[0].length;
        input_ = new NDArray(features);
        output_ = new NDArray(labels);
        hN_ = hiddenNum;
    }

    public AutoEncoder(NDArray features, NDArray labels){
        dN_ = features.rows();
        inN_ = features.columns();
        outN_ = labels.columns();
        input_ = features;
        output_ = features;
        hN_ = hiddenNum;
    }

    public void train(){
        // 网络初始化
        var conf = new NeuralNetConfiguration.Builder()
                .seed(12345)
                .weightInit(WeightInit.XAVIER)
                .updater(new AdaGrad(0.05))
                .activation(Activation.RELU)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT)
                .l2(0.0001)
                .list()
                .layer(0, new DenseLayer.Builder().nIn(inN_).nOut(hN_).build())
                .layer(1, new OutputLayer.Builder().nIn(hN_).nOut(outN_).build())
                .validateOutputLayerConfig(false)
                .build();

        net = new MultiLayerNetwork(conf);
        net.setListeners(new ScoreIterationListener(1));

        for (int i = 0; i < 5; i++){
            net.fit(input_, output_);
        }
    }

    public double[][] predict(double[][] input){
        NDArray in = new NDArray(input);
        INDArray out = net.output(in);
        return out.toDoubleMatrix();
    }

    public double[][] predict(NDArray input){
        INDArray out = net.output(input);
        return out.toDoubleMatrix();
    }
}
