package etmo.metaheuristics.matmy3;

import java.util.List;

import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.learning.config.Sgd;



public class Mlp {
    int epoch;
    float lr;
    List<Integer> hiddenUnits;

    MultiLayerNetwork net;

    public Mlp(int epoch, float lr, List<Integer> hiddenUnits) {
        this.epoch = epoch;
        this.lr = lr;
        this.hiddenUnits = hiddenUnits;

        // var config = new NeuralNetConfiguration.Builder()
        //     .weightInit(WeightInit.ZERO)
        //     .updater(new Sgd(lr))
        //     .

    }

    public void train(double[][] x, double[][] y) {

    }    
}
