package etmo.metaheuristics.matmy3;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

public class test {
    public static void main(String[] args) {
        INDArray x = Nd4j.rand(50, 50);
        INDArray y = Nd4j.rand(50, 50);
        
        Regression rg = new Regression(20, 5e-1, 50, 50, 50);
        rg.train(x.toDoubleMatrix(), y.toDoubleMatrix());
        
        System.out.println("stop");
    }
}