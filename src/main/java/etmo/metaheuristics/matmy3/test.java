package etmo.metaheuristics.matmy3;

import java.io.IOException;
import java.net.URISyntaxException;
import java.text.ParseException;

import smile.data.DataFrame;
import smile.data.formula.Formula;
import smile.io.Read;
import weka.classifiers.trees.RandomForest;

public class test {
    public static void main(String[] args) throws IOException, ParseException, URISyntaxException {
        // INDArray x = Nd4j.rand(50, 50);
        // INDArray y = Nd4j.rand(50, 50);
        
        // Regression rg = new Regression(20, 5e-1, 50, 50, 50);
        // rg.train(x.toDoubleMatrix(), y.toDoubleMatrix());
        
        // var iris = Read.arff("data/weka/iris.arff");
        // var x = iris.drop("class").toArray();
        // var y = iris.column("class").toIntArray();
        // RandomForest model = RandomForest.fit(Formula.lhs("class"), iris);
        // System.out.println(model.metrics().accuracy);
        // System.out.println(model.predict(iris.get(0)));
        // var X = DataFrame.of(x);
        // int[][] yy = new int[y.length][1];
        // for (int i = 0; i < yy.length; i++) {
        //     yy[i][0] = y[i];
        // }
        // var Y = DataFrame.of(yy, "class");
        // var df = X.merge(Y);


    }
}