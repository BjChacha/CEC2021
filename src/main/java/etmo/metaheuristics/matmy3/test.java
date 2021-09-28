package etmo.metaheuristics.matmy3;

import java.io.IOException;
import java.net.URISyntaxException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.WindowConstants;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.deeplearning4j.datasets.iterator.impl.MnistDataSetIterator;
import org.deeplearning4j.optimize.listeners.callbacks.ModelSavingCallback;
import org.deeplearning4j.ui.api.UIServer;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries.XYSeriesRenderStyle;
import org.nd4j.autodiff.samediff.SDVariable;
import org.nd4j.autodiff.samediff.SameDiff;
import org.nd4j.autodiff.samediff.ops.SDValidation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;

import etmo.core.Solution;
import etmo.metaheuristics.matmy3.models.AbstractDistribution;
import etmo.metaheuristics.matmy3.models.CoralClassifier;
import etmo.metaheuristics.matmy3.models.GaussianDistribution;
import etmo.metaheuristics.matmy3.models.MultiVarGaussian;
import etmo.util.math.Matrix;
import smile.classification.GradientTreeBoost;
import smile.data.DataFrame;
import smile.data.formula.Formula;
import smile.manifold.TSNE;
import smile.stat.distribution.MultivariateGaussianDistribution;
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
        // GradientTreeBoost model = GradientTreeBoost.fit(Formula.lhs("class"), iris);
        // System.out.println(model.metrics().accuracy);
        // System.out.println(model.predict(iris.get(0)));
        // var X = DataFrame.of(x);
        // int[][] yy = new int[y.length][1];
        // for (int i = 0; i < yy.length; i++) {
        //     yy[i][0] = y[i];
        // }
        // var Y = DataFrame.of(yy, "class");
        // var df = X.merge(Y);

        // double[] vector = Matrix.randomUnitVector(10);
        // double sum = 0;
        // for (int i = 0; i < vector.length; i++) {
        //     sum += Math.pow(vector[i], 2);
        // }

        // double[][] data = Matrix.getRandomMat(100, 50);

        // // long time = System.currentTimeMillis();
        // // TSNE tSNE = new TSNE(data, 2);
        // // System.out.println(System.currentTimeMillis() - time + " ms.");
        // // Solution s = new Solution();
        // // solutionTest(s);

        // double[] mean = new double[50];
        // Arrays.parallelSetAll(mean, i -> { 
        //     double sum = 0;
        //     for (int j = 0; j < 50; j ++)
        //         sum += data[j][i];
        //     return sum / 50;
        // });

        // System.out.println(Arrays.toString(mean));

        // double[][] data = new double[][]{
        //     {1, 2, 3},
        //     {3, 1, 1}
        // };

        // double[][] sigma = Matrix.getMatSigma(data);
        // RealMatrix m = new Array2DRowRealMatrix(sigma);
        // SingularValueDecomposition svd = new SingularValueDecomposition(m);
        // DecompositionSolver solver = svd.getSolver();
        // double[][] U = svd.getU().getData();
        // double[][] S = svd.getS().getData();
        // double[][] VT = svd.getVT().getData();
        // for (int i = 0; i < S.length; i ++) {
        //     S[i][i] = Math.sqrt(S[i][i]);
        // }
        // double[][] result = Matrix.matMul(Matrix.matMul(U, S), VT);
        // double[][] newSigma = Matrix.matMul(result, result);

        // double[] mean = Matrix.getRandomMat(1, 50)[0];
        // double[][] sigma = Matrix.getMatSigma(data);

        // long startTime = System.currentTimeMillis();
        // commom math 3
        // MultivariateNormalDistribution model = new MultivariateNormalDistribution(mean, sigma);
        // for (int t = 0; t < 2500; t++) {
        //     double[] sample = model.sample();
        // }
        // System.out.println("Commom math: " + (System.currentTimeMillis() - startTime) + " ms.");


        // MultivariateNormalDistribution model[] = new MultivariateNormalDistribution[2500];
        // for (int t = 0; t < 2500; t++) {
        //     model[t] = new MultivariateNormalDistribution(mean, sigma);
        // }
        // System.out.println("Commom math: " + (System.currentTimeMillis() - startTime) + " ms.");

        // startTime = System.currentTimeMillis();
        // // smile
        // for (int t = 0; t < 2500; t++) {
        //     double[] sample = new MultivariateGaussianDistribution(
        //         mean, new smile.math.matrix.Matrix(sigma)).rand();
        // }
        // System.out.println("Smile: " + (System.currentTimeMillis() - startTime) + " ms.");

    //     double[][] data = new double[][]{
    //         {0.5280441909896626, 0.5317092013854541},
    //         {0.5288861365621463, 0.5330338888936633},
    //         {0.5250249651639478, 0.5341493724548296},
    //         {0.5285708205040301, 0.5408522739915376},
    //         {0.5423119526662601, 0.5305014406084265},
    //         {0.5286576909551626, 0.5316726115474976},
    //         {0.5293530126275507, 0.5316972801347171},
    //         {0.5288864128519828, 0.5316972801347171},
    //         {0.5218873904140025, 0.5315691230466589},
    //         {0.5299589219014511, 0.5323659780987745},
    //         {0.5284846323356659, 0.5270780637529834},
    //         {0.5284846323356659, 0.5272037765306486},
    //         {0.5284951044043213, 0.5316946445335982},
    //         {0.5423556145766488, 0.5317124712408449},
    //         {0.5250518230606084, 0.5271637386126421},
    //         {0.5299424282776206, 0.5312682177639744},
    //         {0.5268614636429236, 0.527251759231421},
    //         {0.5250518230606084, 0.5341223723093941},
    //         {0.5430057911237822, 0.5316726115474976},
    //         {0.5426501403609653, 0.5316940919850289},
    //         {0.5280441909896626, 0.5317121303430334},
    //         {0.5454040939299101, 0.531737814637191},
    //         {0.5427531780033176, 0.5304435517520941},
    //         {0.5283925657465497, 0.5316661323797821},
    //         {0.5272885448015285, 0.5316972801347171},
    //         {0.5277781787855917, 0.5316725548417589},
    //         {0.5298016614606409, 0.5273941629624276},
    //         {0.5286740316702699, 0.531662029361059},
    //         {0.5267920562855519, 0.5316777428617483},
    //         {0.5253439048341932, 0.5341493724548296},
    //         {0.5265112886107628, 0.5317091399208997},
    //         {0.5419447723101227, 0.5272037765306486},
    //         {0.5222322913345252, 0.5312857583078014},
    //         {0.5257810243266982, 0.5341223723093941},
    //         {0.5301156419688788, 0.5312627821918681},
    //         {0.5267648849016509, 0.5320012665251262},
    //         {0.5422218797054152, 0.5330180403006023},
    //         {0.5270717838950398, 0.5329586450710462},
    //         {0.5288861365621463, 0.533046614458267},
    //         {0.5428414359103461, 0.5316726115474976},
    //         {0.5288864128519828, 0.5317052009340133},
    //         {0.5420937051259337, 0.5408522739915376},
    //         {0.526792519858837, 0.5317394869504376},
    //         {0.5298027361913874, 0.5273941629624276},
    //         {0.5286752484771348, 0.5315695301872592},
    //         {0.5268493424199171, 0.5317091399208997},
    //         {0.5425863302869529, 0.5274754289175153},
    //         {0.5288952406962738, 0.5316972801347171},
    //         {0.5268016833728537, 0.5316777428617483},
    //         {0.5421127555112851, 0.5330180403006023},
    //         {0.5258309640953651, 0.5341223723093941},
    //         {0.5284846323356659, 0.5272037765306486},
    //         {0.5468235743031855, 0.5316972801347171},
    //         {0.5262745043614756, 0.5317518870226112},
    //         {0.5261486641202386, 0.5323277742331386},
    //         {0.5421127555112851, 0.5330180403006023},
    //         {0.5268470085760563, 0.532395082790657},
    //         {0.5454993825303095, 0.53235868669321},
    //         {0.5250812008131285, 0.5341493724548296},
    //         {0.5210212175051534, 0.5316726115474976},
    //         {0.5425539619726615, 0.5304435517520941},
    //         {0.5220427891520976, 0.5312681431020183},
    //         {0.5284846323356659, 0.5272037765306486},
    //         {0.5269359493433449, 0.5270678739103051},
    //         {0.5287254752272816, 0.5316713033111178},
    //         {0.5430057911237822, 0.5316726115474976},
    //         {0.5268608400988852, 0.5317091399208997},
    //         {0.5217548029215966, 0.5312624336066265},
    //         {0.5285117082232983, 0.5312527371081617},
    //         {0.5281461903334653, 0.5317091399208997},
    //         {0.5285763490134523, 0.5317077821219068},
    //         {0.5250057525791831, 0.5341223723093941},
    //         {0.5218500364632841, 0.5340897979477262},
    //         {0.5268452865488616, 0.5329586450710462},
    //         {0.5209997152653764, 0.533042766937593},
    //         {0.5268608400988852, 0.5317091399208997},
    //         {0.5268016833728537, 0.5316777428617483},
    //         {0.5286827074490844, 0.5317518870226112},
    //         {0.5287969524602701, 0.5321510827329735},
    //         {0.5284846323356659, 0.5272037765306486},
    //         {0.5219151166978585, 0.5312857583078014},
    //         {0.5285117082232983, 0.5312527371081617},
    //         {0.5423904704178418, 0.5272227506151314},
    //         {0.5268019629071169, 0.5270817520886434},
    //         {0.5268379020774062, 0.5311587434516909},
    //         {0.5454040939299101, 0.533042766937593}};

    //     double mean = .5;
    //     double std = .1;

    //     double[] means = Matrix.getMeanOfMat(data);
    //     double[] stds = Matrix.getStdOfMat(data);
    //     double[][] sigma = Matrix.getMatSigma(data);

    //     AbstractDistribution model = new GaussianDistribution(mean, std, 2);
    //     AbstractDistribution model1 = new MultiVarGaussian(means, sigma);
    //     MultivariateNormalDistribution model2 = new MultivariateNormalDistribution(means, sigma);
    //     List<Double> X = new ArrayList<>();
    //     List<Double> Y = new ArrayList<>();
    //     int count = 5000;
    //     for (int i = 0; i < count; i++) {
    //         double[] point = model2.sample();
    //         X.add(point[0]);
    //         Y.add(point[1]);
    //     }

    //     // for (int i = 0; i < data.length; i++) {
    //     //     X.add(data[i][0]);
    //     //     Y.add(data[i][1]);
    //     // }

    //     XYChart chart = new XYChartBuilder().build();
    //     chart.getStyler().setDefaultSeriesRenderStyle(XYSeriesRenderStyle.Scatter);
    //     // chart.getStyler().setXAxisMax(1.0);
    //     // chart.getStyler().setXAxisMin(0.0);
    //     // chart.getStyler().setYAxisMax(1.0);
    //     // chart.getStyler().setYAxisMin(0.0);
    //     chart.addSeries("gaussian_(" + mean + ", " + std + ")", X, Y);
    //     SwingWrapper<XYChart> sw = new SwingWrapper<>(chart);
    //     sw.displayChart().setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);


        // double[][] features = new double[][] {
        //     {1,1,1,1,1},
        //     {1,2,1,2,1},
        //     {2,2,2,2,2},
        //     {5,5,5,4,5},
        //     {4,5,5,4,3}
        // };

        // double[] labels = new double[] {0,0,0,1,1};


        // CoralClassifier model = new CoralClassifier();
        // model.train2(features, null, labels); 

        UIServer server = UIServer.getInstance();

        System.out.println("stop");

    }



    public static void solutionTest(Solution s) {
        s.setFlag(1);
    }
}