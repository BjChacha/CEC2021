package etmo.test;

import java.lang.reflect.InvocationTargetException;

import etmo.util.plot.Plotter;
import smile.math.kernel.PearsonKernel;
import smile.math.kernel.MercerKernel;
import smile.regression.GaussianProcessRegression;

public class Test_Plot {
    public static void main(String[] args) throws InvocationTargetException, InterruptedException{
        double[][] data = {
            {1, 2, 2},
            {2, 3, 4},
        };
        Plotter plotter = new Plotter(data);
        // plotter.plot();
        plotter.animPlotTest2(100);

        double[][] x = {{1,2,3}, {2,3,4}};
        double[] y = {1,2,3};
        MercerKernel<double[]> kernel = new PearsonKernel(0.1, 0.1);
        GaussianProcessRegression.<double[]>fit(x, y, kernel, 0.1);
    }
}
