package etmo.util.plot;

import java.lang.reflect.InvocationTargetException;

import org.knowm.xchart.QuickChart;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;

import etmo.util.errorcheck.Check;
import smile.plot.swing.ScatterPlot;
import smile.plot.swing.Canvas;
import smile.math.kernel.BinarySparseGaussianKernel;

public class Plotter {
    private double[][] plotData;
    private String plotTitle;

    public Plotter(double[][] matrix) { this(matrix, "Front"); }

    public Plotter(double[][] matrix, String title) {
        Check.notNull(matrix);
        Check.that(matrix.length >= 1, "The data matrix is empty.");

        this.plotData = matrix;
        this.plotTitle = title;
    }

    public void plot() throws InvocationTargetException, InterruptedException {
        var canvas = ScatterPlot.of(plotData, '*').canvas();
        canvas.setTitle(plotTitle);
        canvas.window();
    }

    public void animPlotTest(int times) throws InvocationTargetException, InterruptedException{
        double[][] data = plotData.clone();
        Canvas canvas = null;
        long waitSecond = (long)(0.5 * 1000);  
        for (int t = 0; t < times; t++){
            if (canvas == null){
                canvas = ScatterPlot.of(data, '*').canvas();
                canvas.window();
            }else{
                canvas.add(ScatterPlot.of(data, '*'));
            }

            Thread.sleep(waitSecond);

            data[0][1] += 1;

            canvas.clear();
        }
    }

    public void animPlotTest2(int times) throws InvocationTargetException, InterruptedException{
        double[][] data = plotData.clone();
        
        final XYChart chart = QuickChart.getChart("Demo", "x", "y", "P", data[0], data[1]);
        final SwingWrapper<XYChart> sw = new SwingWrapper<XYChart>(chart);
        sw.displayChart();
        for (int t = 0; t < times; t++){
            Thread.sleep(200);
            
            data[0][0] = t % 3;
            javax.swing.SwingUtilities.invokeLater(new Runnable(){
                @Override
                public void run(){
                    chart.updateXYSeries("P", data[0], data[1], null);
                    sw.repaintChart();
                }
            });
        }
    }
    double[] x = {1,2,3};
    double[] y = {1,2,3};
}
