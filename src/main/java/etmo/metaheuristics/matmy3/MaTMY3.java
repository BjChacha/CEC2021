package etmo.metaheuristics.matmy3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.WindowConstants;
import javax.swing.SwingUtilities;

import org.knowm.xchart.XYSeries;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.BitmapEncoder.BitmapFormat;

import etmo.core.MtoAlgorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.sorting.NDSortiong;
import etmo.qualityIndicator.QualityIndicator;

public class MaTMY3 extends MtoAlgorithm{
    private SolutionSet[] population;
    private SolutionSet[] offspring;

    private int populationSize;
    private int taskNum;

    private int evaluations;
    private int maxEvaluations;

    private String XType;

    private Operator crossover;
    private Operator mutation;
    
    // DEBUG: IGD
    String[] pf;
    List<QualityIndicator> indicators;
    double[] igd;

    // DEBUG: PLOT
    XYChart chart;
    List<XYChart> charts;
    SwingWrapper<XYChart> sw;
    int generation;
    List<Integer> generations;
    List<List<Double>> igdPlotValues;

    public MaTMY3(ProblemSet problemSet) {
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initState();
        initPlot();
        while (evaluations < maxEvaluations) {
            iterate();
            // long startTime = System.currentTimeMillis();
            updatePlot();
            // System.out.println("evaluations " + evaluations + "update plot time cost: " + (System.currentTimeMillis() - startTime) + " ms.");
        }
        endPlot();
        return population;
    }

    public void initState() throws JMException, ClassNotFoundException {
        evaluations = 0;
        taskNum = problemSet_.size();

        maxEvaluations = (Integer) this.getInputParameter("maxEvaluations");
        populationSize = (Integer) this.getInputParameter("populationSize");

        XType = (String) this.getInputParameter("XType");

        crossover = operators_.get("crossover");
        mutation = operators_.get("mutation");

        // initialize population
        population = new SolutionSet[taskNum];
        offspring = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
            population[k] = new SolutionSet(populationSize);
            offspring[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++) {
                Solution solution = new Solution(problemSet_);
                solution.setSkillFactor(k);
                problemSet_.get(k).evaluate(solution);
                evaluations++;
                population[k].add(solution);
            }
        }

        // DEBUG: IGD
        pf = new String[taskNum];
        indicators = new ArrayList<>(taskNum);
        for (int k = 0; k < taskNum; k++) {
            pf[k] = "resources/PF/StaticPF/" + problemSet_.get(k).getHType() + "_" + problemSet_.get(k).getNumberOfObjectives() + "D.pf";
            indicators.add(new QualityIndicator(problemSet_.get(k), pf[k]));
        }

        // DEBUG: PLOT
        generation = 0;
        generations = new ArrayList<>();
        igdPlotValues =  new ArrayList<>();
        for (int k = 0; k < taskNum; k++) {
            igdPlotValues.add(new ArrayList<>());
        }
    }

    public void iterate() throws JMException {
        offspringGeneration(XType);
        environmentSelection();
    }

    public void offspringGeneration(String type) throws JMException {
        for (int k = 0; k < taskNum; k++) {
            offspring[k].clear();
            if (type.equalsIgnoreCase("SBX")) {
                for (int i = 0; i < populationSize; i++) {
                    int j = i;
                    while (j == i) 
                        j = PseudoRandom.randInt(0, populationSize - 1);
                    Solution[] parents = new Solution[2];
                    parents[0] = population[k].get(i);
                    parents[1] =population[k].get(j);
                    
                    Solution child = ((Solution[]) crossover.execute(parents))[PseudoRandom.randInt(0,1)];
                    mutation.execute(child);

                    child.setSkillFactor(k);
                    problemSet_.get(k).evaluate(child);
                    evaluations++;
                    offspring[k].add(child);
                }
            }
            else if (type.equalsIgnoreCase("DE")) {
                for (int i = 0; i < populationSize; i++) {
                    int j1 = i, j2 = i;
                    while (j1 == i && j1 == j2){
                        j1 = PseudoRandom.randInt(0, populationSize - 1);
                        j2 = PseudoRandom.randInt(0, populationSize - 1);
                    }
                    Solution[] parents = new Solution[3];
                    parents[0] = population[k].get(j1);
                    parents[1] = population[k].get(j2);
                    parents[2] = population[k].get(i);
                    
                    Solution child = (Solution) crossover.execute(new Object[] {
                        population[k].get(i), parents });
                    mutation.execute(child);

                    child.setSkillFactor(k);
                    problemSet_.get(k).evaluate(child);
                    evaluations++;
                    offspring[k].add(child);
                }
            }
            else {
                System.out.println("Error: unsupported reproduce type: " + type);
                System.exit(1);
            }
        }
    }

    public void environmentSelection() {
        for (int k = 0; k <  taskNum; k++) {
            SolutionSet union = population[k].union(offspring[k]);
            NDSortiong.sort(union, problemSet_, k);
            for (int i = 0; i < populationSize; i++) {
                population[k].replace(i, union.get(i));
            }
        }
    }

    // DEBUG: IGD
    private void calIGD() {
        igd = new double[taskNum];
        for (int k = 0; k < taskNum; k++) {
            igd[k] = indicators.get(k).getIGD(population[k], k);
            igdPlotValues.get(k).add(igd[k]);
        }
        generations.add(generation++);
        // System.out.println("Evaluations " + evaluations + ": " + Arrays.toString(igd));
    }

    public double[] getPlotX() {
        return generations.stream().mapToDouble(d->d).toArray();
    }

    public double[] getPlotX(int maxLength) {
        double[] x = getPlotX();
        return x.length > maxLength ? Arrays.copyOfRange(x, x.length - maxLength, x.length) : x;
    }

    public double[][] getPlotY() {
        double[][] y = new double[taskNum][];
        for (int k = 0; k < taskNum; k++) {
            y[k] = igdPlotValues.get(k).stream().mapToDouble(d->d).toArray();
        }
        return y;
    }

    public double[][] getPlotY(int maxLength) {
        double[][] y = new double[taskNum][];
        for (int k = 0; k < taskNum; k++) {
            double[] tmp = igdPlotValues.get(k).stream().mapToDouble(d->d).toArray();
            y[k] = tmp.length > maxLength ? Arrays.copyOfRange(tmp, tmp.length - maxLength, tmp.length) : tmp;
        }
        return y;
    }
    
    public void initPlot() {
        calIGD();
        double[] x = getPlotX();
        double[][] y = getPlotY();
        chart = new XYChartBuilder()
            .title("Generation: " + generation)
            .xAxisTitle("Generation")
            .yAxisTitle("IGD")
            .build();
        for (int k = 0; k < taskNum; k++) {
            chart.addSeries("Problem " + k, x, y[k]);
        }
    
        chart.getStyler().setYAxisLogarithmic(true);
        // sw = new SwingWrapper<XYChart>(chart);
        // sw.displayChart().setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
    }

    public void updatePlot() {
        calIGD();
        double[] x = getPlotX();
        double[][] y = getPlotY();

        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                chart.setTitle("Generation: " + generation);
                for (int k = 0; k < taskNum; k++){
                    chart.updateXYSeries("Problem " + k, x, y[k], null);
                }
                // sw.repaintChart();
            }
        });
    }

    public void initPlot2() {
        calIGD();
        double[] x = getPlotX();
        double[][] y = getPlotY();

        charts = new ArrayList<XYChart>();
        for (int k = 0; k < taskNum; k++) {
            XYChart c = new XYChartBuilder()
                // .title("Problem " + k + " Generation: " + generation)
                .title(Integer.toString(k))
                .xAxisTitle("Generation")
                .yAxisTitle("IGD")
                .width(300)
                .height(200)
                .build();
            c.getStyler().setYAxisLogarithmic(true);
            XYSeries s = c.addSeries("Problem " + k, x, y[k]);
            charts.add(c);
        }
        sw = new SwingWrapper<XYChart>(charts);
        sw.displayChartMatrix();
    }

    public void updatePlot2() {
        calIGD();
        double[] x = getPlotX();
        double[][] y = getPlotY();

        for (int k = 0; k < taskNum; k++) {
            final int taskID = k;
            SwingUtilities.invokeLater(new Runnable() {
                    @Override
                    public void run() {
                        charts.get(taskID).updateXYSeries("Problem " + taskID, x, y[taskID], null);
                    }
                });
            sw.repaintChart();
        }

        // SwingUtilities.invokeLater(new Runnable() {
        //     @Override
        //     public void run() {
        //         for (int k = 0; k < taskNum; k++){
        //             // charts.get(k).setTitle("Generation: " + generation);
        //             charts.get(k).updateXYSeries("Problem " + k, x, y[k], null);
        //         }
        //         sw.repaintChart();
        //     }
        // });
    }

    public void endPlot() {
        try {
            BitmapEncoder.saveBitmap(chart, "./figs/" + problemSet_.get(0).getName(), BitmapFormat.PNG);
            // VectorGraphicsEncoder.saveVectorGraphic(chart, "./figs/" + problemSet_.get(0).getName(), VectorGraphicsFormat.PDF);
        } catch (IOException e) {
            e.printStackTrace();
        }
        
    }
}
