package etmo.metaheuristics.matmy3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.WindowConstants;

import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.BitmapEncoder.BitmapFormat;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.XYSeries.XYSeriesRenderStyle;
import org.knowm.xchart.style.colors.XChartSeriesColors;
import org.knowm.xchart.style.lines.XChartSeriesLines;
import org.knowm.xchart.style.markers.XChartSeriesMarkers;

import etmo.core.MtoAlgorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet; 
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.sorting.NDSortiong;

public class MaTMY3_Transfer extends MtoAlgorithm {
    private SolutionSet[] population;
    private SolutionSet[] offspring;

    private int generation;
    private int evaluations;
    private int maxEvaluations;

    private int populationSize;
    private int taskNum;
    private int varNum;
    private int[] objStart;
    private int[] objEnd;

    private int[] transferTotalCount;

    private String XType;
    private Operator crossover;
    private Operator mutation;

    private int[] stuckTimes;
    private double[] bestDistances;

    private double transferProbability;

    boolean isMutate;

    // DEBUG: IGD
    String[] pf;
    List<QualityIndicator> indicators;
    double[] igd;

    // DEBUG: PLOT
    boolean isPlot;
    int plotTaskID = 0;
    XYChart chartIGD;
    XYChart chartPF;
    XYChart chartVar;
    List<XYChart> charts;
    SwingWrapper<XYChart> sw;
    SwingWrapper<XYChart> sw2;
    SwingWrapper<XYChart> sw3;
    List<Integer> generations;
    List<List<Double>> igdPlotValues;

    public MaTMY3_Transfer(ProblemSet problemSet) {
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initState();
        if (isPlot)
            initPlot();
        while (evaluations < maxEvaluations) {
            iterate();
            // long startTime = System.currentTimeMillis();
            if (isPlot)
                updatePlot();
            // System.out.println("evaluations " + evaluations + "update plot time cost: " +
            // (System.currentTimeMillis() - startTime) + " ms.");

            // System.out.println(evaluations + ": " + Arrays.toString(stuckTimes));
        }
        // if (isPlot)
        // endPlot();

        return population;
    }

    public void initState() throws JMException, ClassNotFoundException {
        generation = 0;
        evaluations = 0;
        taskNum = problemSet_.size();
        varNum = problemSet_.getMaxDimension();

        maxEvaluations = (Integer) this.getInputParameter("maxEvaluations");
        populationSize = (Integer) this.getInputParameter("populationSize");

        // // DEBUG partly run
        // maxEvaluations /= taskNum;
        // problemSet_ = problemSet_.getTask(0);
        // taskNum = 1;

        XType = (String) this.getInputParameter("XType");
        isPlot = (Boolean) this.getInputParameter("isPlot");
        isMutate = (Boolean) this.getInputParameter("isMutate");

        transferProbability = (Double) this.getInputParameter("transferProbability");

        crossover = operators_.get("crossover");
        mutation = operators_.get("mutation");

        objStart = new int[taskNum];
        objEnd = new int[taskNum];
        bestDistances = new double[taskNum];
        stuckTimes = new int[taskNum];
        Arrays.fill(bestDistances, Double.MAX_VALUE);
        Arrays.fill(stuckTimes, 0);

        population = new SolutionSet[taskNum];
        offspring = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
            objStart[k] = problemSet_.get(k).getStartObjPos();
            objEnd[k] = problemSet_.get(k).getEndObjPos();

            population[k] = new SolutionSet(populationSize);
            offspring[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++) {
                Solution solution = new Solution(problemSet_);
                solution.setSkillFactor(k);
                // double[] features = new double[]{PseudoRandom.randDouble(), 0.49802832, 0.53675057, 0.45468875, 0.52329801, 0.488946  ,0.51170569, 0.45461273, 0.50697831, 0.47060217, 0.49723718,0.53323252, 0.50307031, 0.46455359, 0.5280332 , 0.50182858,0.46595516, 0.52930945, 0.52279787, 0.4700964 , 0.46489979,0.50416938, 0.51328641, 0.49796944, 0.51451951, 0.46973632,0.54227379, 0.54462272, 0.52398501, 0.53820799, 0.4534272 ,0.46881809, 0.50857003, 0.54416919, 0.50364375, 0.45100508,0.48236421, 0.450201  , 0.50381148, 0.54630262, 0.49303332,0.47703986, 0.54664697, 0.51835294, 0.50294303, 0.47206974,0.4557836 , 0.48946816, 0.53982868, 0.51843245};
                // solution.setDecisionVariables(features);
                problemSet_.get(k).evaluate(solution);
                evaluations++;
                population[k].add(solution);
            }

            updateBestDistances(k);
        }

        // DEBUG: IGD
        pf = new String[taskNum];
        indicators = new ArrayList<>(taskNum);
        for (int k = 0; k < taskNum; k++) {
            pf[k] = "resources/PF/StaticPF/" + problemSet_.get(k).getHType() + "_"
                    + problemSet_.get(k).getNumberOfObjectives() + "D.pf";
            indicators.add(new QualityIndicator(problemSet_.get(k), pf[k]));
        }

        // DEBUG: PLOT
        generations = new ArrayList<>();
        igdPlotValues = new ArrayList<>();
        for (int k = 0; k < taskNum; k++) {
            igdPlotValues.add(new ArrayList<>());
        }

        transferTotalCount = new int[taskNum];
        Arrays.fill(transferTotalCount, 0);
    }

    public void iterate() throws JMException, ClassNotFoundException {
        offspringGeneration(XType);
        environmentSelection();
        generation++;
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
                    int k2 = k;
                    if (PseudoRandom.randDouble() < transferProbability) {
                        while (k2 == k) {
                            k2 = PseudoRandom.randInt(0, taskNum - 1);
                        }
                    }
                    parents[0] = population[k].get(i);
                    parents[1] = population[k2].get(j);

                    Solution child = ((Solution[]) crossover.execute(parents))[PseudoRandom.randInt(0, 1)];
                    if (isMutate)
                        mutation.execute(child);

                    child.setSkillFactor(k);
                    child.setFlag(1);
                    problemSet_.get(k).evaluate(child);
                    evaluations++;
                    offspring[k].add(child);
                }
            } else if (type.equalsIgnoreCase("DE")) {
                for (int i = 0; i < populationSize; i++) {
                    int j1 = i, j2 = i;
                    while (j1 == i && j1 == j2) {
                        j1 = PseudoRandom.randInt(0, populationSize - 1);
                        j2 = PseudoRandom.randInt(0, populationSize - 1);
                    }
                    Solution[] parents = new Solution[3];
                    int k2 = k;
                    if (PseudoRandom.randDouble() < transferProbability) {
                        while (k2 == k) {
                            k2 = PseudoRandom.randInt(0, taskNum - 1);
                        }
                    }
                    parents[0] = population[k2].get(j1);
                    parents[1] = population[k2].get(j2);
                    parents[2] = population[k].get(i);

                    Solution child = (Solution) crossover.execute(new Object[] { population[k].get(i), parents });

                    if (isMutate) {
                        mutation.execute(child);
                    }

                    child.setSkillFactor(k);
                    if (k2 == k){
                        child.setFlag(1);
                    }
                    else {
                        child.setFlag(2);
                        transferTotalCount[k] ++;
                    }
                    problemSet_.get(k).evaluate(child);
                    evaluations++;
                    offspring[k].add(child);
                }
            } else {
                System.out.println("Error: unsupported reproduce type: " + type);
                System.exit(1);
            }
        }
    }

    public void environmentSelection() throws ClassNotFoundException, JMException {
        for (int k = 0; k < taskNum; k++) {
            SolutionSet union = population[k].union(offspring[k]);
            NDSortiong.sort(union, problemSet_, k);

            int normalBetter = 0, transferBetter = 0;
            for (int i = 0; i < populationSize; i++) {
                if (union.get(i).getFlag() == 1)
                    normalBetter ++;
                else if (union.get(i).getFlag() == 2)
                    transferBetter ++;
                
                union.get(i).setFlag(0);
                population[k].replace(i, union.get(i));
                
            }
            
            System.out.println("G_T: " + generation + "_" + k);
            System.out.println("Normal offspring survival Rate: " + (double)normalBetter / population[k].size());
            System.out.println("Transfer offspring survival Rate: " + (double)transferBetter / population[k].size());
            System.out.println("Transfer Success Rate: " + (double)transferBetter / transferTotalCount[k]);
            
            transferTotalCount[k] = 0;

            updateBestDistances(k);

            // 灾变
            // catastrophe(k, 0.1, 50);
        }
    }

    public void catastrophe(int taskID, double survivalRate, int threshold) throws ClassNotFoundException, JMException {
        if (stuckTimes[taskID] >= threshold
                && evaluations < maxEvaluations - 2 * threshold * taskNum * populationSize) {
            // System.out.println(evaluations + ": task " + k +" : reset.");
            stuckTimes[taskID] = 0;
            int[] perm = PseudoRandom.randomPermutation(populationSize, (int) (populationSize * (1 - survivalRate)));
            for (int i = 0; i < perm.length; i++) {
                Solution solution = new Solution(problemSet_);
                solution.setSkillFactor(taskID);
                problemSet_.get(taskID).evaluate(solution);
                // evaluations++;
                population[taskID].replace(perm[i], solution);
            }
        }
    }

    public void updateBestDistances(int taskID) {
        boolean updated = false;
        double avgDistance = 0;
        for (int j = 0; j < population[taskID].size(); j++) {
            double distance = 0;
            for (int i = objStart[taskID]; i <= objEnd[taskID]; i++) {
                distance += Math.pow(population[taskID].get(j).getObjective(i), 2);
            }
            distance = Math.sqrt(distance);
            avgDistance += distance;
        }
        avgDistance /= population[taskID].size();

        if (Math.abs(avgDistance - bestDistances[taskID]) > bestDistances[taskID] * 5e-4) {
            if (avgDistance < bestDistances[taskID]) {
                updated = true;
            }
            bestDistances[taskID] = avgDistance;
        }

        // // DEBUG
        // if (taskID == plotTaskID) {
        // if (updated)
        // System.out.println(avgDistance);
        // else
        // System.out.println("stucking " + stuckTimes[plotTaskID] + " ...");
        // }

        if (updated) {
            stuckTimes[taskID] = 0;
        } else {
            stuckTimes[taskID]++;
        }
    }

    // DEBUG: IGD
    private void calIGD() {
        igd = new double[taskNum];
        for (int k = 0; k < taskNum; k++) {
            igd[k] = indicators.get(k).getIGD(population[k], k);
            igdPlotValues.get(k).add(igd[k]);
        }
        generations.add(generation);
        // System.out.println("Evaluations " + evaluations + ": " +
        // Arrays.toString(igd));
    }

    public double[] getPlotX() {
        return generations.stream().mapToDouble(d -> d).toArray();
    }

    public double[] getPlotX(int maxLength) {
        double[] x = getPlotX();
        return x.length > maxLength ? Arrays.copyOfRange(x, x.length - maxLength, x.length) : x;
    }

    public double[][] getPlotY() {
        double[][] y = new double[taskNum][];
        for (int k = 0; k < taskNum; k++) {
            y[k] = igdPlotValues.get(k).stream().mapToDouble(d -> d).toArray();
        }
        return y;
    }

    public double[][] getPlotY(int maxLength) {
        double[][] y = new double[taskNum][];
        for (int k = 0; k < taskNum; k++) {
            double[] tmp = igdPlotValues.get(k).stream().mapToDouble(d -> d).toArray();
            y[k] = tmp.length > maxLength ? Arrays.copyOfRange(tmp, tmp.length - maxLength, tmp.length) : tmp;
        }
        return y;
    }

    public void initPlot() throws JMException {
        calIGD();
        double[] x = getPlotX();
        double[][] y = getPlotY();
        chartIGD = new XYChartBuilder().title("Generation: " + generation).xAxisTitle("Generation").yAxisTitle("IGD")
                .build();
        // for (int k = 0; k < taskNum; k ++) {
        // chartIGD.addSeries("Problem " + k, x, y[k]);
        // }
        chartIGD.addSeries("Problem " + plotTaskID, x, y[plotTaskID]);
        chartIGD.getStyler().setYAxisLogarithmic(true);

        chartPF = new XYChartBuilder().title("PF: " + generation).xAxisTitle("x").yAxisTitle("y").build();
        chartPF.getStyler().setDefaultSeriesRenderStyle(XYSeriesRenderStyle.Scatter);
        SolutionSet trueParetoFront = new etmo.qualityIndicator.util.MetricsUtil()
                .readNonDominatedSolutionSet(pf[plotTaskID]);
        double[] truePFX = trueParetoFront.getObjectiveVec(problemSet_.get(plotTaskID).getStartObjPos());
        double[] truePFY = trueParetoFront.getObjectiveVec(problemSet_.get(plotTaskID).getEndObjPos());
        chartPF.addSeries("TruePF", truePFX, truePFY);

        double[] PFX = population[plotTaskID].getObjectiveVec(problemSet_.get(plotTaskID).getStartObjPos());
        double[] PFY = population[plotTaskID].getObjectiveVec(problemSet_.get(plotTaskID).getEndObjPos());
        chartPF.addSeries("PF", PFX, PFY);

        chartVar = new XYChartBuilder().title("Var: " + generation).xAxisTitle("Dimension").yAxisTitle("value")
                .width(1024).height(512).build();
        chartVar.getStyler().setLegendVisible(false);
        double[] varX = new double[problemSet_.get(plotTaskID).getNumberOfVariables()];
        for (int i = 0; i < varX.length; i++) {
            varX[i] = i + 1;
        }
        for (int i = 0; i < population[plotTaskID].size(); i++) {
            XYSeries s = chartVar.addSeries("Solution " + i, varX,
                    population[plotTaskID].get(i).getDecisionVariablesInDouble());
            s.setLineColor(XChartSeriesColors.BLUE);
            s.setLineStyle(XChartSeriesLines.SOLID);
            s.setMarker(XChartSeriesMarkers.NONE);
        }

        sw = new SwingWrapper<XYChart>(chartVar);
        sw.displayChart().setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
        sw2 = new SwingWrapper<XYChart>(chartPF);
        sw2.displayChart().setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
        sw3 = new SwingWrapper<XYChart>(chartIGD);
        sw3.displayChart().setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
    }

    public void updatePlot() {
        calIGD();
        double[] x = getPlotX();
        double[][] y = getPlotY();

        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                chartIGD.setTitle("Generation: " + generation);
                // for (int k = 0; k < taskNum; k ++){
                // chartIGD.updateXYSeries("Problem " + k, x, y[k], null);
                // }
                chartIGD.updateXYSeries("Problem " + plotTaskID, x, y[plotTaskID], null);

                chartPF.setTitle("PF: " + generation);
                double[] PFX = population[plotTaskID].getObjectiveVec(problemSet_.get(plotTaskID).getStartObjPos());
                double[] PFY = population[plotTaskID].getObjectiveVec(problemSet_.get(plotTaskID).getEndObjPos());
                chartPF.updateXYSeries("PF", PFX, PFY, null);

                chartVar.setTitle("Var: " + generation);
                double[] varX = new double[problemSet_.get(plotTaskID).getNumberOfVariables()];
                for (int i = 0; i < varX.length; i++) {
                    varX[i] = i + 1;
                }
                for (int i = 0; i < population[0].size(); i++) {
                    try {
                        chartVar.updateXYSeries("Solution " + i, varX,
                                population[plotTaskID].get(i).getDecisionVariablesInDouble(), null);
                    } catch (JMException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }

                sw.repaintChart();
                sw2.repaintChart();
                sw3.repaintChart();
            }
        });
    }

    public void endPlot() {
        try {
            BitmapEncoder.saveBitmap(chartIGD, "./figs/" + problemSet_.get(0).getName(), BitmapFormat.PNG);
            // VectorGraphicsEncoder.saveVectorGraphic(chart, "./figs/" +
            // problemSet_.get(0).getName(), VectorGraphicsFormat.PDF);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}