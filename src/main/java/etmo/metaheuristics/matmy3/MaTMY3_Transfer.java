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
import etmo.util.math.Distance;
import etmo.util.math.Probability;
import etmo.util.math.Random;
import etmo.util.math.Vector;
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
    private String TXType;
    private Operator DECrossover;
    private Operator SBXCrossover;
    private Operator mutation;

    private int[] stuckTimes;
    private double[] bestDistances;
    
    private double transferProbability;
    private double[][] tP;
    private double[] transferBaseline;
    private double[][] distances;

    boolean isMutate;

    double[] center;
    double[][] bias;
    double[][] momentum;
    // momentum effect factor
    double alpha;
    // momentum update factor
    double beta;

    // TODO: duplicate with bias
    double[][] means;
    double[][] stds;

    // DEBUG: IGD
    String[] pf;
    List<QualityIndicator> indicators;
    double[] igd;

    // DEBUG: PLOT
    boolean isPlot;
    int plotTaskID;
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
        // System.out.println(igdPlotValues.get(0).toString());
        return population;
    }

    void initState() throws JMException, ClassNotFoundException {
        generation = 0;
        evaluations = 0;
        taskNum = problemSet_.size();
        varNum = problemSet_.getMaxDimension();

        maxEvaluations = (Integer) this.getInputParameter("maxEvaluations");
        populationSize = (Integer) this.getInputParameter("populationSize");

        XType = (String) this.getInputParameter("XType");
        TXType = (String) this.getInputParameter("TXType");
        isPlot = (Boolean) this.getInputParameter("isPlot");
        isMutate = (Boolean) this.getInputParameter("isMutate");

        transferProbability = (Double) this.getInputParameter("transferProbability");

        tP = new double[taskNum][taskNum];

        // DEBUG: IGD PLOTTING
        plotTaskID = (Integer) this.getInputParameter("plotTaskID");
 
        DECrossover = operators_.get("DECrossover");
        SBXCrossover = operators_.get("SBXCrossover");
        mutation = operators_.get("mutation");

        objStart = new int[taskNum];
        objEnd = new int[taskNum];
        bestDistances = new double[taskNum];
        stuckTimes = new int[taskNum];
        Arrays.fill(bestDistances, Double.MAX_VALUE);
        Arrays.fill(stuckTimes, 0);

        center = new double[varNum];
        Arrays.fill(center, 0.5);
        bias = new double[taskNum][varNum];
        momentum = new double[taskNum][varNum];
        alpha = 0.2;
        beta = 0.7;
        
        distances = new double[taskNum][taskNum];

        means = new double[taskNum][varNum];
        stds = new double[taskNum][varNum];

        population = new SolutionSet[taskNum];
        offspring = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
            objStart[k] = problemSet_.get(k).getStartObjPos();
            objEnd[k] = problemSet_.get(k).getEndObjPos();

            Arrays.fill(distances[k], 0);
            Arrays.fill(means[k], 0);
            Arrays.fill(stds[k], 0);

            Arrays.fill(bias[k], 0);
            Arrays.fill(momentum[k], 0);
            Arrays.fill(tP[k], transferProbability);

            population[k] = new SolutionSet(populationSize);
            offspring[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++) {
                Solution solution = new Solution(problemSet_);
                evaluate(solution, k);
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

    void iterate() throws JMException, ClassNotFoundException {
        offspringGeneration();
        environmentSelection();
        generation++;
    }

    void offspringGeneration() throws JMException {
        if (generation % 10 == 1)
            updateDistances();
        
        for (int k = 0; k < taskNum; k++) {
            offspring[k].clear();
            int[] perm = PseudoRandom.randomPermutation(population[k].size(), population[k].size());
            for (int i = 0; i < populationSize; i ++) {
                if (PseudoRandom.randDouble() < transferProbability) {
                    int k2 = getAssistTaskID(k);
                    transferGenerating(k, k2, perm[i], TXType);
                } else {
                    evolutionaryGenerating(k, perm[i], XType);
                }
            }
        }
    }

    void transferGenerating(int taskID, int assistTaskID, int i, String type) throws JMException {
        int j = PseudoRandom.randInt(0, population[assistTaskID].size() - 1);
        
        // // explicit
        // Solution child = offsetIndividual(population[assistTaskID].get(j), assistTaskID, taskID);
        
        // SBX implicit
        // Solution[] parents = new Solution[2];
        // parents[0] = population[taskID].get(i);
        // parents[1] = offsetIndividual(population[assistTaskID].get(j), assistTaskID, taskID);
        // Solution child = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];
        // mutateIndividual(taskID, child);
        
        // combine with Gaussian Distribution
        Solution[] parents = new Solution[2];
        parents[0] = new Solution(population[taskID].get(i));
        parents[0].setDecisionVariables(Probability.sampleByNorm(means[taskID], stds[taskID]));
        parents[1] = offsetIndividual(population[assistTaskID].get(j), assistTaskID, taskID);
        Solution child = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];
        mutateIndividual(taskID, child);
        
        evaluate(child, taskID);
        offspring[taskID].add(child);
    }

    void evolutionaryGenerating(int taskID, int i, String type) throws JMException {
        Solution child;
        if (type.equalsIgnoreCase("SBX")) {
            child = SBXChildGenerating(taskID, taskID, i);
            offspring[taskID].add(child);
        } else if (type.equalsIgnoreCase("DE")) {
            child = DEChildGenerating(taskID, taskID, i);
            offspring[taskID].add(child);
        }
        else {
            System.out.println("Error: unsupported reproduce type: " + type);
            System.exit(1);
        }
    }

    void environmentSelection() throws ClassNotFoundException, JMException {
        for (int k = 0; k < taskNum; k++) {
            SolutionSet union = population[k].union(offspring[k]);
            NDSortiong.sort(union, problemSet_, k);

            // int normalBetter = 0, transferBetter = 0;
            for (int i = 0; i < populationSize; i++) {
                // if (union.get(i).getFlag() == 1)
                //     normalBetter ++;
                // else if (union.get(i).getFlag() == 2)
                //     transferBetter ++;
                
                union.get(i).setFlag(0);
                population[k].replace(i, union.get(i));
            }
            
            // System.out.println("Task " + k + ": " + Arrays.toString(bias[k]));
            // System.out.println("G_T: " + generation + "_" + k);
            // System.out.println("Mutation offspring survival Rate: " + (double)normalBetter / population[k].size());
            // System.out.println("Normal offspring survival Rate: " + (double)normalBetter / population[k].size());
            // System.out.println("Transfer offspring survival Rate: " + (double)transferBetter / population[k].size());
            // System.out.println("Transfer Success Rate: " + (double)transferBetter / transferTotalCount[k]);
            
            // transferTotalCount[k] = 0;
            
            updateBias(k);
            updateBestDistances(k);

            // 灾变
            catastrophe(k, 0.5, 75);
        }
    }

    Solution DEChildGenerating(int taskID, int assistTaskID, int i) throws JMException {
        int j1 = i, j2 = i;
        while (j1 == i && j1 == j2) {
            j1 = PseudoRandom.randInt(0, populationSize - 1);
            j2 = PseudoRandom.randInt(0, populationSize - 1);
        }
        Solution[] parents = new Solution[3];
        parents[0] = population[assistTaskID].get(j1);
        parents[1] = population[assistTaskID].get(j2);
        parents[2] = population[taskID].get(i);

        Solution child = (Solution) DECrossover.execute(new Object[] { population[taskID].get(i), parents });
        mutateIndividual(taskID, child);
        evaluate(child, taskID);

        return child;
    }

    Solution SBXChildGenerating(int taskID, int assistTaskID, int i) throws JMException {
        int j = i;
        while (j == i)
            j = PseudoRandom.randInt(0, populationSize - 1);
        Solution[] parents = new Solution[2];
        parents[0] = population[taskID].get(i);
        parents[1] = population[assistTaskID].get(j);

        Solution child = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];
        mutateIndividual(taskID, child);
        evaluate(child, taskID);

        return child;
    }

    void mutateIndividual(int taskID, Solution individual) throws JMException {
        // if (isMutate)
        //     mutation.execute(individual);

        // if (PseudoRandom.randDouble() < stuckTimes[taskID] * 0.15)
        //     mutation.execute(individual);

        if (PseudoRandom.randDouble() < 0.5) {
            mutation.execute(individual);
            // individual.setFlag(1);
        }
    }

    void updateBias(int taskID) throws JMException {
        // update bias distance
        double[] oldBias = bias[taskID].clone();
        Arrays.fill(bias[taskID], 0);
        for (int i = 0; i < populationSize / 2; i ++) {
            Vector.vecAdd_(bias[taskID], Vector.vecSub(population[taskID].get(i).getDecisionVariablesInDouble(), center));
        }
        Vector.vecElemMul_(bias[taskID], 2.0 / populationSize);
        
        // update momentum
        if (generation >= 2) {
            // momentum[k] = Vector.vecSub(bias[k], oldBias);
            double[] newMomentum = Vector.vecSub(bias[taskID], oldBias);
            momentum[taskID] = Vector.vecAdd(
                Vector.vecElemMul(momentum[taskID], 1 - beta),
                Vector.vecElemMul(newMomentum, beta)
            );
        }
    }

    void updateDistributions(double partition) {
        for (int k = 0; k < taskNum; k++) {
            int size = (int)(population[k].size() * partition);
            SolutionSet tmpSet = new SolutionSet(size);
            for (int i = 0; i < size; i++) {
                tmpSet.add(population[k].get(i));
            }
            means[k] = tmpSet.getMean();
            stds[k] = tmpSet.getStd();
        }
    }

    int getAssistTaskID(int taskID) throws JMException {
        int assistTaskID = taskID;

        // while (assistTaskID == taskID) {
        //     assistTaskID = PseudoRandom.randInt(0, taskNum - 1);
        // }

        assistTaskID = Random.rouletteWheel(distances[taskID], taskID);

        return assistTaskID;
    }
  
    Solution offsetIndividual(Solution solution, int fromTaskID, int toTaskID) throws JMException {
        double[] features = solution.getDecisionVariablesInDouble();
        Vector.vecSub_(features, bias[fromTaskID]);
        Vector.vecAdd_(features, bias[toTaskID]);
        Vector.vecAdd_(features, Vector.vecElemMul(momentum[fromTaskID], alpha));
        Vector.vecClip_(features, 0.0, 1.0);
        Solution newSolution = new Solution(solution);
        newSolution.setDecisionVariables(features);
        return newSolution;
    }
      
    void updateBestDistances(int taskID) {
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
        
        if (Math.abs(avgDistance - bestDistances[taskID]) > bestDistances[taskID] * 1e-3) {
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
    
    void updateDistances() throws JMException {
        for (int k = 0; k < taskNum; k++) {
            for (int kk = k + 1; kk < taskNum; kk++) {
                distances[k][kk] = distances[kk][k] = 1.0 / Distance.getWassersteinDistance(
                    population[kk].getMat(),
                    population[k].getMat());
            }
        }
    }
        
    void evaluate(Solution solution, int taskID) throws JMException {
        solution.setSkillFactor(taskID);
        problemSet_.get(taskID).evaluate(solution);
        evaluations ++;
    }

    void catastrophe(int taskID, double survivalRate, int threshold) throws ClassNotFoundException, JMException {
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
    
    // DEBUG: IGD
    void calIGD() {
        igd = new double[taskNum];
        for (int k = 0; k < taskNum; k++) {
            igd[k] = indicators.get(k).getIGD(population[k], k);
            igdPlotValues.get(k).add(igd[k]);
        }
        generations.add(generation);
        // System.out.println("Evaluations " + evaluations + ": " +
        // Arrays.toString(igd));
    }
        
    double[] getPlotX() {
        return generations.stream().mapToDouble(d -> d).toArray();
    }
    
    double[] getPlotX(int maxLength) {
        double[] x = getPlotX();
        return x.length > maxLength ? Arrays.copyOfRange(x, x.length - maxLength, x.length) : x;
    }
    
    double[][] getPlotY() {
        double[][] y = new double[taskNum][];
        for (int k = 0; k < taskNum; k++) {
            y[k] = igdPlotValues.get(k).stream().mapToDouble(d -> d).toArray();
        }
        return y;
    }

    double[][] getPlotY(int maxLength) {
        double[][] y = new double[taskNum][];
        for (int k = 0; k < taskNum; k++) {
            double[] tmp = igdPlotValues.get(k).stream().mapToDouble(d -> d).toArray();
            y[k] = tmp.length > maxLength ? Arrays.copyOfRange(tmp, tmp.length - maxLength, tmp.length) : tmp;
        }
        return y;
    }

    void initPlot() throws JMException {
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
        double[] truePFX = trueParetoFront.getObjectiveVec(0);
        double[] truePFY = trueParetoFront.getObjectiveVec(1);
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

    void updatePlot() {
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

    void endPlot() {
        try {
            BitmapEncoder.saveBitmap(chartIGD, "./figs/" + problemSet_.get(0).getName(), BitmapFormat.PNG);
            // VectorGraphicsEncoder.saveVectorGraphic(chart, "./figs/" +
            // problemSet_.get(0).getName(), VectorGraphicsFormat.PDF);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
