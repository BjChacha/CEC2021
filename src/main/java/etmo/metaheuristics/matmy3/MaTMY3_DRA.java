package etmo.metaheuristics.matmy3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.WindowConstants;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.SingularMatrixException;
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
import etmo.util.math.Random;
import etmo.util.math.Vector;
import etmo.util.sorting.NDSortiong;
import etmo.util.sorting.SortingIdx;

public class MaTMY3_DRA extends MtoAlgorithm {
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

    private double[] decisionDistances;
    private double[] objectiveDistances;
    private double[] improvements;
    private double[] resources;

    private double mutationProbability;
    private double transferProbability;
    private double[][] tP;

    boolean isMutate;


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
    XYChart chartDRA;
    List<XYChart> charts;
    SwingWrapper<XYChart> sw;
    SwingWrapper<XYChart> sw2;
    SwingWrapper<XYChart> sw3;
    SwingWrapper<XYChart> sw4;
    List<Integer> generations;
    List<List<Double>> igdPlotValues;

    public MaTMY3_DRA(ProblemSet problemSet) {
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initState();
        if (isPlot)
            initPlot();
        while (evaluations < maxEvaluations) {
            iterate();
            if (isPlot) {
                updatePlot();
            }
            resetFlag();
        }
        // if (isPlot)
        // endPlot();
        // System.out.println(igdPlotValues.get(0).toString());

        System.out.println(Arrays.toString(resources));
        return population;
    }

    void initState() throws JMException, ClassNotFoundException {
        generation = 0;
        evaluations = 0;
        taskNum = problemSet_.size();
        varNum = problemSet_.getMaxDimension();

        maxEvaluations = (Integer) this.getInputParameter("maxEvaluations");
        populationSize = (Integer) this.getInputParameter("populationSize");

        // // DEBUG partly run
        // maxEvaluations = maxEvaluations / taskNum * 2;
        // // problemSet_ = problemSet_.getTask(0);
        // taskNum = 2;

        XType = (String) this.getInputParameter("XType");
        TXType = (String) this.getInputParameter("TXType");
        isPlot = (Boolean) this.getInputParameter("isPlot");
        isMutate = (Boolean) this.getInputParameter("isMutate");

        transferProbability = (Double) this.getInputParameter("transferProbability");

        // DEBUG: IGD PLOTTING
        plotTaskID = (Integer) this.getInputParameter("plotTaskID");

        objStart = new int[taskNum];
        objEnd = new int[taskNum];
        tP = new double[taskNum][taskNum];

        DECrossover = operators_.get("DECrossover");
        SBXCrossover = operators_.get("SBXCrossover");
        mutation = operators_.get("mutation");

        decisionDistances = new double[taskNum];
        objectiveDistances = new double[taskNum];
        improvements = new double[taskNum];
        resources = new double[taskNum];

        mutationProbability = 0.5;

        population = new SolutionSet[taskNum];
        offspring = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
            objStart[k] = problemSet_.get(k).getStartObjPos();
            objEnd[k] = problemSet_.get(k).getEndObjPos();

            Arrays.fill(tP[k], transferProbability);

            population[k] = new SolutionSet(populationSize);
            offspring[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++) {
                Solution solution = new Solution(problemSet_);
                solution.setSkillFactor(k);
                evaluate(solution, k);
                population[k].add(solution);
            }
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
        resourceAllocation(1);
        updateState();
    }

    void offspringGeneration() throws JMException {
        for (int k = 0; k < taskNum; k++) {
            offspringGeneration(k);
        }
    }

    void offspringGeneration(int taskID) throws JMException {
        int k = taskID;
        offspring[k].clear();
        Solution child;

        int[] perm = PseudoRandom.randomPermutation(population[k].size(), population[k].size());
        for (int i = 0; i < populationSize && evaluations < maxEvaluations; i ++) {
            child = null;
            if (PseudoRandom.randDouble() < transferProbability) {
                int k2 = getAssistTaskID(k);
                child = transferGenerating(k, k2, perm[i], TXType);
            } else {
                child = evolutionaryGenerating(k, perm[i], XType);
            }

            evaluate(child, k);
            offspring[k].add(child);
        }
    }

    void environmentSelection() throws ClassNotFoundException, JMException {
        for (int k = 0; k < taskNum; k++) {
            environmentSelection(k);
        }
    }

    void environmentSelection(int taskID) throws ClassNotFoundException, JMException {
        int k = taskID;
        SolutionSet union = population[k].union(offspring[k]);
        NDSortiong.sort(union, problemSet_, k);

        double[][] oldDecisionPosition = population[k].getMat();
        double[][] oldObjectivePosition = population[k].writeObjectivesToMatrix(objStart[k], objEnd[k]);

        for (int i = 0; i < populationSize; i++) {
            // union.get(i).setFlag(0);
            union.get(i).setSkillFactor(k);
            population[k].replace(i, union.get(i));
        }

        double[][] newDecisionPosition = population[k].getMat();
        double[][] newObjectivePosition = population[k].writeObjectivesToMatrix(objStart[k], objEnd[k]);

        decisionDistances[k] = Distance.getWassersteinDistance(
            oldDecisionPosition,
            newDecisionPosition);
        objectiveDistances[k] = Distance.getWassersteinDistance(
            oldObjectivePosition, 
            newObjectivePosition);
        improvements[k] = objectiveDistances[k] / decisionDistances[k];
        resources[k] += 1.0;
    }

    void resourceAllocation(int num) throws JMException, ClassNotFoundException {
        int[] order = SortingIdx.sort(Vector.vecElemDiv(
            objectiveDistances, decisionDistances), false);
        
        for (int i = 0; i < num; i ++) {
            offspringGeneration(order[i]);
            environmentSelection(order[i]);
        }
    }

    Solution transferGenerating(int taskID, int assistTaskID, int i, String type) throws JMException {
        int j = PseudoRandom.randInt(0, population[assistTaskID].size() - 1);
        
        // // explicit
        // Solution child = null;
        // child = new Solution(population[assistTaskID].get(j));

        // SBX implicit
        Solution[] parents = new Solution[2];
        parents[0] = population[taskID].get(i);
        parents[1] = population[assistTaskID].get(j);
        Solution child = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];
        mutateIndividual(taskID, child);
        
        // // combine with Gaussian Distribution
        // Solution[] parents = new Solution[2];
        // parents[0] = new Solution(population[taskID].get(i));
        // double[] sampleFeatures = Probability.sampleByNorm(means[taskID], stds[taskID]);
        // Vector.vecClip_(sampleFeatures, 0.0, 1.0);
        // parents[0].setDecisionVariables(sampleFeatures);
        // parents[1] = offsetIndividual(population[assistTaskID].get(j), assistTaskID, taskID);
        // Solution child = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];
        // mutateIndividual(taskID, child);

        // // combine with Gaussian Distribution advanced
        // Solution[] parents = new Solution[2];
        // parents[0] = new Solution(population[taskID].get(i));
        // double[] sampleFeatures = Probability.sampleByNorm(means[taskID], stds[taskID]);
        // Vector.vecClip_(sampleFeatures, 0.0, 1.0);
        // parents[0].setDecisionVariables(sampleFeatures);
        // parents[1] = offsetIndividual(evolutionaryGenerating(assistTaskID, j, XType), assistTaskID, taskID);
        // Solution child = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];
        // mutateIndividual(taskID, child);
        
        // // explicit: sample then transform
        // Solution child = new Solution(population[assistTaskID].get(j));
        // double[] sampleFeatures = Probability.sampleByNorm(means[assistTaskID], stds[assistTaskID]);
        // child.setDecisionVariables(sampleFeatures);
        // child = offsetIndividual(child, assistTaskID, taskID);
        child.setFlag(1);
        return child;
    }

    Solution evolutionaryGenerating(int taskID, int i, String type) throws JMException {
        Solution child = null;
        if (type.equalsIgnoreCase("SBX")) {
            child = SBXChildGenerating(taskID, taskID, i);
        } else if (type.equalsIgnoreCase("DE")) {
            child = DEChildGenerating(taskID, taskID, i);
        }
        else {
            System.out.println("Error: unsupported reproduce type: " + type);
            System.exit(1);
        }
        child.setFlag(2);
        return child;
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

        return child;
    }

    void mutateIndividual(int taskID, Solution individual) throws JMException {
        // if (isMutate)
        //     mutation.execute(individual);

        // if (PseudoRandom.randDouble() < stuckTimes[taskID] * 0.15)
        //     mutation.execute(individual);

        if (PseudoRandom.randDouble() < mutationProbability) {
            mutation.execute(individual);
            // individual.setFlag(1);
        }
    }

    int getAssistTaskID(int taskID) throws JMException {
        int assistTaskID = taskID;

        // // random
        // while (assistTaskID == taskID) {
        //     assistTaskID = PseudoRandom.randInt(0, taskNum - 1);
        // }

        assistTaskID = Random.rouletteWheel(resources);

        // double[] scores = new double[taskNum];
        // for (int i = 0; i < scores.length; i ++) {
        //     scores[i] = 1.0 / distances[taskID][i];
        // }
        // assistTaskID = Random.rouletteWheel(scores, taskID);

        return assistTaskID;
    }
        
    void evaluate(Solution solution, int taskID) throws JMException {
        // solution.setSkillFactor(taskID);
        problemSet_.get(taskID).evaluate(solution);
        evaluations ++;
    }

    void updateState() {
        generation ++;
        resetFlag();
    }

    void resetFlag() {
        for (int k = 0; k < taskNum; k++) {
            for (int i = 0; i < population[k].size(); i++) {
                population[k].get(i).setFlag(0);
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

    double[][][] getPFPlotData(int taskID) {
        SolutionSet[] subPops = new SolutionSet[3];
        for (int i = 0; i < 3; i++) {
            subPops[i] = new SolutionSet();
        }
        for (int i = 0; i < population[taskID].size(); i++) {
            int flag = population[taskID].get(i).getFlag();
            subPops[flag].add(population[taskID].get(i));
        }
        double[][][] data = new double[subPops.length][2][];
        for (int i = 0; i < data.length; i ++) {
            data[i][0] = subPops[i].getObjectiveVec(objStart[taskID]);
            data[i][1] = subPops[i].getObjectiveVec(objEnd[taskID]);
        }
        return data;
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
        XYSeries s = chartPF.addSeries("TruePF", truePFX, truePFY);
        s.setMarkerColor(XChartSeriesColors.BLACK);

        double[] PFX = population[plotTaskID].getObjectiveVec(problemSet_.get(plotTaskID).getStartObjPos());
        double[] PFY = population[plotTaskID].getObjectiveVec(problemSet_.get(plotTaskID).getEndObjPos());
        s = chartPF.addSeries("Parents", PFX, PFY);
        s.setMarker(XChartSeriesMarkers.CIRCLE);

        chartVar = new XYChartBuilder().title("Var: " + generation).xAxisTitle("Dimension").yAxisTitle("value")
                .width(1024).height(512).build();
        chartVar.getStyler().setLegendVisible(false);
        double[] varX = new double[varNum];
        for (int i = 0; i < varX.length; i++) {
            varX[i] = i + 1;
        }
        for (int i = 0; i < population[plotTaskID].size(); i++) {
            s = chartVar.addSeries("Solution " + i, varX,
                    population[plotTaskID].get(i).getDecisionVariablesInDouble());
            s.setLineColor(XChartSeriesColors.BLUE);
            s.setLineStyle(XChartSeriesLines.SOLID);
            s.setMarker(XChartSeriesMarkers.NONE);
        }

        double[] varXX = new double[taskNum];
        for (int i = 0; i < varX.length; i++) {
            varXX[i] = i + 1;
        }
        double[] varYY = Vector.vecElemDiv(objectiveDistances, decisionDistances);
        chartDRA = new XYChartBuilder().title("DRA: " + generation).build();
        chartDRA.getStyler().setLegendVisible(false);
        s = chartDRA.addSeries("DRA", varXX, varYY);
        s.setLineColor(XChartSeriesColors.BLUE);
        s.setLineStyle(XChartSeriesLines.SOLID);
        s.setMarker(XChartSeriesMarkers.NONE);
        sw4 = new SwingWrapper<XYChart>(chartDRA);
        sw4.displayChart().setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);

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

        double[][][] data = getPFPlotData(plotTaskID);
        String[] names = new String[] {"Parents", "tChlid", "nChild"};
        
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                chartIGD.setTitle("Generation: " + generation);
                // for (int k = 0; k < taskNum; k ++){
                // chartIGD.updateXYSeries("Problem " + k, x, y[k], null);
                // }
                chartIGD.updateXYSeries("Problem " + plotTaskID, x, y[plotTaskID], null);

                chartPF.setTitle("PF: " + generation);
                // double[] PFX = population[plotTaskID].getObjectiveVec(problemSet_.get(plotTaskID).getStartObjPos());
                // double[] PFY = population[plotTaskID].getObjectiveVec(problemSet_.get(plotTaskID).getEndObjPos());
                // chartPF.updateXYSeries("PF", PFX, PFY, null);
                for (int i = 0; i < 3; i++) {
                    if (data[i][0].length > 0) {
                        if (chartPF.getSeriesMap().containsKey(names[i])) {
                            chartPF.updateXYSeries(names[i], data[i][0], data[i][1], null);
                        } else {
                            XYSeries ss = chartPF.addSeries(names[i], data[i][0], data[i][1]);
                            ss.setMarker(XChartSeriesMarkers.CIRCLE);
                        }
                    }
                }

                chartVar.setTitle("Var: " + generation);
                double[] varX = new double[varNum];
                for (int i = 0; i < varX.length; i++) {
                    varX[i] = i + 1;
                }
                for (int i = 0; i < population[0].size(); i++) {
                    try {
                        chartVar.updateXYSeries("Solution " + i, varX,
                                population[plotTaskID].get(i).getDecisionVariablesInDouble(), null);
                    } catch (JMException e) {
                        e.printStackTrace();
                    }
                }
                
                double[] varXX = new double[taskNum];
                Arrays.parallelSetAll(varXX, i -> i + 1);
                chartDRA.setTitle("DRA: " + generation);
                chartDRA.updateXYSeries("DRA", varXX, Vector.vecElemDiv(
                    objectiveDistances, decisionDistances), null);
                sw4.repaintChart();

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
