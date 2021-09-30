package etmo.metaheuristics.matmy3;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
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

import breeze.linalg.kron;
import etmo.core.MtoAlgorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.metaheuristics.matmy3.models.AbstractDistribution;
import etmo.metaheuristics.matmy3.models.GaussianDistribution;
import etmo.metaheuristics.matmy3.models.MultiVarGaussian;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.Configuration;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.math.Distance;
import etmo.util.math.Matrix;
import etmo.util.math.Probability;
import etmo.util.math.Random;
import etmo.util.math.Vector;
import etmo.util.sorting.NDSortiong;
import spire.optional.intervalValuePartialOrder;

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

    private double mutationProbability;
    private double transferProbability;
    private double[] tP;

    boolean isMutate;

    double[][] elitePositions;
    double[][] eliteMovements;

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
    XYChart chartCluster;
    List<XYChart> charts;
    SwingWrapper<XYChart> sw;
    SwingWrapper<XYChart> sw2;
    SwingWrapper<XYChart> sw3;
    SwingWrapper<XYChart> sw4;
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
            if (isPlot) {
                updatePlot();
            }
            resetFlag();
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
 
        DECrossover = operators_.get("DECrossover");
        SBXCrossover = operators_.get("SBXCrossover");
        mutation = operators_.get("mutation");

        tP = new double[taskNum];
        Arrays.fill(tP, transferProbability);
        mutationProbability = 0.5;

        objStart = new int[taskNum];
        objEnd = new int[taskNum];

        elitePositions = new double[taskNum][varNum];
        eliteMovements = new double[taskNum][varNum];

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
                evaluate(solution, k);
                population[k].add(solution);
            }
            NDSortiong.sort(population[k], problemSet_, k);
            computeElitePosition(k);
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

    void computeElitePosition(int taskID) {
        SolutionSet eliteSet = new SolutionSet();
        for (int i = 0; i < population[taskID].size(); i++) {
            if (population[taskID].get(i).getRank() == 0)
                eliteSet.add(population[taskID].get(i));
            else
                break;
        }
        elitePositions[taskID] = eliteSet.getMean();
    }

    void iterate() throws JMException, ClassNotFoundException {
        offspringGeneration();
        environmentSelection();
        // writePopulationVariablesMatrix(plotTaskID, generation);
        generation++;
    }

    void offspringGeneration() throws JMException {
        for (int k = 0; k < taskNum; k++) {
            offspring[k].clear();
            Solution child;

            int[] perm = PseudoRandom.randomPermutation(population[k].size(), population[k].size());
            for (int i = 0; i < populationSize && !offspring[k].isFull(); i ++) {
                child = null;
                if (PseudoRandom.randDouble() < tP[k]) {
                    int k2 = getAssistTaskID(k);
                    child = transferGenerating(k, k2, perm[i], TXType);
                } else {
                    child = evolutionaryGenerating(k, perm[i], XType);
                }
                evaluate(child, k);
                offspring[k].add(child);
            }
        }
    }

    Solution generatingChild(int taskID, int i) throws JMException {
        Solution child = null;
        int[] perm = PseudoRandom.randomPermutation(population[taskID].size(), population[taskID].size());
        if (PseudoRandom.randDouble() < transferProbability) {
            int assistTaskID = getAssistTaskID(taskID);
            child = transferGenerating(taskID, assistTaskID, perm[i], TXType);
        } else {
            child = evolutionaryGenerating(taskID, perm[i], XType);
        }
        return child;
    }

    Solution transferGenerating(int taskID, int assistTaskID, int i, String type) throws JMException {
        int j = PseudoRandom.randInt(0, population[assistTaskID].size() - 1);
        Solution child = null;
        
        // explicit
        child = new Solution(population[assistTaskID].get(j));
        mutateIndividual(taskID, child);
        
        // // SBX implicit
        // Solution[] parents = new Solution[2];
        // parents[0] = population[taskID].get(i);
        // parents[1] = population[assistTaskID].get(j);
        // child = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];
        // mutateIndividual(taskID, child);
        
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

    void environmentSelection() throws ClassNotFoundException, JMException {
        for (int k = 0; k < taskNum; k++) {
            SolutionSet union = population[k].union(offspring[k]);
            NDSortiong.sort(union, problemSet_, k);

            double[] oldPosition = elitePositions[k].clone();

            for (int i = 0; i < populationSize; i++) {
                union.get(i).setSkillFactor(k);
                population[k].replace(i, union.get(i));
            }

            computeElitePosition(k);
            eliteMovements[k] = Vector.vecSub(elitePositions[k], oldPosition);
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
    
    int getAssistTaskID(int taskID) {
        int assistTaskID = taskID;

        // // random
        // while (assistTaskID == taskID)
        //     assistTaskID = PseudoRandom.randInt(0, taskNum - 1);
        double factor = 1.0;
        double norm = Vector.vecModule(eliteMovements[taskID]);
        double[] score = new double[taskNum];
        boolean found = false;
        while (factor <= 3.0) {
            for (int k = 0; k < taskNum; k ++) {
                if (k == taskID) continue;
                if (Vector.vecModule(Vector.vecSub(elitePositions[taskID], elitePositions[k])) > norm) continue;
                score[k] = 1 - Distance.getCosineSimilarity(elitePositions[taskID], elitePositions[k]);
                found = true;
            }
            if (found) break;
            factor += 1.0;
        }
        assistTaskID = Random.rouletteWheel(score);

        return assistTaskID;
    }

    void evaluate(Solution solution, int taskID) throws JMException {
        // solution.setSkillFactor(taskID);
        problemSet_.get(taskID).evaluate(solution);
        evaluations ++;
    }

    void resetFlag() {
        for (int k = 0; k < taskNum; k++) {
            for (int i = 0; i < population[k].size(); i++) {
                population[k].get(i).setFlag(0);
            }
        }
    }

    void writePopulationVariablesMatrix(int taskID, int generation) throws JMException {
        String algoName = "Gaussian_exp"; 
        String folderPath = "./data/variables";
        File folder = new File(folderPath);
        if (!folder.exists()){
            folder.mkdirs();
        }
        String filePath = folderPath + "/" + algoName + "_" + generation + ".txt";
        double[][] data = population[taskID].getMat();
        try {
            FileOutputStream fos = new FileOutputStream(filePath);
            OutputStreamWriter osw = new OutputStreamWriter(fos);
            BufferedWriter bw = new BufferedWriter(osw);

            for (double[] line: data) {
                String sLine = Arrays.toString(line)
                    .replace("[", "")
                    .replace("]", "")
                    .replace(",", "")
                    .strip();
                bw.write(sLine);
                bw.newLine();
            }
            bw.close();
        } catch (IOException e) {
            Configuration.logger_.severe("Error acceding to the file");
            e.printStackTrace();
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

        // // debug: tSNE performance
        // SolutionSet total = new SolutionSet(populationSize * taskNum / 2);
        // for (int k = 0; k < taskNum; k++) {
        //     for (int i = 0; i < population[k].size() / 2; i++){
        //         total.add(population[k].get(i));
        //     }
        // }
        // TSNE tSNE = null;
        // try {
        //     tSNE = new TSNE(total.getMat(), 2);
        // } catch (JMException e) {
        //     // TODO Auto-generated catch block
        //     e.printStackTrace();
        // }
        // double[][] coordinates = tSNE.coordinates;
        // double[][] X = new double[taskNum][populationSize / 2];
        // double[][] Y = new double[taskNum][populationSize / 2];
        // for (int i = 0; i < coordinates.length; i++) {
        //     X[i/populationSize][i%(populationSize/2)] = coordinates[i][0];
        //     Y[i/populationSize][i%(populationSize/2)] = coordinates[i][1];
        // }

        // chartCluster = new XYChartBuilder().title("t-SNE: " + generation).xAxisTitle("x").yAxisTitle("y").build();
        // chartCluster.getStyler().setDefaultSeriesRenderStyle(XYSeriesRenderStyle.Scatter);

        // for (int k = 0; k < taskNum; k++) {
        //     chartCluster.addSeries("P " + k, X[k], Y[k]);
        // }

        // sw4 = new SwingWrapper<XYChart>(chartCluster);
        // sw4.displayChart().setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);

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
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }

                // if (generation % 20 == 0) {
                //     SolutionSet total = new SolutionSet(populationSize * taskNum / 2);
                //     for (int k = 0; k < taskNum; k++) {
                //         for (int i = 0; i < population[k].size() / 2; i++){
                //             total.add(population[k].get(i));
                //         }
                //     }
                //     TSNE tSNE = null;
                //     try {
                //         tSNE = new TSNE(total.getMat(), 2);
                //     } catch (JMException e) {
                //         // TODO Auto-generated catch block
                //         e.printStackTrace();
                //     }
                //     double[][] coordinates = tSNE.coordinates;
                //     double[][] X = new double[taskNum][populationSize / 2];
                //     double[][] Y = new double[taskNum][populationSize / 2];
                //     for (int i = 0; i < coordinates.length; i++) {
                //         X[i/populationSize][i%(populationSize/2)] = coordinates[i][0];
                //         Y[i/populationSize][i%(populationSize/2)] = coordinates[i][1];
                //     }
                //     chartCluster.setTitle("t-SNE: " + generation);
                //     for (int k = 0; k < taskNum; k++) {
                //         chartCluster.updateXYSeries("P " + k, X[k], Y[k], null);
                //     }
                //     sw4.repaintChart();
                // }
                
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
