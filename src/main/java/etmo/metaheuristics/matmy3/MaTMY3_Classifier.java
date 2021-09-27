package etmo.metaheuristics.matmy3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;

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
import etmo.metaheuristics.matmy3.models.Classifier;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.math.Probability;
import etmo.util.sorting.NDSortiong;
import smile.classification.DataFrameClassifier;
import smile.classification.LogisticRegression;
import smile.classification.OnlineClassifier;
import smile.classification.RandomForest;
import smile.data.DataFrame;
import smile.data.formula.Formula;

enum EvolutionMode {
    MODEL_ASSIST_SEARCH,
    TRANSFER_EVOLUTION,
    NORMAL_EVOLUTION,
    FINAL_CONVERGE
}

public class MaTMY3_Classifier extends MtoAlgorithm {
    private SolutionSet[] population;
    private SolutionSet[] offspring;
    private SolutionSet[] archive;

    private int populationSize;
    private int archiveSize;
    private int taskNum;
    private int varNum;

    private int generation;
    private int evaluations;
    private int maxEvaluations;

    private String XType;
    private String TXType;

    private Operator DECrossover;
    private Operator SBXCrossover;
    private Operator BLXAlphaCrossover;
    private Operator mutation;

    int[] objStart;
    int[] objEnd;

    double[] bestDistances;
    int[] stuckTimes;
    int[] modelFailTimes;

    boolean isMutate;

    // model
    int epoch = 200;
    double lr = 2e-3;
    Classifier[] models;
    OnlineClassifier<double[]>[] onlineModels;
    Properties randomForestProperties;

    double transferProbability;
    double archiveReplaceRate;
    EvolutionMode[] states; 

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

    public MaTMY3_Classifier(ProblemSet problemSet) {
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initState();
        if (isPlot)
            initPlot();
        while (evaluations < maxEvaluations) {
            iterate();
            if (isPlot)
                updatePlot();
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
        // maxEvaluations = maxEvaluations / taskNum * 2;
        // // problemSet_ = problemSet_.getTask(0);
        // taskNum = 2;

        XType = (String) this.getInputParameter("XType");
        TXType = (String) this.getInputParameter("TXType");
        isPlot = (Boolean) this.getInputParameter("isPlot");
        isMutate = (Boolean) this.getInputParameter("isMutate");
        transferProbability = (Double) this.getInputParameter("transferProbability");

        DECrossover = operators_.get("DECrossover");
        SBXCrossover = operators_.get("SBXCrossover");
        BLXAlphaCrossover = operators_.get("BLXAlphaCrossover");
        mutation = operators_.get("mutation");

        objStart = new int[taskNum];
        objEnd = new int[taskNum];
        bestDistances = new double[taskNum];
        stuckTimes = new int[taskNum];
        modelFailTimes = new int[taskNum];
        Arrays.fill(bestDistances, Double.MAX_VALUE);
        Arrays.fill(stuckTimes, 0);
        Arrays.fill(modelFailTimes, 0);

        models = new Classifier[taskNum];
        randomForestProperties = new Properties();
        randomForestProperties.setProperty("smile.random.forest.max.depth", "16");
        randomForestProperties.setProperty("smile.random.forest.trees", "300");
        randomForestProperties.setProperty("smile.random.forest.mtry", "7");
        randomForestProperties.setProperty("smile.random.forest.max.nodes", String.valueOf(Integer.MAX_VALUE));

        onlineModels = new OnlineClassifier[taskNum];

        archive = new SolutionSet[taskNum];
        archiveSize = 300;

        archiveReplaceRate = 0.5;

        states = new EvolutionMode[taskNum];
        Arrays.fill(states, EvolutionMode.MODEL_ASSIST_SEARCH);

        population = new SolutionSet[taskNum];
        offspring = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
            objStart[k] = problemSet_.get(k).getStartObjPos();
            objEnd[k] = problemSet_.get(k).getEndObjPos();

            population[k] = new SolutionSet(populationSize);
            offspring[k] = new SolutionSet(populationSize);
            archive[k] = new SolutionSet(archiveSize);

            // models[k] = new Classifier(epoch, lr, varNum, varNum);

            for (int i = 0; i < populationSize; i++) {
                Solution solution = new Solution(problemSet_);
                solution.setSkillFactor(k);
                problemSet_.get(k).evaluate(solution);
                evaluations++;
                population[k].add(solution);
            }

            updateBestDistances(k);
        }

        generation ++;

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
    }

    public void iterate() throws JMException, ClassNotFoundException {
        offspringGeneration();
        environmentSelection();
        generation++;
    }

    public void offspringGeneration() throws JMException {
        for (int k = 0; k < taskNum; k++) {
            offspring[k].clear();
            if (generation > 1 && generation < 800 && models[k].getF1() >= 0.57) {
                modelAssistGenerating(k);
                evolutionaryGenerating(k, XType);
                // transferGenerating(k, XType, TXType);
            }
            else if (states[k] == EvolutionMode.TRANSFER_EVOLUTION) {
                transferGenerating(k, XType, TXType);
            }
            else if (states[k] == EvolutionMode.FINAL_CONVERGE) {
                sampleGenerating(k);
            }
            else {
                evolutionaryGenerating(k, XType);
            }
        }
    }

    public void environmentSelection() throws ClassNotFoundException, JMException {
        for (int k = 0; k < taskNum; k++) {
            // System.out.println("G_T: " + generation + "_" + k);

            int transferCount = 0;
            int reproduceCount = 0;
            int notModelCount = 0;
            for (int i = 0; i < offspring[k].size(); i++) {
                if (offspring[k].get(i).getFlag() == 2) {
                    transferCount ++;
                }
                else if (offspring[k].get(i).getFlag() == 3) {
                    notModelCount ++;
                }
                else {
                    reproduceCount ++;
                }
            }

            SolutionSet union = population[k].union(offspring[k]);
            NDSortiong.sort(union, problemSet_, k);
            
            int filtedBetter = 0;
            int normalBetter = 0;
            int notModelBetter = 0;
            for (int i = 0; i < populationSize; i++) {
                if (union.get(i).getFlag() == 2) 
                    filtedBetter ++;
                else if (union.get(i).getFlag() == 1)
                    normalBetter ++;
                else if (union.get(i).getFlag() == 3)
                    notModelBetter ++;
                union.get(i).setFlag(0);
                population[k].replace(i, union.get(i));
            }
            // System.out.println("Normal offspring survival Rate: " + (double)normalBetter / population[k].size());
            // System.out.println("Filted offspring survival Rate: " + (double)filtedBetter / population[k].size());
            // System.out.println(k+" Normal Transfer Success Rate: " + normalBetter + "/" + reproduceCount);
            // System.out.println(k+" Classifier Transfer Success Rate: " + filtedBetter + "/" + transferCount);
            // System.out.println(k+" Not Classifier Transfer Success Rate: " + notModelBetter + "/" + notModelCount);

            if (generation < 800)
                trainModel(k, union.getMat());
        }
    }

    public void updateOnlineModel(int taskID, SolutionSet union) throws JMException {
        int pLength = union.size() / 2;
        int nLength = pLength;
        double[][] features = new double[pLength+nLength][];
        int[] labels = new int[features.length];
        for (int i = 0; i < union.size(); i++) {
            features[i] = union.get(i).getDecisionVariablesInDouble();
            labels[i] = i < pLength ? 0 : 1;
        }

        if (onlineModels[taskID] == null) {
            onlineModels[taskID] = LogisticRegression.fit(features, labels);
        } else {
            onlineModels[taskID].update(features, labels);
        }

        // test
        int[] pre = onlineModels[taskID].predict(features);
        int TP = 0, TN = 0, FP = 0, FN = 0;
        for (int i = 0; i < pre.length; i++) {
            if (pre[i] == 0 && labels[i] == 0) TP ++;
            else if (pre[i] == 0 && labels[i] == 1) FP ++;
            else if (pre[i] == 1 && labels[i] == 0) FN ++;
            else if (pre[i] == 1 && labels[i] == 1) TN ++;
        }
        System.out.println("Classifier " + taskID + ":");
        System.out.println("\t" + TP + "\t" + FN);
        System.out.println("\t" + FP + "\t" + TN);

    }

    public void sampleGenerating(int taskID) throws JMException {
        sampleGenerating(taskID, populationSize);
    }

    public void sampleGenerating(int taskID, int num) throws JMException  {
        double[] mean = population[taskID].getMean();
        double[] std = population[taskID].getStd();
        for (int i = 0; i < num && !offspring[taskID].isFull(); i++) {
            Solution newSolution = new Solution(population[taskID].get(i));
            double[] newFeatures = Probability.sampleByNorm(mean, std);
            newSolution.setDecisionVariables(newFeatures);
            newSolution.setSkillFactor(taskID);
            problemSet_.get(taskID).evaluate(newSolution);
            evaluations ++;
            offspring[taskID].add(newSolution);
        }
    }

    public void evolutionaryGenerating(int taskID, String type) throws JMException {
        evolutionaryGenerating(taskID, taskID, type);
    }

    public void evolutionaryGenerating(int taskID, int assistTaskID, String type) throws JMException {
        if (type.equalsIgnoreCase("SBX")) {
            SBXGenerating(taskID, assistTaskID);
        } else if (type.equalsIgnoreCase("DE")) {
            DEGenerating(taskID, assistTaskID);
        } else if (type.equalsIgnoreCase("BLXAlpha")) {
            BLXGenerating(taskID, assistTaskID);
        }
        else {
            System.out.println("Error: unsupported reproduce type: " + type);
            System.exit(1);
        }
    }

    public Solution generatingChild(int taskID, int assistTaskID, int i, String type) throws JMException {
        Solution child = new Solution();
        if (type.equalsIgnoreCase("SBX")) {
            child = SBXChildGenerating(taskID, assistTaskID, i);
        } else if (type.equalsIgnoreCase("DE")) {
            child = DEChildGenerating(taskID, assistTaskID, i);
        } else if (type.equalsIgnoreCase("BLXAlpha")) {
            child = BLXChildGenerating(taskID, assistTaskID, i);
        }
        else {
            System.out.println("Error: unsupported reproduce type: " + type);
            System.exit(1);
        }

        return child;
    }

    public void modelAssistGenerating(int taskID) throws JMException {
        modelAssistGenerating(taskID, populationSize / 2);
    }

    public void modelAssistGenerating(int taskID, int num) throws JMException {
        for (int i = 0; i < num && !offspring[taskID].isFull(); i++) {
            int[] taskOrder = PseudoRandom.randomPermutation(taskNum, taskNum);
            for (int assistTaskID: taskOrder) {
                if (assistTaskID == taskID) continue;
                if (offspring[taskID].isFull()) break;
                if (judgeIndividual(taskID, population[assistTaskID].get(i))) {
                    // SBX
                    Solution[] parents = new Solution[2];
                    parents[0] = population[taskID].get(i);
                    parents[1] = population[assistTaskID].get(i);
                    Solution newSolution = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];
                        
                    // // DE
                    // int[] pids = PseudoRandom.randomPermutation(populationSize, 2);
                    // int j1 = pids[0], j2 = pids[1];
                    // Solution[] parents = new Solution[3];
                    // parents[0] = population[taskID].get(j1);
                    // parents[1] = population[taskID].get(j2);
                    // parents[2] = population[assistTaskID].get(i);
                    // Solution newSolution = (Solution) DECrossover.execute(new Object[] { population[assistTaskID].get(i), parents });

                    // // explicit
                    // Solution newSolution = new Solution(population[assistTaskID].get(i));

                    newSolution.setSkillFactor(taskID);
                    newSolution.setFlag(2);
                    newSolution.resetObjective();
                    problemSet_.get(taskID).evaluate(newSolution);
                    evaluations ++;
                    offspring[taskID].add(newSolution);
                } 
                // // DEBUG
                // else {
                //     // explicit
                //     Solution newSolution = new Solution(population[assistTaskID].get(i));

                //     newSolution.setSkillFactor(taskID);
                //     newSolution.setFlag(3);
                //     newSolution.resetObjective();
                //     problemSet_.get(taskID).evaluate(newSolution);
                //     evaluations ++;
                //     offspring[taskID].add(newSolution);
                // }
            }
        }
    }

    public void transferGenerating(int taskID, String normalType, String transferType) throws JMException {
        int assistTaskID = findAssistTask(taskID);

        for (int i = 0; i < populationSize && !offspring[taskID].isFull(); i++) {
            if (PseudoRandom.randDouble() < transferProbability) {
                offspring[taskID].add(generatingChild(taskID, assistTaskID, i, transferType));
            } else {
                offspring[taskID].add(generatingChild(taskID, taskID, i, normalType));
            }
        }
    }

    public int findAssistTask(int taskID) {
        int assistTaskID = taskID;
        if (PseudoRandom.randDouble() < transferProbability) {
            while (assistTaskID == taskID) {
                assistTaskID = PseudoRandom.randInt(0, taskNum - 1);
            }
        }
        return assistTaskID;
    }

    public void DEGenerating(int taskID) throws JMException {
        DEGenerating(taskID, taskID);
    }

    public void DEGenerating(int taskID, int assistTaskID) throws JMException {
        for (int i = 0; i < populationSize && !offspring[taskID].isFull(); i++) {
            offspring[taskID].add(DEChildGenerating(taskID, assistTaskID, i));
        }
    }

    public Solution DEChildGenerating(int taskID, int assistTaskID, int i) throws JMException {
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

        mutateIndividual(child);

        child.setSkillFactor(taskID);
        child.setFlag(1);
        problemSet_.get(taskID).evaluate(child);
        evaluations++;

        return child;
    }

    public void SBXGenerating(int taskID) throws JMException {
        SBXGenerating(taskID, taskID);
    }

    public void SBXGenerating(int taskID, int assistTaskID) throws JMException {
        for (int i = 0; i < populationSize && !offspring[taskID].isFull(); i++) {
            offspring[taskID].add(SBXChildGenerating(taskID, assistTaskID, i));
        }
    }

    public Solution SBXChildGenerating(int taskID, int assistTaskID, int i) throws JMException {
        int j = i;
        while (j == i)
            j = PseudoRandom.randInt(0, populationSize - 1);
        Solution[] parents = new Solution[2];
        parents[0] = population[taskID].get(i);
        parents[1] = population[assistTaskID].get(j);

        Solution child = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];

        mutateIndividual(child);

        child.setSkillFactor(taskID);
        child.setFlag(1);
        problemSet_.get(taskID).evaluate(child);
        evaluations++;

        return child;
    }

    public void BLXGenerating(int taskID) throws JMException {
        BLXGenerating(taskID, taskID);
    }

    public void BLXGenerating(int taskID, int assistTaskID) throws JMException {
        for (int i = 0; i < populationSize && !offspring[taskID].isFull(); i++) {
            offspring[taskID].add(BLXChildGenerating(taskID, assistTaskID, i));
        }
    }

    public Solution BLXChildGenerating(int taskID, int assistTaskID, int i) throws JMException {
        int j = i;
        while (j == i)
            j = PseudoRandom.randInt(0, populationSize - 1);
        Solution[] parents = new Solution[2];
        parents[0] = population[taskID].get(i);
        parents[1] = population[assistTaskID].get(j);

        Solution child = ((Solution[]) BLXAlphaCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];

        // mutateIndividual(child);

        child.setSkillFactor(taskID);
        child.setFlag(1);
        problemSet_.get(taskID).evaluate(child);
        evaluations++;
        
        return child;
    }

    public void mutateIndividual(Solution individual) throws JMException {
        // if (isMutate)
        //     mutation.execute(individual);

        if (PseudoRandom.randDouble() < 0.5)
            mutation.execute(individual);
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

    public void trainModel(int taskID, double[][] data) throws JMException {
        int length = data.length;
        double[][] pTrain = new double[length / 2][];
        double[][] nTrain = new double[length / 2][];

        Arrays.parallelSetAll(pTrain, i -> data[i]);
        Arrays.parallelSetAll(nTrain, i -> data[i + length / 2]);

        models[taskID] = new Classifier(epoch, lr, varNum, varNum);
        models[taskID].train(pTrain, nTrain);;
    }

    public boolean judgeIndividual(int taskID, Solution individual) throws JMException {
        double[] features = individual.getDecisionVariablesInDouble();
        return models[taskID].judge(features, 0.5);
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

    public void updateArchive(int taskID, SolutionSet union) {
        for (int i = union.size() / 2; i < union.size(); i++) {
            if (!archive[taskID].isFull()){
                archive[taskID].add(new Solution(union.get(i)));
            }
            else if (PseudoRandom.randDouble() < archiveReplaceRate) {
                archive[taskID].replace(
                    PseudoRandom.randInt(0, archive[taskID].size() - 1), 
                    new Solution(union.get(i)));
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
        double[] varX = new double[varNum];
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
