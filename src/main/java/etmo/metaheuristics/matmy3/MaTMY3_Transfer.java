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
import etmo.metaheuristics.matmy3.models.AbstractDistribution;
import etmo.metaheuristics.matmy3.models.GaussianDistribution;
import etmo.metaheuristics.matmy3.models.MultiVarGaussian;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.math.Distance;
import etmo.util.math.Matrix;
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
    
    private double mutationProbability;
    private double transferProbability;
    private double[][] tP;
    private double[] transferBaseline;
    private double[][] distances;
    private double[][] confidences;
    private int[][] transferredCounts;
    private double[][] lastTransferSuccessRate;

    boolean isMutate;

    double biasPartition;
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
    double[][][] sigmas;
    AbstractDistribution[] models;
    boolean[] isSingular;

    // DEBUG
    int selectVariableID = 1;

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

        mutationProbability = 0.5;
        biasPartition = 0.5;
        center = new double[varNum];
        Arrays.fill(center, 0.5);
        bias = new double[taskNum][varNum];
        momentum = new double[taskNum][varNum];
        alpha = 0.5;
        beta = 0.8;
        
        confidences = new double[taskNum][taskNum];
        transferredCounts = new int[taskNum][taskNum];
        lastTransferSuccessRate = new double[taskNum][taskNum];

        distances = new double[taskNum][taskNum];

        means = new double[taskNum][varNum];
        stds = new double[taskNum][varNum];
        sigmas = new double[taskNum][varNum][varNum];
        models = new AbstractDistribution[taskNum];
        isSingular = new boolean[taskNum];

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
            Arrays.fill(confidences[k], 0.01);
            Arrays.fill(transferredCounts[k], 0);
            Arrays.fill(lastTransferSuccessRate[k], 0.0);

            population[k] = new SolutionSet(populationSize);
            offspring[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++) {
                Solution solution = new Solution(problemSet_);
                solution.setSkillFactor(k);
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
        updateDistributions(0.5);
        // updateDistances();
        // updateBias();
        // if (generation % 50 == 1)

        for (int k = 0; k < taskNum; k++) {
            offspring[k].clear();
            Solution child;

            Arrays.fill(models, null);
            for (int kk = 0; kk < taskNum; kk++) {
                // if (kk == k) continue;
                // try {
                //     models[kk] = new MultiVarGaussian(means[k], sigmas[kk]);
                // }
                // catch (org.apache.commons.math3.linear.SingularMatrixException e) {
                //     models[kk] = new GaussianDistribution(means[k], stds[kk]);
                // }
                models[kk] = new GaussianDistribution(means[k], stds[kk]);
            }

            int[] perm = PseudoRandom.randomPermutation(population[k].size(), population[k].size());
            for (int i = 0; i < populationSize; i ++) {
                child = null;
                if (generation <= 1000) {
                    if (PseudoRandom.randDouble() < transferProbability) {
                        int k2 = getAssistTaskID(k);
                        child = transferGenerating(k, k2, perm[i], TXType);
                        transferredCounts[k][k2] ++;
                    } else {
                        child = evolutionaryGenerating(k, perm[i], XType);
                    }
                } else {
                    child = sampleGenerating(k, perm[i]);
                }
                evaluate(child, k);
                offspring[k].add(child);
            }
        }
    }

    Solution sampleGenerating(int taskID, int i) throws JMException {
        Solution child = null;
        child = new Solution(population[taskID].get(i));
        double[] newFeatures = Probability.sampleByNorm(means[taskID], stds[taskID]);
        Vector.vecClip_(newFeatures, 0.0, 1.0);
        child.setDecisionVariables(newFeatures);
        child.resetObjective();
        return child;
    }

    Solution transferGenerating(int taskID, int assistTaskID, int i, String type) throws JMException {
        int j = PseudoRandom.randInt(0, population[assistTaskID].size() - 1);
        
        // explicit
        Solution child = null;
        if (PseudoRandom.randDouble() < 1.0) {
            child = new Solution(population[assistTaskID].get(j));
            // double[] newFeatures = Probability.sampleByNorm(means[taskID], sigmas[assistTaskID]);
            // double[] newFeatures = Probability.sampleByNorm(means[taskID], stds[assistTaskID]);
            // double[] newFeatures = new MultivariateNormalDistribution(
            //     means[taskID], sigmas[assistTaskID]).sample();
            double[] newFeatures = null;
            newFeatures = models[assistTaskID].sample();
            Vector.vecClip_(newFeatures, 0.0, 1.0);
            child.setDecisionVariables(newFeatures);
            child.resetObjective();
        } else {
            child = offsetIndividual(population[assistTaskID].get(j), assistTaskID, taskID);
        }
        // SBX implicit
        // Solution[] parents = new Solution[2];
        // parents[0] = population[taskID].get(i);
        // parents[1] = offsetIndividual(population[assistTaskID].get(j), assistTaskID, taskID);
        // Solution child = ((Solution[]) SBXCrossover.execute(parents))[PseudoRandom.randInt(0, 1)];
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

            int normalBetter = 0;
            int[] transferBetter = new int[taskNum];
            for (int i = 0; i < populationSize; i++) {
                if (union.get(i).getFlag() == 1)
                    normalBetter ++;
                else if (union.get(i).getFlag() == 2)
                    transferBetter[union.get(i).getSkillFactor()] ++;
                
                // union.get(i).setFlag(0);
                union.get(i).setSkillFactor(k);
                population[k].replace(i, union.get(i));
            }
            
            // System.out.println("Task " + k + ": " + Arrays.toString(bias[k]));
            // System.out.println("G_T: " + generation + "_" + k);
            // System.out.println("Mutation offspring survival Rate: " + (double)normalBetter / population[k].size());
            // System.out.println("Normal offspring survival Rate: " + (double)normalBetter / population[k].size());
            // System.out.println("Transfer offspring survival Rate: " + (double)transferBetter / population[k].size());
            // System.out.println("Transfer Success Rate: " + (double)transferBetter / transferTotalCount[k]);
            
            // transferTotalCount[k] = 0;

            // double rate = 0;
            // for (int kk = 0; kk < taskNum; kk++) {
            //     if (transferredCounts[k][kk] == 0) continue;
            //     rate = transferBetter[kk] * 1.0 / transferredCounts[k][kk]; 
            //     if (rate < lastTransferSuccessRate[k][kk]) {
            //         confidences[k][kk] = Math.max(0, confidences[k][kk] - 0.1);
            //     }
            //     lastTransferSuccessRate[k][kk] = rate;
            // }
            Arrays.fill(transferredCounts[k], 0);

            // if (k == plotTaskID) 
            //     // System.out.println(Arrays.toString(lastTransferSuccessRate[k]));
            //     System.out.println(Arrays.toString(confidences[k]));
            //     // System.out.println(Arrays.toString(lastTransferSuccessRate[k]));

            updateBestDistances(k);

            // // 灾变
            // catastrophe(k, 0, 8);
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

    void updateBias() throws JMException {
        for (int k = 0; k < taskNum; k++) {
            updateBias(k, biasPartition);
        }
    }

    void updateBias(int taskID, double partition) throws JMException {
        // update bias distance
        double[] oldBias = bias[taskID].clone();
        Arrays.fill(bias[taskID], 0);
        for (int i = 0; i < (int)(populationSize * partition); i ++) {
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

    void updateDistributions(double partition) throws JMException {
        double[] weights = null;
        SolutionSet[] tmpSet = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
            int size = (int)(population[k].size() * partition);
            weights = new double[size];
            tmpSet[k] = new SolutionSet(size);
            for (int i = 0; i < size; i++) {
                tmpSet[k].add(population[k].get(i));
                weights[i] = 1 / (population[k].get(i).getRank() + 1.0);
            }
            // // means[k] = tmpSet.getWeightedMean(weights);
            // means[k] = tmpSet[k].getMean();
            // stds[k] = tmpSet[k].getStd();
            // sigmas[k] = Matrix.getMatSigma(tmpSet[k].getMat());
        }
        Arrays.parallelSetAll(means, k -> tmpSet[k].getMean());
        Arrays.parallelSetAll(stds, k -> tmpSet[k].getStd());
        Arrays.parallelSetAll(sigmas, k -> {
            double[][] output = null;
            try {
				output = Matrix.getMatSigma(tmpSet[k].getMat());
			} catch (JMException e) {
				e.printStackTrace();
			}
            return output;
		});
    }

    int getAssistTaskID(int taskID) throws JMException {
        int assistTaskID = taskID;

        while (assistTaskID == taskID) {
            assistTaskID = PseudoRandom.randInt(0, taskNum - 1);
        }

        // double[] scores = new double[taskNum];
        // for (int i = 0; i < scores.length; i ++) {
        //     scores[i] = 1.0 / distances[taskID][i];
        // }
        // assistTaskID = Random.rouletteWheel(scores, taskID);

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
        
        if (Math.abs(avgDistance - bestDistances[taskID]) > bestDistances[taskID] * 5e-3) {
            if (avgDistance < bestDistances[taskID]) {
                updated = true;
            }
            bestDistances[taskID] = avgDistance;
        }
        
        // DEBUG
        // if (taskID == plotTaskID) {
        //     if (updated)
        //     System.out.println(avgDistance);
        //     else
        //     System.out.println("stucking " + stuckTimes[plotTaskID] + " ...");
        // }
            
        if (updated) {
            stuckTimes[taskID] = 0;
        } else {

            stuckTimes[taskID]++;
        }
    }
    
    void updateDistances() throws JMException {
        // // t-SNE
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
        //     e.printStackTrace();
        // }
        // double[][] points = tSNE.coordinates;
        // double[][][] clusters = new double[taskNum][][];
        // for (int k = 0; k < taskNum; k++) {
        //     clusters[k] = Matrix.matSlice(points, 
        //     k * (populationSize / 2), (k + 1) * (populationSize / 2) - 1, 
        //     0, points[0].length - 1);
        // }

        for (int k = 0; k < taskNum; k++) {
            for (int kk = k + 1; kk < taskNum; kk++) {
                // WD
                // distances[k][kk] = distances[kk][k] = Distance.getWassersteinDistance(
                //     population[kk].getMat(),
                //     population[k].getMat());

                // Cosine
                // distances[k][kk] = distances[kk][k] = Distance.getCosineSimilarity(bias[k], bias[kk]);
            
                // // t-SNE + WD
                // distances[k][kk] = distances[kk][k] = 
                //     Distance.getWassersteinDistance(clusters[k], clusters[kk]);

                // co-variance matrix similarity
                // Coral Loss
                double dist = 0;
                for (int i = 0; i < varNum; i ++) {
                    for (int j = 0; j < varNum; j ++) {
                        dist += Math.pow(sigmas[k][i][j] - sigmas[kk][i][j], 2);
                    }
                }
                dist = Math.sqrt(dist) / (4.0 * Math.pow(varNum, 2));
            }
        }
    }
        
    void evaluate(Solution solution, int taskID) throws JMException {
        // solution.setSkillFactor(taskID);
        problemSet_.get(taskID).evaluate(solution);
        evaluations ++;
    }

    void catastrophe(int taskID, double survivalRate, int threshold) throws ClassNotFoundException, JMException {
        if (stuckTimes[taskID] >= threshold
                && evaluations < maxEvaluations - 100) {
            // System.out.println(evaluations + ": task " + taskID +" : reset on " + selectVariableID);
            int[] perm = PseudoRandom.randomPermutation(populationSize, (int) (populationSize * (1 - survivalRate)));
            double[] ms = population[taskID].getMean();
            double[] ss = population[taskID].getStd();
            for (int i = 0; i < perm.length; i++) {
                // // totally random individual
                // Solution solution = new Solution(problemSet_);
                // solution.setSkillFactor(taskID);
                // problemSet_.get(taskID).evaluate(solution);
                // // evaluations++;
                // population[taskID].replace(perm[i], solution);
                
                // partily random vaiable
                Vector.vecElemMul_(ss, selectVariableID);
                Vector.vecClip_(ss, 0.0, 0.5);
                Solution tmp = population[taskID].get(perm[i]);
                double[] newFeatures = Probability.sampleByNorm(ms, ss);
                Vector.vecClip_(newFeatures, 0.0, 1.0);
                tmp.setDecisionVariables(newFeatures);
                problemSet_.get(taskID).evaluate(tmp);
                population[taskID].replace(perm[i], tmp);
            }
            if (taskID == plotTaskID){
                // XType = "DE";
                selectVariableID *= 2;
                System.out.println(generation + ": " + selectVariableID);
            }
            stuckTimes[taskID] = 0;
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

    void resetFlag() {
        for (int k = 0; k < taskNum; k++) {
            for (int i = 0; i < population[k].size(); i++) {
                population[k].get(i).setFlag(0);
            }
        }
    }
}
