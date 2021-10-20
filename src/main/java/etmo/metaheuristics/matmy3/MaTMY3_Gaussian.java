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

import etmo.core.MtoAlgorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.metaheuristics.matmy3.models.AbstractDistribution;
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

public class MaTMY3_Gaussian extends MtoAlgorithm {
    private SolutionSet[] population;
    private Solution[][] offspring;
    private SolutionSet[] union;

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
    private double[] previousBestDistances;
    
    private double mutationProbability;
    private double transferProbability;
    private double[] tP;
    private double[][] distances;
    private double[][] distances2;
    private double[][] confidences;
    private int[][] transferredCounts;
    private int[][] transferredEliteCounts;
    private double[][] lastTransferSuccessRate;

    private double elitePart;

    boolean isMutate;

    double[][] means;
    double[][] stds;
    double[][][] sigmas;
    AbstractDistribution[] models;
    boolean[] isSingular;
    double[][] eliteDirections;

    // parallel runner
    Object[] runner;
    final Object lock = new Object();

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

    public MaTMY3_Gaussian(ProblemSet problemSet) {
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
        // System.out.println(igdPlotValues.get(plotTaskID).toString());

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
        mutationProbability = (Double) this.getInputParameter("mutationProbability");
        
        elitePart = (Double) this.getInputParameter("elitePartition");

        tP = new double[taskNum];
        Arrays.fill(tP, transferProbability);

        // DEBUG: IGD PLOTTING
        plotTaskID = (Integer) this.getInputParameter("plotTaskID");

        DECrossover = operators_.get("DECrossover");
        SBXCrossover = operators_.get("SBXCrossover");
        mutation = operators_.get("mutation");

        objStart = new int[taskNum];
        objEnd = new int[taskNum];
        bestDistances = new double[taskNum];
        previousBestDistances = new double[taskNum];
        stuckTimes = new int[taskNum];
        Arrays.fill(bestDistances, Double.MAX_VALUE);
        Arrays.fill(previousBestDistances, Double.MAX_VALUE);
        Arrays.fill(stuckTimes, 0);
        
        confidences = new double[taskNum][taskNum];
        transferredCounts = new int[taskNum][taskNum];
        lastTransferSuccessRate = new double[taskNum][taskNum];
        transferredEliteCounts = new int[taskNum][taskNum];

        distances = new double[taskNum][taskNum];
        distances2 = new double[taskNum][taskNum];

        means = new double[taskNum][varNum];
        stds = new double[taskNum][varNum];
        sigmas = new double[taskNum][varNum][varNum];
        models = new AbstractDistribution[taskNum];
        isSingular = new boolean[taskNum];
        eliteDirections = new double[taskNum][varNum];


        population = new SolutionSet[taskNum];
        union = new SolutionSet[taskNum];
        // offspring = new SolutionSet[taskNum];
        offspring = new Solution[taskNum][populationSize];
        for (int k = 0; k < taskNum; k++) {
            objStart[k] = problemSet_.get(k).getStartObjPos();
            objEnd[k] = problemSet_.get(k).getEndObjPos();

            Arrays.fill(distances[k], 0);
            Arrays.fill(distances2[k], 0);
            Arrays.fill(means[k], 0);
            Arrays.fill(stds[k], 0);
            Arrays.fill(eliteDirections[k], 0);

            Arrays.fill(confidences[k], 0.01);
            Arrays.fill(transferredCounts[k], 0);
            Arrays.fill(lastTransferSuccessRate[k], 0.0);

            population[k] = new SolutionSet(populationSize);
            union[k] = new SolutionSet();
            for (int i = 0; i < populationSize; i++) {
                // if (k == 0) {
                Solution solution = new Solution(problemSet_);
                solution.setSkillFactor(k);
                evaluate(solution, k);
                population[k].add(solution);
                // } else {
                //     Solution solution = new Solution(population[0].get(i));
                //     solution.setSkillFactor(k);
                //     evaluate(solution, k);
                //     population[k].add(solution);
                // }
            }
            NDSortiong.sort(population[k], problemSet_, k);
            for (int i = 0; i < populationSize; i ++) {
                population[k].get(i).setFlag2(population[k].get(i).getRank());
            }
            // updateBestDistances(k);
        }

        runner = new Object[taskNum];

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
        // long time = System.currentTimeMillis();
        offspringGeneration();
        // System.out.println("offspringGeneration: " + (System.currentTimeMillis() - time));
        // time = System.currentTimeMillis();
        environmentSelection();
        // System.out.println("environmentSelection: " + (System.currentTimeMillis() - time));
        // writePopulationVariablesMatrix(plotTaskID, generation);
        generation++;
    }

    void offspringGeneration() throws JMException {
        updateDistributions(elitePart);
        updateDistances();

        // parallel
        Arrays.parallelSetAll(runner, k -> {
            Object res = null;
            try {
                res = generatingOffspring(k);
            } catch (JMException e) {
                e.printStackTrace();
            }
            return res;
        });

        // for (int k = 0; k < taskNum; k++) {
        //     Arrays.fill(offspring[k], null);
        //     Solution child;

        //     int[] perm = PseudoRandom.randomPermutation(population[k].size(), population[k].size());
        //     for (int i = 0; i < populationSize; i ++) {
        //         child = null;
        //         // if (PseudoRandom.randDouble() < transferProbability) {
        //         //     int k2 = getAssistTaskID(k);
        //         //     child = transferGenerating(k, k2, perm[i], TXType);
        //         // } else {
        //         //     child = evolutionaryGenerating(k, perm[i], XType);
        //         // }

        //         evaluate(child, k);
        //         offspring[k][i] = child;
        //     }
        // }
    }

    Object generatingOffspring(int taskID) throws JMException {
        Arrays.fill(offspring[taskID], null);
        int[] perm = PseudoRandom.randomPermutation(population[taskID].size(), population[taskID].size());
        
        Solution child;
        for (int i = 0; i < populationSize; i ++) {
            child = null;
            if (PseudoRandom.randDouble() < tP[taskID]) {
                int k2 = getAssistTaskID(taskID);
            // if (i < populationSize / 2) {
            //     int k2 = i;
                child = transferGenerating(taskID, k2, perm[i], TXType);
            //  int[] assistTaskList = getAssistTaskIDs(taskID, 10);
            //  double[] 
            } else {
                child = evolutionaryGenerating(taskID, perm[i], XType);
            }
            evaluate(child, taskID);
            child.setFlag2(-1);
            offspring[taskID][i] = child;
        }
        return null;
    }

    Solution transferGenerating(int taskID, int assistTaskID, int i, String type) throws JMException {
        int j = PseudoRandom.randInt(0, population[assistTaskID].size() - 1);
        
        // explicit
        Solution child = null;
        child = new Solution(population[assistTaskID].get(j));

        double[] tmpMean = population[assistTaskID].getMean();
        Vector.vecSub_(tmpMean, means[assistTaskID]);
        Vector.vecAdd_(tmpMean, means[taskID]);

        // double[] tmpMean = means[taskID];
        // Vector.vecSub_(tmpMean, population[assistTaskID].getMean());
        // Vector.vecAdd_(tmpMean, means[assistTaskID]);

        // double[] tmpStd = stds[assistTaskID];
        double[] tmpStd = population[assistTaskID].getStd();
        // double[] tmpStd = population[taskID].getStd();

        // double[] newFeatures = Probability.sampleByNorm(tmpMean, tmpStd);
        double[][] tmpSigma = Matrix.getMatSigma(population[assistTaskID].getMat());
        double[] newFeatures = Probability.sampleByNorm(tmpMean, tmpSigma);

        Vector.vecClip_(newFeatures, 0.0, 1.0);
        child.setDecisionVariables(newFeatures);
        mutateIndividual(taskID, child);
        child.resetObjective();

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
        // TODO: failed to parallelize
        // Arrays.parallelSetAll(runner,  i -> environmentSelection(i));
        
        for (int taskID = 0; taskID < taskNum; taskID ++) {
            SolutionSet offspringSet = new SolutionSet(offspring[taskID]);
            SolutionSet union = population[taskID].union(offspringSet);
            NDSortiong.sort(union, problemSet_, taskID);
    
            // double improvement = 0;
            // for (int i = 0; i < union.size(); i ++) {
            //     if (union.get(i).getFlag2() >= 0) {
            //         improvement += (union.get(i).getRank() - union.get(i).getFlag2());
            //     }
            // }
            // improvement /= populationSize;
            // tP[taskID] = improvement > 1 ? 0.5 : 0.5 * (2 - improvement);
            // tP[taskID] = improvement > 1 ? 0.5 : 0.5 * improvement;
            // if (taskID == plotTaskID) {
            //     System.out.println(generation + ": " + improvement);
            //     System.out.println(generation + ": " + tP[taskID]);
            // }
            // Arrays.fill(transferredEliteCounts[taskID], 0);
            // int eliteCount = 0;
            // int transferredEliteCount = 0;
            for (int i = 0; i < populationSize; i++) {
                // if (union.get(i).getRank() <= 1) {
                // if (i < populationSize / 2) {
                //     eliteCount ++;
                //     int sf = union.get(i).getSkillFactor();
                //     if (sf != taskID) {
                //         transferredEliteCounts[taskID][sf] = 1;
                //         transferredEliteCount ++;
                //     }
                // }
                union.get(i).setFlag2(union.get(i).getRank());
                union.get(i).setSkillFactor(taskID);
                population[taskID].replace(i, union.get(i));
            }

            // tP[taskID] = 0.9 * tP[taskID] + 0.1 * (transferredEliteCount * 1.0 / eliteCount / tP[taskID]);
            // // if (taskID == plotTaskID) {
            // //     System.out.println(tP[taskID]);
            // //     System.out.println(generation + ": " + transferredEliteCount + " / " + eliteCount);
            // //     System.out.println(Arrays.toString(transferredEliteCounts[taskID]));
            // // }
            // tP[taskID] = Math.max(0.1, tP[taskID]);
            // tP[taskID] = Math.min(0.9, tP[taskID]);

            // updateBestDistances(taskID);
        }
    }

    Object environmentSelection(int taskID) {
        // SolutionSet offspringSet = new SolutionSet(offspring[taskID]);
        // SolutionSet union = population[taskID].union(offspringSet);
        union[taskID] = population[taskID].union(offspring[taskID]);
        NDSortiong.sort(union[taskID], problemSet_, taskID);

        synchronized(lock) {
            for (int i = 0; i < populationSize; i++) {
                union[taskID].get(i).setSkillFactor(taskID);
                population[taskID].replace(i, union[taskID].get(i));
            }
        }

        return null;
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

    void updateDistributions(double partition) throws JMException {
        // double[] weights = null;
        SolutionSet[] tmpSet = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {

             int size = (int)(population[k].size() * partition);
            //  weights = new double[size];
             tmpSet[k] = new SolutionSet(size);
             for (int i = 0; i < size; i++) {
                 tmpSet[k].add(population[k].get(i));
                //  weights[i] = 1 / (population[k].get(i).getRank() + 1.0);
             }

        //    tmpSet[k] = new SolutionSet();
        //    for (int i = 0; i < population[k].size(); i ++) {
        //        if ((population[k].get(i).getRank() == 0 || tmpSet[k].size() < 10) && i < 50) {
        //            tmpSet[k].add(population[k].get(i));
        //        }
        //        else {
        //            break;
        //        }
        //    }

            // means[k] = tmpSet[k].getWeightedMean(weights);
            // stds[k] = tmpSet[k].getWeightedStd(weights);
            // means[k] = tmpSet[k].getMean();
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
        // Arrays.setAll(eliteDirections, k -> Vector.vecSub(means[k], population[k].getMean()));
    }

    int[] getAssistTaskIDs(int taskID, int num) {
        int[] idx = new int[num];
        Arrays.fill(idx, -1);
        for (int k = 0; k < taskNum; k++) {
            if (k == taskID) continue;
            for (int i = 0; i < num; i ++) {
                if (idx[i] == -1) {
                    idx[i] = k;
                }
                else if (distances[taskID][idx[i]] > distances[taskID][k]) {
                    for (int j = i + 1; j < num; j++) {
                        if (idx[j] == -1) {
                            idx[j] = idx[j-1];
                            break;
                        }
                        idx[j] = idx[j-1];
                    }
                    idx[i] = k;
                }
            }
        }
        return idx;
    }

    int getAssistTaskID(int taskID) throws JMException {
        int assistTaskID = taskID;

    //    // random
    //    while (assistTaskID == taskID) {
    //        assistTaskID = PseudoRandom.randInt(0, taskNum - 1);
    //    }

        // double[] scores = new double[taskNum];
        // Arrays.setAll(scores, i -> transferredEliteCounts[taskID][i]);
        // assistTaskID = Random.rouletteWheel(scores, taskID);

        // double[] scores = new double[taskNum];
        // Arrays.setAll(scores, i -> distances[taskID][i]);
        // assistTaskID = Random.rouletteWheel(scores, taskID);

        double[] scores = new double[taskNum];
        // CMD
        Arrays.setAll(scores, i -> distances[taskID][i]);
        int res1 = Random.rouletteWheel(scores, taskID);
        // elite mean distance
        Arrays.setAll(scores, i -> 1 / distances2[taskID][i]);
        int res2 = Random.rouletteWheel(scores, taskID);
        assistTaskID = PseudoRandom.randDouble() < 0.5 ? res1 : res2;

        // int length = 10;
        // double[] scores = new double[length];
        // while (assistTaskID == taskID) {
        //     int[] perm = PseudoRandom.randomPermutation(taskNum, length);
        //     Arrays.setAll(scores, i -> distances[taskID][perm[i]]);
        //     assistTaskID = perm[Random.rouletteWheel(scores)];
        // }

        // // coral distance
        // int[] perm = PseudoRandom.randomPermutation(taskNum, taskNum);
        // double minDist = Double.MAX_VALUE;
        // int selectTaskID = -1;
        // // TODO: hard code
        // for (int k = 0; k < 10; k ++) {
        //     if (distances[taskID][k] < minDist) {
        //         minDist = distances[taskID][perm[k]];
        //         selectTaskID = perm[k];
        //     }
        // }
        // assistTaskID = selectTaskID == -1 ? taskID : selectTaskID;
        
        // assistTaskID = Random.rouletteWheel(distances[taskID], taskID);

        return assistTaskID;
    }
  
    void updateDistances() throws JMException {
        for (int k = 0; k < taskNum; k++) {
           final int srcTaskID = k;
           Arrays.parallelSetAll(distances[k], trgTaskID -> {
                double dist = 0;
                if (trgTaskID > srcTaskID) {
                        // int d = sigmas[srcTaskID].length;
                        // for (int i = 0; i < d; i ++) {
                        //     for (int j = 0; j < d; j ++) {
                        //         dist += Math.pow(sigmas[srcTaskID][i][j] - sigmas[trgTaskID][i][j], 2);
                        //     }
                        // }
                        // dist = Math.sqrt(dist);
                        dist = Distance.getCorrelationMatrixDistance(sigmas[srcTaskID], sigmas[trgTaskID]);    
                        // dist = Distance.getCosineSimilarity(eliteDirections[srcTaskID], eliteDirections[trgTaskID]);
                        // dist = (1 + dist) / 2 * -1;
                    } else {
                        dist = distances[trgTaskID][srcTaskID];
                    }
                return dist;
            });
            Arrays.parallelSetAll(distances2[k], trgTaskID -> {
            double dist = 0;
            if (trgTaskID > srcTaskID) {
                    // int d = sigmas[srcTaskID].length;
                    // for (int i = 0; i < d; i ++) {
                    //     for (int j = 0; j < d; j ++) {
                    //         dist += Math.pow(sigmas[srcTaskID][i][j] - sigmas[trgTaskID][i][j], 2);
                    //     }
                    // }
                    // dist = Math.sqrt(dist);
                    dist = Distance.getDistance(means[srcTaskID], means[trgTaskID]);    
                } else {
                    dist = distances2[trgTaskID][srcTaskID];
                }
            return dist;
            });
        }
    }
      
    void updateBestDistances(int taskID) {
        // boolean updated = false;
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

        if (generation >= 2) {
            double improve1 = Math.max(0, previousBestDistances[taskID] - bestDistances[taskID]);
            double improve2 = Math.max(0, bestDistances[taskID] - avgDistance);
            tP[taskID] = improve2 / (improve1 + improve2 + 1e-13);
            // DEBUG
            if (taskID == plotTaskID)
                System.out.println(taskID + ": " + tP[taskID]);
        }

        previousBestDistances[taskID] = bestDistances[taskID];
        bestDistances[taskID] = avgDistance;

        
        // if (Math.abs(avgDistance - bestDistances[taskID]) > bestDistances[taskID] * 5e-3) {
        //     if (avgDistance < bestDistances[taskID]) {
        //         updated = true;
        //     }
        //     bestDistances[taskID] = avgDistance;
        // }
            
        // if (updated) {
        //     stuckTimes[taskID] = 0;
        // } else {
        //     stuckTimes[taskID]++;
        // }
    }

    void evaluate(Solution solution, int taskID) throws JMException {
        synchronized (lock){
            // solution.setSkillFactor(taskID);
            problemSet_.get(taskID).evaluate(solution);
            evaluations ++;
        }
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
