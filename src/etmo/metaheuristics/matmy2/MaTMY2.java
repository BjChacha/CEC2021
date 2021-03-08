package etmo.metaheuristics.matmy2;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.util.JMException;
import etmo.metaheuristics.matmy2.libs.*;
import etmo.util.PseudoRandom;
import etmo.util.logging.LogPopulation;
import jdk.jshell.execution.Util;

import javax.management.openmbean.OpenMBeanParameterInfo;
import java.util.Arrays;
import java.util.HashMap;

public class MaTMY2 extends MtoAlgorithm {
    MaTAlgorithm[] optimizers;
    SolutionSet[] populations;

    Operator crossover;
    Operator selection;

    private int populationSize;
    private int maxEvaluations;
    private int evaluations;
    int taskNum;

    int[] transferSourceIndexes;
    int[] runTimes;
    int transferVolume;
    int baseRunTime;
    int minRunTime;

    boolean[] inProgress;
    int inProgressTaskNum;

    double[][] scores;
    double scoreIncrement;
    double scoreDecreaseRate;

    double[][] convergeDirections;
    double[][] populationCenterPositions;

    double forceTransferRate;

    boolean isDRA;

    private String algoName;

    // TODO: DEBUG
    int[] debugRunTimes;
    int goodTransferCnt = 0;
    int badTransferCnt = 0;

    public MaTMY2(ProblemSet problemSet) {
        super(problemSet);
    }


    private void initState(){
        populationSize = (Integer) getInputParameter("populationSize");
        maxEvaluations = (Integer) getInputParameter("maxEvaluations");
        transferVolume = (Integer) getInputParameter("transferVolume");
        baseRunTime = (Integer) getInputParameter("baseRunTime");

        scoreIncrement = (Double) getInputParameter("scoreIncrement");
        scoreDecreaseRate = (Double) getInputParameter("scoreDecreaseRate");

        forceTransferRate = (Double) getInputParameter("forceTransferRate");

        isDRA = (Boolean) getInputParameter("isDRA");

        algoName = (String) getInputParameter("algoName");

        crossover = operators_.get("crossover");
        selection = operators_.get("selection");

        taskNum = problemSet_.size();
        transferSourceIndexes = new int[taskNum];
        Arrays.fill(transferSourceIndexes, -1);

        inProgress = new boolean[taskNum];
        Arrays.fill(inProgress, true);

        inProgressTaskNum = taskNum;
        runTimes = new int[taskNum];
        minRunTime = 1;

        scores = new double[taskNum][taskNum];
        for (int k = 0; k < taskNum; k++){
            Arrays.fill(scores[k], 1);
            scores[k][k] = 0;
        }

        convergeDirections = new double[taskNum][];
        for (int k = 0; k < taskNum; k++){
            convergeDirections[k] = new double[problemSet_.get(k).getNumberOfVariables()];
        }

        populationCenterPositions = new double[taskNum][];
        for (int k = 0; k < taskNum; k++){
            populationCenterPositions[k] = new double[problemSet_.get(k).getNumberOfVariables()];
        }

        // TODO: DEBUG
        debugRunTimes = new int[taskNum];
        Arrays.fill(debugRunTimes, 0);
    }

    private void initPopulations() throws ClassNotFoundException, JMException {
        populations = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++){
            populations[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++){
                Solution newSolution = new Solution(problemSet_.getTask(k));
                problemSet_.get(k).evaluate(newSolution);
                populations[k].add(newSolution);
            }
        }
        evaluations += (taskNum * populationSize);
        Arrays.fill(debugRunTimes, 1);
        LogPopulation.LogPopulation("MaTMY2", populations, problemSet_, evaluations);
    }

    private void initOptimizers() throws JMException, ClassNotFoundException {
        optimizers = new MaTAlgorithm[taskNum];
        initPopulations();

        if (algoName.equalsIgnoreCase("MOEAD")){
            Operator crossover;
            Operator mutation;

            HashMap parameters;

            for (int k = 0; k < taskNum; k++) {
                ProblemSet pS = problemSet_.getTask(k);
                optimizers[k] = new MOEAD(pS, populations[k]);
                optimizers[k].setInputParameter("populationSize", 100);
                optimizers[k].setInputParameter("maxEvaluations", 1000 * 100);

                optimizers[k].setInputParameter("dataDirectory", "D:\\_r\\EA\\ETMO\\MTO-cec2021-\\resources\\weightVectorFiles\\moead");

                optimizers[k].setInputParameter("T", 20);
                optimizers[k].setInputParameter("delta", 0.9);
                optimizers[k].setInputParameter("nr", 2);

                parameters = new HashMap();
                parameters.put("CR", 1.0);
                parameters.put("F", 0.5);
                crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / pS.get(0).getNumberOfVariables());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                optimizers[k].addOperator("crossover", crossover);
                optimizers[k].addOperator("mutation", mutation);

                optimizers[k].initState();
            }
        }
        else if (algoName.equalsIgnoreCase("MaOEAC")){
            Operator crossover;
            Operator mutation;
            Operator selection;

            HashMap parameters;

            for (int k = 0; k < taskNum; k++) {
                ProblemSet pS = problemSet_.getTask(k);
                optimizers[k] = new MaOEAC(pS, populations[k]);

                optimizers[k].setInputParameter("populationSize",100);
                optimizers[k].setInputParameter("maxGenerations",1000);

                parameters = new HashMap();
                parameters.put("probability", 1.0);
                parameters.put("distributionIndex", 30.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / pS.get(0).getNumberOfVariables());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                parameters = null;
                selection = SelectionFactory.getSelectionOperator("RandomSelection", parameters);

                optimizers[k].addOperator("crossover", crossover);
                optimizers[k].addOperator("mutation", mutation);
                optimizers[k].addOperator("selection", selection);

                optimizers[k].initState();
            }
        }
        else if (algoName.equalsIgnoreCase("NSGAII")){
            Operator crossover;
            Operator mutation;
            Operator selection;

            HashMap parameters;

            for (int k = 0; k < taskNum; k++){
                ProblemSet pS = problemSet_.getTask(k);
                optimizers[k] = new NSGAII(pS, populations[k]);

                optimizers[k].setInputParameter("populationSize", 100);
                optimizers[k].setInputParameter("maxEvaluations", 100 * 1000);

                parameters = new HashMap();
                parameters.put("probability", 0.9);
                parameters.put("distributionIndex", 20.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / pS.getMaxDimension());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                parameters = null;
                selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters);

                optimizers[k].addOperator("crossover", crossover);
                optimizers[k].addOperator("mutation", mutation);
                optimizers[k].addOperator("selection", selection);

                optimizers[k].initState();
            }
        }
        else{
            throw new JMException("Exception in " + algoName + ": Not such algorithm.");
        }
    }

    private void Converge(int times) throws JMException, ClassNotFoundException {
        int[] tmpTimes = new int[taskNum];
        Arrays.fill(tmpTimes, times);
        Converge(tmpTimes);
    }

    private void Converge(int[] times) throws JMException, ClassNotFoundException {
        for (int k = 0; k < taskNum; k++) {
            for (int t = 0; t < times[k]; t++) {
                if (inProgress[k]) {
                    inProgress[k] = optimizers[k].step();
                    evaluations += populationSize;
                    if (!inProgress[k])
                        inProgressTaskNum--;

//                    // TODO: DEBUG
//                    debugRunTimes[k] ++;
//                    System.out.println(Arrays.toString(debugRunTimes));
//                    System.out.println(evaluations);
//                    System.out.println(IntStream.of(debugRunTimes).sum() * populationSize);
                }

                if (evaluations % (20 * taskNum * populationSize) == 0) {
//                System.out.println(totalCount+"x"+taskNum+"x"+populationSize+"="+totalCount * taskNum * populationSize);
                    LogPopulation.LogPopulation("MaTMY2", populations, problemSet_, evaluations);
                }
            }
        }
    }

    private void calculativeConverge() throws JMException, ClassNotFoundException {
        for (int k = 0; k < taskNum; k++){
            populationCenterPositions[k] = populations[k].getAverageDecisionVector();
        }
        Converge(3);
        for (int k = 0; k < taskNum; k++){
            double[] newPosition = populations[k].getAverageDecisionVector();
            convergeDirections[k] = Utils.vectorMinus(newPosition, populationCenterPositions[k]);
            populationCenterPositions[k] = newPosition;
        }
    }

    private void tentativeTransfer() throws JMException, ClassNotFoundException {
        transferSourceIndexes = getTransferSourceIndexes();
        transferIndividuals();

        double[][] firstBestVectors = new double[taskNum][];
        for (int k = 0; k < taskNum; k++){
            firstBestVectors[k] = populations[k].getBestObjectiveVector();
//            firstBestVectors[k] = populations[k].getAverageObjectiveVector();
        }

        Converge(1);

        double[][] secondBestVectors = new double[taskNum][];
        for (int k = 0; k < taskNum; k++){
            secondBestVectors[k] = populations[k].getBestObjectiveVector();
//            secondBestVectors[k] = populations[k].getAverageObjectiveVector();
        }

        processConvergeImprovement(firstBestVectors, secondBestVectors);
    }

    private void selectiveTransfer(int[] times) throws JMException, ClassNotFoundException {
        transferIndividuals();
        Converge(times);
    }

    private int[] getTransferSourceIndexes() {
        int[] srcTask = new int[taskNum];

//        // 随机挑选
//        for (int k = 0; k < taskNum; k++) {
//            if (inProgress[k]) {
//                srcTask[k] = k;
//                while (srcTask[k] == k)
//                    srcTask[k] = PseudoRandom.randInt(0, taskNum - 1);
//            }
//            else
//                srcTask[k] = -1;
//        }

//        // 根据scores轮盘赌
//        for (int k = 0; k < taskNum; k++) {
//            if (inProgress[k]) {
//                double[] softmax = Utils.softMaxOnlyPositive(scores[k]);
//                srcTask[k] = Utils.roulette(softmax);
//            } else {
//                srcTask[k] = -1;
//            }
//        }

        // 根据收敛方向
        for (int k = 0; k < taskNum; k++){
            int zeroCnt = 0;
            for (int i = 0; i < convergeDirections[k].length; i++){
                if (convergeDirections[k][i] == 0.0)
                    zeroCnt ++;
            }
            if (zeroCnt == convergeDirections[k].length){
                if (PseudoRandom.randDouble() < forceTransferRate) {
                    int tmp = k;
                    while (tmp == k)
                        tmp = PseudoRandom.randInt(0, taskNum-1);
                    srcTask[k] = tmp;
                }else{
                    srcTask[k] = -1;
                }
                continue;
            }

            double[] finalScore = new double[taskNum];
            double minAngle = Double.MAX_VALUE;
            int minIdx = -1;
            double[] angles = new double[taskNum];
            for (int kk = 0; kk < taskNum; kk++){
                if (k == kk) {
                    finalScore[kk] = 0;
                    continue;
                }

                double[] directionBetweenPops = Utils.vectorMinus(populationCenterPositions[kk], populationCenterPositions[k]);
                angles[kk] = Utils.calVectorAngle(convergeDirections[k], directionBetweenPops);

                if (angles[kk] < minAngle){
                    minAngle = angles[kk];
                    minIdx = kk;
                }

                if (angles[kk] > 1.57){
                    finalScore[kk] = 0;
                }else{
                    finalScore[kk] = 1.57 - angles[kk];
                    finalScore[kk] = (finalScore[kk] + scores[k][kk]) / finalScore[kk];
                }
            }

            if (PseudoRandom.randDouble() < forceTransferRate) {
                double[] softmax = Utils.softMaxOnlyPositive(finalScore);
                int tmp = k;
                while (tmp == k)
                    tmp = Utils.roulette(softmax);
                srcTask[k] = tmp;
            }else{
                srcTask[k] = minIdx;
            }
        }
        return srcTask;
    }

    private SolutionSet[] prepareTransferIndividuals() throws JMException {
        SolutionSet[] subPop = new SolutionSet[taskNum];
//        // 随机挑选个体
//        for (int k = 0; k < taskNum; k++){
//            int[] permutation = new int[transferVolume];
//            Utils.randomPermutation(permutation, transferVolume, populations[k].size());
//
//            subPop[k] = new SolutionSet(transferVolume);
//            for (int i = 0; i < transferVolume; i++){
//                Solution indiv = new Solution(populations[k].get(permutation[i]));
//                subPop[k].add(indiv);
//            }
//        }

        // 使用2-联赛选择算子
        for (int k = 0; k < taskNum; k++){
            subPop[k] = new SolutionSet(transferVolume);
            for (int i = 0; i < transferVolume; i++){
                Solution indiv = new Solution((Solution) selection.execute(populations[k]));
                subPop[k].add(indiv);
            }
        }

        return subPop;
    }

    private void transferIndividuals() throws JMException {
        SolutionSet[] subPop = prepareTransferIndividuals();

        for (int k = 0; k < taskNum; k++){
            if (transferSourceIndexes[k] >= 0) {
                int[] permutation = new int[transferVolume];
                Utils.randomPermutation(permutation, transferVolume, populations[k].size());

//                // 直接替换
//                for (int i = 0; i < transferVolume; i++) {
//                    Solution newIndiv = new Solution(subPop[transferSourceIndexes[k]].get(i));
//                    newIndiv.setProblemSet_(problemSet_.getTask(k));
//                    populations[k].replace(permutation[i], newIndiv);
//                }
                // 先交叉再替换
                for (int i = 0; i < transferVolume; i++){
                    Solution[] parents = new Solution[2];
                    parents[0] = subPop[transferSourceIndexes[k]].get(i);
                    parents[1] = populations[k].get(permutation[i]);
                    Solution[] children = (Solution[]) crossover.execute(parents);
                    Solution newIndiv = new Solution(children[0]);
                    newIndiv.setProblemSet_(problemSet_.getTask(k));
                    populations[k].replace(permutation[i], newIndiv);
                }
            }
        }
    }

    private void processConvergeImprovement(double[][] firstBestVectors, double[][] secondBestVectors){
        double[] improveModulus = new double[taskNum];
        Arrays.fill(improveModulus, 0);
        int betterCount = inProgressTaskNum;
        for (int k = 0; k < taskNum; k++) {
            double[] difference = Utils.vectorMinus(secondBestVectors[k], firstBestVectors[k]);
//            // 1. 宽松更优
//            boolean isBetter = false;
//            for (int i = 0; i < problemSet_.get(k).getNumberOfObjectives(); i++) {
//                if (difference[i] > 0) {
//                    isBetter = true;
//                    break;
//                }
//            }
//            // 2. 严格更优
//            boolean isBetter = true;
//            for (int i = 0; i < problemSet_.get(k).getNumberOfObjectives(); i++) {
//                if (difference[i] > 0) {
//                    isBetter = false;
//                    break;
//                }
//            }

            // 3. 差值向量投影值最优
            for (int i = 0; i < difference.length; i++) {
                improveModulus[k] += difference[i] * -1;
            }
            boolean isBetter;
            isBetter = improveModulus[k] > 0;

            if (!isBetter) {
                if (transferSourceIndexes[k] >= 0)
                    scores[k][transferSourceIndexes[k]] *= scoreDecreaseRate;
                transferSourceIndexes[k] = -1;
                betterCount --;


            }else{
                if (transferSourceIndexes[k] >= 0)
                    scores[k][transferSourceIndexes[k]] += scoreIncrement;
            }

            // TODO:DEBUG
            if (transferSourceIndexes[k] >= 0) {
                if (isBetter)
                    goodTransferCnt ++;
                else
                    badTransferCnt ++;
            }

        }
        if (isDRA)
            dynamicResourceAllocating(improveModulus, betterCount);
        else {
            Arrays.fill(runTimes, baseRunTime);
        }
    }

    void dynamicResourceAllocating(double[] improveModulus, int betterCount){
        double[] improveSoftmax = Utils.softMaxOnlyPositive(improveModulus);
        long totalRunTimes = betterCount * baseRunTime;
        for (int k = 0; k < taskNum; k++){
            if (inProgress[k]) {
                if (improveModulus[k] > 0) {
                    runTimes[k] = Math.max((int) Math.round(improveSoftmax[k] * totalRunTimes), 1);
                } else {
                    runTimes[k] = minRunTime;
                }
            }
            else{
                runTimes[k] = 0;
            }
        }
//        System.out.println("");
        // TODO: DEBUG
        for (int k = 0; k < taskNum; k++){
            if (inProgress[k] && runTimes[k] == 0)
                System.out.println("Exist task " + k + " run time is 0.");
        }
    }


    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        evaluations = 0;
        initState();
        initOptimizers();
//        Converge(1);
        while (evaluations < maxEvaluations) {
//            if (evaluations < maxEvaluations * 0.9) {
//                tentativeTransfer();
//                selectiveTransfer(runTimes);
//            } else{
//                Converge(1);
//            }
            calculativeConverge();
            tentativeTransfer();
            selectiveTransfer(runTimes);
//            System.out.println(evaluations + "/" + maxEvaluations);
        }
//        System.out.println("Good transfer rate: " + (double)goodTransferCnt / (goodTransferCnt + badTransferCnt));
        return populations;
    }
}