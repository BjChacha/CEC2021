package etmo.metaheuristics.matmy2;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.util.*;
import etmo.metaheuristics.matmy2.libs.*;
import etmo.util.logging.LogPopulation;

import java.util.Arrays;
import java.util.HashMap;

public class MaTMY2 extends MtoAlgorithm {
    MaTAlgorithm[] optimizers;
    SolutionSet[] populations;

    Operator crossover_;
    Operator selection_;

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

    int cntGAN = 0;
    int cntCross = 0;

    // MaTDE
    double[][] probability;

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

        crossover_ = operators_.get("crossover");
        selection_ = operators_.get("selection");

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

        // MaTDE
        probability = new double[problemSet_.size()][problemSet_.size()];
    }

    private void initPopulations() throws ClassNotFoundException, JMException {
        populations = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++){
            populations[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++){
//                Solution newSolution = new Solution(problemSet_.getTask(k));
                Solution newSolution = new Solution(problemSet_);
                newSolution.setSkillFactor(k);
                problemSet_.get(k).evaluate(newSolution);
                evaluations ++;
                populations[k].add(newSolution);
            }

//            // MaTDE
//            for (int kk = 0; kk < problemSet_.size(); kk++){
//                probability[k][kk] = 0.0;
//                scores[k][kk] = 1.0;
//            }
        }
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
//                ProblemSet pS = problemSet_.getTask(k);
                optimizers[k] = new MOEAD(problemSet_, populations[k], k);
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
                parameters.put("probability", 1.0 / problemSet_.get(k).getNumberOfVariables());
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
//                ProblemSet pS = problemSet_.getTask(k);
//                optimizers[k] = new MaOEAC(pS, populations[k]);
                optimizers[k] = new MaOEAC(problemSet_, populations[k], k);

                optimizers[k].setInputParameter("populationSize",100);
                optimizers[k].setInputParameter("maxGenerations",1000);

                parameters = new HashMap();
                parameters.put("probability", 1.0);
                parameters.put("distributionIndex", 30.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet_.get(k).getNumberOfVariables());
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

                if (((evaluations/100)*100) % (20 * taskNum * populationSize) == 0) {
//                System.out.println(totalCount+"x"+taskNum+"x"+populationSize+"="+totalCount * taskNum * populationSize);
                    LogPopulation.LogPopulation("MaTMY2", populations, problemSet_, evaluations);
                }
            }
        }
    }

    private void calculativeConverge() throws JMException, ClassNotFoundException {
        for (int k = 0; k < taskNum; k++){
//            NondominationSorting(k);
//            populationCenterPositions[k] = populations[k].getNDWeightedDecisionVector();
            populationCenterPositions[k] = populations[k].getAverageDecisionVector();
        }
        Converge(3);
        for (int k = 0; k < taskNum; k++){
//            NondominationSorting(k);
//            double[] newPosition = populations[k].getNDWeightedDecisionVector();
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

    private int[] getTransferSourceIndexes() throws JMException {
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
//                    finalScore[kk] = (1/(1 + Math.log(1 + angles[kk])) + scores[k][kk]) / (1 + Math.log(1 + angles[kk]));
                }
            }

            double transferRate = Arrays.stream(scores[k]).sum()/scores[k].length;

            if (PseudoRandom.randDouble(0, 1) < transferRate) {
                double[] softmax = Utils.softMaxOnlyPositive(finalScore);
                int tmp = k;
                int cnt = 0;
                while (tmp == k && cnt < 100) {
                    tmp = Utils.roulette(softmax);
                    cnt ++;
//                    System.out.println("Count: "+cnt);
                }
                if (cnt < 100)
                    srcTask[k] = tmp;
                else
                    srcTask[k] = -1;
            }else{
                srcTask[k] = -1;
            }
        }
//        System.out.println(evaluations + ": " + Arrays.toString(srcTask));

//        // 采用MaTDE方法
//        for (int k = 0; k < taskNum; k++){
//            srcTask[k] = findAssistTask(k);
//        }
        return srcTask;
    }

    // MaTDE
    private int findAssistTask(int task) throws JMException {
        KLD kldCalculator = new KLD(problemSet_, populations);
        double[] kld = kldCalculator.getKDL(task);
        double sum = 0;
        for (int k = 0; k < problemSet_.size(); k++){
            if (k == task)
                continue;
            probability[task][k] = 0.8 * probability[task][k] + scores[task][k] / (1 + Math.log(1 + kld[k]));
            sum += probability[task][k];
        }
        double s = 0;
        double p = PseudoRandom.randDouble();
        int idx;
        // 轮盘赌算法
        for (idx = 0; idx < problemSet_.size(); idx++) {
            if (idx == task)
                continue;
            s += probability[task][idx] / sum;
            if (s >= p)
                break;
        }
        if (idx >= problemSet_.size())
            idx = problemSet_.size() - 1;
        return idx;
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
                Solution indiv = new Solution((Solution) selection_.execute(populations[k]));
                subPop[k].add(indiv);
            }
        }

//        // 选择最优个体
//        for (int k = 0; k < taskNum; k++){
//            subPop[k] = new SolutionSet(transferVolume);
//            NondominationSorting(k);
//            for (int i = 0; i < transferVolume; i++){
//                subPop[k].add(populations[k].get(i));
//            }
//        }

        return subPop;
    }

    private void transferIndividuals() throws JMException {
        boolean[] sorted = new boolean[taskNum];

//        // LTR
//        double[][][] Ms = LTR.LTR(populations, populations[0].get(0).getDecisionVariables().length);

        SolutionSet[] subPop = prepareTransferIndividuals();

        for (int k = 0; k < taskNum; k++){
            if (transferSourceIndexes[k] >= 0) {
                int[] permutation = new int[transferVolume];
                Utils.randomPermutation(permutation, transferVolume, populations[k].size());

//                // 直接随机替换
//                for (int i = 0; i < transferVolume; i++) {
//                    Solution newIndiv = new Solution(subPop[transferSourceIndexes[k]].get(i));
//                    newIndiv.setSkillFactor(k);
//                    populations[k].replace(permutation[i], newIndiv);
//                }

                // 先交叉再随机替换
                for (int i = 0; i < transferVolume; i++){
                    Solution[] parents = new Solution[2];
                    parents[0] = subPop[transferSourceIndexes[k]].get(i);
                    parents[1] = populations[k].get(permutation[i]);
                    Solution[] children = (Solution[]) crossover_.execute(parents);
                    Solution newIndiv = new Solution(children[PseudoRandom.randInt(0, children.length - 1)]);
//                    newIndiv.setProblemSet_(problemSet_.getTask(k));
                    newIndiv.setSkillFactor(k);
                    problemSet_.get(k).evaluate(newIndiv);
                    evaluations ++;

                    populations[k].replace(permutation[i], newIndiv);
                }

//                // PCA映射
//                for (int i = 0; i < transferVolume; i++){
//                    double[][] mappingPops = Utils.MappingViaPCA(populations[transferSourceIndexes[k]].getMat(), populations[k].getMat());
//                    Solution newIndiv = new Solution(subPop[transferSourceIndexes[k]].get(i));
//                    newIndiv.setDecisionVariables(mappingPops[permutation[i]]);
//                    newIndiv.setProblemSet_(problemSet_.getTask(k));
////                    //TODO:DEBUG
////                    problemSet_.get(k).evaluate(newIndiv);
////                    if (Double.isNaN(newIndiv.getObjective(0)) || Double.isNaN(newIndiv.getObjective(1)))
////                        System.out.println("Error! Nan Transfer");
//                    populations[k].replace(permutation[i], newIndiv);
//                }

//                // AutoEncoder
//                for (int i = 0; i < transferVolume; i++){
//                    if (!sorted[k]){
//                        NondominationSorting(k);
//                        sorted[k] = true;
//                    }
//                    if (!sorted[transferSourceIndexes[k]]){
//                        NondominationSorting(transferSourceIndexes[k]);
//                        sorted[k] = true;
//                    }
//
//                    double[][] mappingPops = Utils.MappingViaAE(populations[transferSourceIndexes[k]].getMat(), populations[k].getMat(), transferVolume);
//
//                    Solution newIndiv = new Solution(populations[k].get(populationSize - 1 - i));
//                    newIndiv.setDecisionVariables(mappingPops[i]);
//                    populations[k].replace(populationSize - 1 - i, newIndiv);
//                }

//                // PCA + AutoEncoder
//                for (int i = 0; i < transferVolume; i++){
//                    // 非支配排序
//                    if (!sorted[k]){
//                        NondominationSorting(k);
//                        sorted[k] = true;
//                    }
//                    if (!sorted[transferSourceIndexes[k]]){
//                        NondominationSorting(transferSourceIndexes[k]);
//                        sorted[k] = true;
//                    }
//
//                    double[][] mappingPops = Utils.MappingViaPE(populations[transferSourceIndexes[k]].getMat(), populations[k].getMat(), transferVolume);
//
//                    Solution newIndiv = new Solution(populations[k].get(populationSize - 1 - i));
//                    newIndiv.setDecisionVariables(mappingPops[i]);
//                    //TODO:DEBUG
//                    problemSet_.get(k).evaluate(newIndiv);
//                    if (Double.isNaN(newIndiv.getObjective(0)) || Double.isNaN(newIndiv.getObjective(1)))
//                        System.out.println("Error! Nan Transfer");
//                    populations[k].replace(populationSize - 1 - i, newIndiv);
//                }


//                // GAN
//                boolean[] isGAN = {false};
//                double[][] mappingPops = Utils.MappingViaGANMixCrossover(
//                        populations[transferSourceIndexes[k]].getMat(),
//                        populations[k].getMat(),
//                        isGAN,
//                        transferVolume);
////                    System.out.println("Eval: " + evaluations + "\tGAN: " + isGAN[0]);
//                for (int i = 0; i < transferVolume; i++){
//                    Solution newIndiv;
//                    if (isGAN[0]) {
//                        newIndiv = new Solution(populations[k].get(permutation[i]));
//                        newIndiv.setDecisionVariables(mappingPops[i]);
//                        cntGAN ++;
//                    }else{
//                        cntCross ++;
//                        Solution[] parents = new Solution[2];
//                        parents[0] = new Solution(subPop[transferSourceIndexes[k]].get(i));
//                        parents[1] = new Solution(populations[k].get(permutation[i]));
//                        Solution[] children = (Solution[]) crossover_.execute(parents);
//                        newIndiv = new Solution(children[PseudoRandom.randInt(0, children.length - 1)]);
//                    }
//                    newIndiv.setSkillFactor(k);
//                    problemSet_.get(k).evaluate(newIndiv);
//                    evaluations++;
//                    populations[k].replace(permutation[i], newIndiv);
//                }

//                // GAN & Crossover
//                for (int i = 0; i < transferVolume; i++){
//                    double[] features = Utils.MappingViaGAN(
//                            subPop[transferSourceIndexes[k]].get(i), populations[k].getMat()
//                    );
//
//                    Solution newIndiv = new Solution(populations[k].get(i));
//                    newIndiv.setDecisionVariables(features);
//
//                    Solution[] parents = new Solution[3];
//                    parents[0] = newIndiv;
//                    parents[1] = new Solution(subPop[transferSourceIndexes[k]].get(i));
//                    parents[2] = new Solution(subPop[transferSourceIndexes[k]].get(i));
//
//                    // DE
////                    Solution offspring = (Solution) crossover_.execute(new Object[] {subPop[transferSourceIndexes[k]].get(i), parents});
//
//                    // SBX
//                    Solution offspring = ((Solution[]) crossover_.execute(parents))[0];
//
//                    problemSet_.get(k).evaluate(offspring);
//                    evaluations ++;
//                    populations[k].replace(i, offspring);
//                }
            }
        }
    }

    private void NondominationSorting(int k){
        Distance distance = new Distance();
        Ranking ranking = new Ranking(populations[k]);
        int remain = populations[k].size();
        int index = 0;
        SolutionSet front = null;
        populations[k].clear();

        while (remain > 0){
            front = ranking.getSubfront(index);
            distance.crowdingDistanceAssignment(front, problemSet_.get(k).getNumberOfObjectives());
            for (int i = 0; i < front.size(); i++){
                populations[k].add(front.get(i));
            }
            remain -= front.size();
            index ++;
        }
    }

    private void processConvergeImprovement(double[][] firstBestVectors, double[][] secondBestVectors){
        double[] improveModulus = new double[taskNum];
        Arrays.fill(improveModulus, 0);
        int betterCount = inProgressTaskNum;
        for (int k = 0; k < taskNum; k++) {
//            if (evaluations < maxEvaluations / 3)
//                System.out.println("DEBUG: 638");
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
                if (transferSourceIndexes[k] >= 0) {
                    scores[k][transferSourceIndexes[k]] *= scoreDecreaseRate;

                    goodTransferCnt++;
//                    System.out.println(evaluations + "/" + maxEvaluations +
//                            ": Bad transfer:" + transferSourceIndexes[k] + "->" + k +
//                            "\tTransfer Rate: " + Arrays.stream(scores[k]).sum()/scores[k].length) ;

                }
                transferSourceIndexes[k] = -1;
                betterCount --;

            }else{
                if (transferSourceIndexes[k] >= 0) {
                    scores[k][transferSourceIndexes[k]] += scoreIncrement;

                    badTransferCnt++;
//                    System.out.println(evaluations + "/" + maxEvaluations +
//                            ": Good transfer:" + transferSourceIndexes[k] + "->" + k +
//                            "\tTransfer Rate: " + Arrays.stream(scores[k]).sum()/scores[k].length);
                }
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
        System.out.println("GAN count: " + cntGAN +"/"+ (cntGAN + cntCross));
        //        System.out.println("Good transfer rate: " + (double)goodTransferCnt / (goodTransferCnt + badTransferCnt));
        return populations;
    }
}