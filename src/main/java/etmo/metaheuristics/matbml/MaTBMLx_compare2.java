package etmo.metaheuristics.matbml;

import etmo.core.*;
import etmo.metaheuristics.matbml.libs.MOEAD;
import etmo.metaheuristics.matbml.libs.MaOEAC;
import etmo.metaheuristics.matbml.libs.MaTAlgorithm;
import etmo.metaheuristics.matbml.libs.NSGAII;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.crossover.TransferDECrossover;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.util.*;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.DominanceComparator;
import etmo.util.comparators.LocationComparator;
import etmo.util.sorting.SortingIdx;

import java.util.Arrays;
import java.util.HashMap;
import java.util.stream.IntStream;

public class MaTBMLx_compare2 extends MtoAlgorithm{
    MaTAlgorithm[] optimizers;
    SolutionSet[] populations;
    Operator crossover;
    Operator mutation;
    Operator selection;

    int populationSize;
    int maxEvaluations;
    int evaluations;
    int taskNum;
    int objNum;

    int[] objStart;
    int[] objEnd;
    boolean[] isFinished;
    double[] ideals;
    double[] nadirs;
    double[][] scores;
    double[][] distances;
    double[] improvements;
    double[][] transferProbability;

    int initScore;
    double betterThreshold;

    int[] groups;
    int[] leaders;
    int[] fails;
    int[] skips;

    int minGroupNum;

    int solelyConvergeTimes;
    int transferConvergeTimes;
    int convergeStep;
    double stepShrinkRate;
    double transferConvergeStep;

    double transferScale;
    int distanceType;
    int environmentSelectionType;

    String algoName;

    public MaTBMLx_compare2(ProblemSet problemSet){
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initState();
        initPopulations();
        initOptimizers();
        while(evaluations < maxEvaluations){
            solelyConverge(convergeStep);
            transferConverge(convergeStep);
        }

        return populations;
    }

    private void solelyConverge(int times) throws JMException, ClassNotFoundException {
        double[] oldIdeal = ideals.clone();
        Converge(times);
        updateIdealPoint();
        updateNadirPoint();

        Arrays.fill(leaders, -1);
        Arrays.fill(groups, -1);

        improvements = new double[taskNum];
        for (int k = 0; k < taskNum; k++){
            for (int j = objStart[k]; j <= objEnd[k]; j++)
                improvements[k] += (oldIdeal[j] - ideals[j]);
        }

        if (Arrays.stream(improvements).sum() > 0) {
            int[] idxs = SortingIdx.SortingIdx(improvements, true);
            for (int i = 0; i < leaders.length; i++)
                leaders[i] = idxs[i];

            UpdateDistances();
            for (int k = 0; k < taskNum; k++) {
                int finalK = k;
                if (IntStream.of(leaders).anyMatch(x -> x == finalK))
                    continue;
                double[] finalScore = new double[leaders.length];
                for (int i = 0; i < finalScore.length; i++) {
                    finalScore[i] = scores[k][leaders[i]] * Math.exp(1 / (1 + distances[k][leaders[i]]) - 1);
                }

                if (Arrays.stream(finalScore).sum() == 0)
                    continue;
                else {
                    groups[k] = leaders[Utils.rouletteExceptZero(finalScore)];
                }
            }
        }else{
            for (int k = 0; k < taskNum; k++){
                if (PseudoRandom.randDouble() < 0.1){
                    int idx = k;
                    while (idx == k)
                        idx = PseudoRandom.randInt(0, taskNum - 1);
                    groups[k] = idx;
                }
            }
        }
    }

    private void transferConverge(int times) throws JMException, ClassNotFoundException {
        for (int k = 0; k < taskNum; k++){
            if (groups[k] < 0 || isFinished[k])
                continue;
            transfer(k, groups[k]);
        }

        updateIdealPoint();
        updateNadirPoint();
    }

    private void transfer(int targetTask, int sourceTask) throws JMException {
        SolutionSet transferSet = new SolutionSet(populationSize);
        Solution transferIndividual = null;
        int[] perm = Utils.permutation(populations[sourceTask].size(), populations[sourceTask].size());

        for (int i = 0; i < populations[sourceTask].size() * transferScale; i++) {
            int flag = populations[sourceTask].get(i).getFlag();
            int r1 = perm[i];

//            // 标记法
//            if (flag == 1) {
//                // Implicit
//                Solution[] parents = new Solution[2];
//                parents[0] = new Solution(populations[targetTask].get(i));
//                parents[1] = new Solution(populations[sourceTask].get(r1));
//                Solution[] offsprings = (Solution[]) crossover.execute(parents);
//                transferIndividual = offsprings[PseudoRandom.randInt(0, 1)];
//                mutation.execute(transferIndividual);
//                transferIndividual.setFlag(flag * 10);
//            }
//            else if (flag == 2){
//                //Explicit
//                transferIndividual = new Solution(populations[sourceTask].get(r1));
//                transferIndividual.setFlag(flag * 10);
//            }

            // 平均法
            populations[sourceTask].get(i).setFlag(0);
            if (i < populations[sourceTask].size() * transferScale / 2) {
                // Implicit
                Solution[] parents = new Solution[2];
                parents[0] = new Solution(populations[targetTask].get(i));
                parents[1] = new Solution(populations[sourceTask].get(r1));
                Solution[] offsprings = (Solution[]) crossover.execute(parents);
                transferIndividual = offsprings[PseudoRandom.randInt(0, 1)];
                mutation.execute(transferIndividual);
                transferIndividual.setFlag(1);
            } else {
                //Explicit
                transferIndividual = new Solution(populations[sourceTask].get(r1));
                transferIndividual.setFlag(1);
            }

            transferIndividual.setSkillFactor(targetTask);
            transferIndividual.resetObjective();
            problemSet_.get(targetTask).evaluate(transferIndividual);
            evaluations++;

            // Environment Selection
            if (environmentSelectionType == 1){
                // NS
                transferSet.add(transferIndividual);
            } else {
                // DR
                int betterFlag = new DominanceComparator().compare(populations[targetTask].get(i), transferIndividual);
                if (betterFlag == 1) {
                    populations[targetTask].replace(i, transferIndividual);
                }
            }

        }
        int betterCount = 0;
        if (environmentSelectionType == 1) {
            SolutionSet union = populations[targetTask].union(transferSet);
            rankSolutionOnTask(union, targetTask, true);
            for (int i = 0; i < populations[targetTask].size(); i++) {
                if (union.get(i).getFlag() == 1) {
                    betterCount += 1;
                }
                populations[targetTask].replace(i, union.get(i));
            }
        } else{
            for (int i = 0; i < populations[targetTask].size(); i++) {
                if (populations[targetTask].get(i).getFlag() == 1) {
                    betterCount += 1;
                }
            }
        }

        if (betterCount < populationSize * betterThreshold)
            scores[targetTask][sourceTask] = 0;

//        System.out.println("Explicit: " + eCount + "\tImplicit: " + iCount);
    }

    private void initState() {
        evaluations = 0;
        populationSize = (Integer) getInputParameter("populationSize");
        maxEvaluations = (Integer) getInputParameter("maxEvaluations");
        solelyConvergeTimes = (Integer) getInputParameter("solelyConvergeTimes");
        transferConvergeTimes = (Integer) getInputParameter("transferConvergeTimes");
        distanceType = (Integer) getInputParameter("distanceType");
        environmentSelectionType = (Integer) getInputParameter("environmentSelectionType");
        transferScale = (Double) getInputParameter("transferScale");
        algoName = (String) getInputParameter("algoName");

        initScore = (Integer) getInputParameter("initScore");
        betterThreshold = (Double) getInputParameter("betterThreshold");

        convergeStep = (Integer) getInputParameter("convergeStep");
        stepShrinkRate = 0.99;
        transferConvergeStep = 2 * convergeStep;

        crossover = operators_.get("crossover");
        mutation = operators_.get("mutation");
        selection = operators_.get("selection");

        taskNum = problemSet_.size();
        objNum= problemSet_.getTotalNumberOfObjs();
        minGroupNum = (int) Math.sqrt(taskNum);
        isFinished = new boolean[taskNum];
        ideals = new double[objNum];
        nadirs = new double[objNum];
        objStart = new int[taskNum];
        objEnd = new int[taskNum];

        scores = new double[taskNum][taskNum];
        distances = new double[taskNum][taskNum];

        transferProbability = new double[taskNum][taskNum];

        fails = new int[taskNum];
        skips = new int[taskNum];

        leaders = new int[minGroupNum];
        groups = new int[taskNum];

        Arrays.fill(isFinished, false);
        Arrays.fill(ideals, Double.MAX_VALUE);
        Arrays.fill(nadirs, 0);


        for (int k = 0; k < taskNum; k++){
            objStart[k] = problemSet_.get(k).getStartObjPos();
            objEnd[k] = problemSet_.get(k).getEndObjPos();

            Arrays.fill(scores[k], initScore);
            Arrays.fill(transferProbability[k], 1);
        }
    }

    private void initPopulations() throws ClassNotFoundException, JMException {
        populations = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++){
            populations[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++){
                Solution individual = new Solution(problemSet_);
//                int flag = i < populationSize / 2 ? 1 : 2;
//                individual.setFlag(flag);
                individual.setSkillFactor(k);
                problemSet_.get(k).evaluate(individual);
                evaluations ++;
                populations[k].add(individual);
            }
        }
        updateIdealPoint();
        updateNadirPoint();
    }

    private void initOptimizers() throws JMException, ClassNotFoundException {
        optimizers = new MaTAlgorithm[taskNum];
        if (algoName.equalsIgnoreCase("MOEAD")){
            Operator crossover;
            Operator mutation;

            HashMap parameters;

            for (int k = 0; k < taskNum; k++) {
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
                optimizers[k] = new NSGAII(problemSet_, populations[k], k);

                optimizers[k].setInputParameter("populationSize", 100);
                optimizers[k].setInputParameter("maxEvaluations", 100 * 1000);

                parameters = new HashMap();
                parameters.put("probability", 0.9);
                parameters.put("distributionIndex", 20.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet_.getMaxDimension());
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

    private void updateIdealPoint() {
        for (int k = 0; k < taskNum; k++){
            for (int j = objStart[k]; j <= objEnd[k]; j++){
                for (int i = 0; i < populations[k].size(); i++){
                    ideals[j] = Math.min(ideals[j], populations[k].get(i).getObjective(j));
                }
            }
        }
    }

    private void updateNadirPoint() {
        for (int k = 0; k < taskNum; k++){
            for (int j = objStart[k]; j <= objEnd[k]; j++){
                for (int i = 0; i < populations[k].size(); i++){
                    nadirs[j] = Math.max(nadirs[j], populations[k].get(i).getObjective(j));
                }
            }
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
                if (!isFinished[k]) {
                    isFinished[k] = !optimizers[k].step();
                    evaluations += populationSize;
                }
            }
        }
    }

    private void UpdateDistances() throws JMException {
        switch (distanceType) {
            case 0:
                // random
                for (int i = 0; i < taskNum; i++)
                    Arrays.fill(distances[i], 0);
                break;
            case 1:
                // Wasserstein Distance (psedo)
                for (int i = 0; i < taskNum - 1; i++) {
                    for (int j = i + 1; j < taskNum; j++) {
                        double d1 = WassersteinDistance.getWD(populations[i].getMat(), populations[j].getMat());
                        distances[i][j] = distances[j][i] = d1;
                    }
                }
                break;
            case 2:
                // KL Diversity
                KLD kld = new KLD(problemSet_, populations);
                for (int i = 0; i < taskNum; i++)
                    distances[i] = kld.getKDL(i);
                break;
        }
    }

    private void rankSolutionOnTask(SolutionSet pop, int taskId, boolean sorting) {
        for (int i = 0; i < pop.size(); i++)
            pop.get(i).setLocation(Integer.MAX_VALUE);

        boolean selec[] = new boolean[problemSet_.getTotalNumberOfObjs()];

        for (int i = 0; i < selec.length; i++) {
            if (i < objStart[taskId] || i > objEnd[taskId])
                selec[i] = false;
            else
                selec[i] = true;
        }
        Distance distance = new Distance();
        PORanking pr = new PORanking(pop, selec);
        int loc = 0;
        for (int i = 0; i < pr.getNumberOfSubfronts(); i++) {
            SolutionSet front = pr.getSubfront(i);
            distance.crowdingDistanceAssignment(front, problemSet_.getTotalNumberOfObjs(), selec);
            front.sort(new CrowdingComparator());
            for (int j = 0; j < front.size(); j++) {
                if (loc < front.get(j).getLocation())
                    front.get(j).setLocation(loc);
                loc++;
            }
        }

        if (sorting)
            pop.sort(new LocationComparator());
    }
}
