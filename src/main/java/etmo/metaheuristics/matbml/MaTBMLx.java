package etmo.metaheuristics.matbml;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.stream.IntStream;

import etmo.core.MtoAlgorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.metaheuristics.matbml.libs.MOEAD;
import etmo.metaheuristics.matbml.libs.MaOEAC;
import etmo.metaheuristics.matbml.libs.MaTAlgorithm;
import etmo.metaheuristics.matbml.libs.NSGAII;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.util.*;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.DominanceComparator;
import etmo.util.comparators.LocationComparator;
import etmo.util.sorting.SortingIdx;

public class MaTBMLx extends MtoAlgorithm{
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

    int[] groups;
    int[] leaders;
    int[] fails;
    int[] skips;

    int minGroupNum;
    int solelyConvergeTimes;
    int transferConvergeTimes;
    int implicitTransferNum;
    int explicitTransferNum;

    boolean isExplicit;
    boolean isImplicit;
    boolean isNS;

    String algoName;
    
    public MaTBMLx(ProblemSet problemSet){
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initState();
        initPopulations();
        initOptimizers();
        while(evaluations < maxEvaluations){
            solelyConverge(solelyConvergeTimes);
            transferConverge(transferConvergeTimes);
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
            if (improvements[k] < 1e-4)
                improvements[k] = 0;
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
                    finalScore[i] = 0.5 * scores[k][leaders[i]] * Math.exp(1 / (1 + distances[k][leaders[i]]) - 1);
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
        // Leader先收敛k2代。
        int[] runTimes = new int[taskNum];
        for (int leader: leaders)
            if (leader >= 0)
                runTimes[leader] = times;
        Converge(runTimes);

        for (int k = 0; k < taskNum; k++){
            if (groups[k] < 0 || isFinished[k])
                continue;
            transfer(k, groups[k]);
        }

        updateIdealPoint();
        updateNadirPoint();
    }

    private void initState() {
        evaluations = 0;
        populationSize = (Integer) getInputParameter("populationSize");
        maxEvaluations = (Integer) getInputParameter("maxEvaluations");
        solelyConvergeTimes = (Integer) getInputParameter("solelyConvergeTimes");
        transferConvergeTimes = (Integer) getInputParameter("transferConvergeTimes");
        implicitTransferNum = (Integer) getInputParameter("implicitTransferNum");
        explicitTransferNum = (Integer) getInputParameter("explicitTransferNum");
        algoName = (String) getInputParameter("algoName");

        isImplicit = (Boolean) getInputParameter("isImplicit");
        isExplicit = (Boolean) getInputParameter("isExplicit");
        isNS = (Boolean) getInputParameter("isNS");

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

            Arrays.fill(scores[k], 3);
        }
    }

    private void initPopulations() throws ClassNotFoundException, JMException {
        populations = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++){
            populations[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++){
                Solution individual = new Solution(problemSet_);
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
        // Wasserstein Distance (psedo)
        for (int i = 0; i < taskNum - 1; i++){
            for (int j = i + 1; j < taskNum; j++){
                double d1 = WassersteinDistance.getWD(populations[i].getMat(), populations[j].getMat());
                distances[i][j] = distances[j][i] = d1;
            }
        }

//        // KL Diversity
//        KLD kld = new KLD(problemSet_, populations_);
//        for (int i = 0; i < taskNum_; i++)
//            distances[i] = kld.getKDL(i);

//        // random
//        for (int i = 0; i < taskNum_; i++)
//            Arrays.fill(distances[i], 0);
    }


    private void transfer(int targetTask, int sourceTask) throws JMException {
        int totalSize = 0;
        if (isExplicit)
            totalSize += explicitTransferNum;
        if (isImplicit)
            totalSize += implicitTransferNum;
        SolutionSet toTransfer = new SolutionSet(totalSize);

        if (isExplicit) {
            int[] toTransferIdx = new int[explicitTransferNum];
            etmo.metaheuristics.matmy2.Utils.randomPermutation(toTransferIdx, explicitTransferNum, populations[sourceTask].size());
            for (int i = 0; i < explicitTransferNum; i++) {
                Solution toTransferSolution = new Solution(populations[sourceTask].get(toTransferIdx[i]));
                toTransferSolution.setSkillFactor(targetTask);
                toTransferSolution.resetObjective();
                problemSet_.get(targetTask).evaluate(toTransferSolution);

                if (isNS) {
                    toTransfer.add(toTransferSolution);
                }else{
                    int flag = (new DominanceComparator()).compare(toTransferSolution, populations[targetTask].get(toTransferIdx[i]));
                    if (flag == -1){
                        populations[targetTask].replace(toTransferIdx[i], toTransferSolution);
                    }
                }
            }
        }
        if (isImplicit){
            int[] toTransferIdx = new int[explicitTransferNum];
            etmo.metaheuristics.matmy2.Utils.randomPermutation(toTransferIdx, explicitTransferNum, populations[sourceTask].size());
            for (int i = 0; i < explicitTransferNum; i++) {
                Solution[] parent = new Solution[2];
                parent[0] = populations[targetTask].get(i);
                parent[1] = populations[sourceTask].get(toTransferIdx[i]);
                Solution toTransferSolution = ((Solution[]) crossover.execute(parent))[0];
                toTransferSolution.setSkillFactor(targetTask);
                toTransferSolution.resetObjective();
                problemSet_.get(targetTask).evaluate(toTransferSolution);

                if (isNS) {
                    toTransfer.add(toTransferSolution);
                }else{
                    int flag = (new DominanceComparator()).compare(toTransferSolution, populations[targetTask].get(toTransferIdx[i]));
                    if (flag == -1){
                        populations[targetTask].replace(toTransferIdx[i], toTransferSolution);
                    }
                }
            }
        }

        if (isNS) {
            SolutionSet union = populations[targetTask].union(toTransfer);
            rankSolutionOnTask(populations[targetTask], targetTask, true);

            for (int i = 0; i < populations[targetTask].size(); i++) {
                populations[targetTask].replace(i, union.get(i));
            }
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
