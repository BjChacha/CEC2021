package etmo.metaheuristics.matbml;

import etmo.core.*;
import etmo.metaheuristics.matbml.libs.MOEAD;
import etmo.metaheuristics.matbml.libs.MaOEAC;
import etmo.metaheuristics.matbml.libs.MaTAlgorithm;
import etmo.metaheuristics.matbml.libs.NSGAII;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PORanking;
import etmo.util.PseudoRandom;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.LocationComparator;
import etmo.util.sorting.SortingIdx;

import java.util.*;
import java.util.stream.IntStream;

import static etmo.metaheuristics.matbml.libs.Utils.randomPermutation;

public class MaTBML22 extends MtoAlgorithm {
    MaTAlgorithm[] optimizers_;
    SolutionSet[] populations_;
    Operator crossover_;
    Operator mutation_;
    Operator selection_;

    String algoName_;
    int populationSize_;
    int maxEvaluations_;
    int evaluations_;
    int taskNum_;
    int objNum_;
    int minGroupNum_;
    int k1_;
    int k2_;
    int[][] transferNum;
    double P_;
    int transferCount;

    int[] objStart_;
    int[] objEnd_;
    boolean[] isFinished_;
    double[] ideals_;
    double[] nadirs_;
    double[][] scores;
    double[][] distances;
    double[][] transferScale;
    double[] improvements;

    int[] groups_;
    int[] leaders_;
    int[] fails;
    int[] skips;

    // DEBUG
    Map<String, Integer> CBD = new HashMap<>();
    Map<String, Integer> PBD = new HashMap<>();
    Map<String, Integer> NDB = new HashMap<>();
    Map<String, Integer> NDW = new HashMap<>();
    Map<String, Integer> WD = new HashMap<>();
    double[][] proceed;
    int proceedIdx;

    int[] leaderTimes;
    int[] transferredTimes;

    int CBTimes;
    int PBTimes;
    int NDTimes;
    int WDTimes;
    int savedEvalTimes;

    ArrayList<Integer> originalIndividualCount = new ArrayList<>();
    ArrayList<Integer> explicitTransferCount = new ArrayList<>();
    ArrayList<Integer> implicitTransferCount = new ArrayList<>();
    ArrayList<Integer> otherTransferCount = new ArrayList<>();

    public MaTBML22(ProblemSet problemSet){
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        // Initialization
        initState();
        initPopulations();
        initOptimizers();

        // Algorithm execution
        while (evaluations_ < maxEvaluations_){
            solelyConverge(k1_);
            transferConverge(k2_);
        }

        //DEBUG
//        for (int i = 0; i < taskNum_; i++) {
//            System.out.println(i + ": " + Arrays.toString(proceed[i]));
//        }
//        System.out.println("origin individual: " + originalIndividualCount.toString());
//        System.out.println("explicit transfer: " + explicitTransferCount.toString());
//        System.out.println("implicit transfer: " + implicitTransferCount.toString());
//        System.out.println("other  individual: " + otherTransferCount.toString());

        return populations_;
    }

    private void solelyConverge(int k1) throws JMException, ClassNotFoundException {
        double[] oldIdeal = ideals_.clone();
        Converge(k1);
        updateIdealPoint();
        updateNadirPoint();

        Arrays.fill(leaders_, -1);
        Arrays.fill(groups_, -1);
////         1. 先分组后选leader
//        Arrays.fill(groups_, -1);
//        List<List<Integer>> clusters = Utils.WDGrouping(populations_, minGroupNum_, taskNum_);
//
//////        // DEBUG
////        for (int c = 0; c < clusters.size(); c++)
////            System.out.println(evaluations_ + ": " + clusters.get(c));
//
//        for (int g = 0; g < clusters.size(); g++) {
//            double maxImprovement = 0;
//            leaders_[g] = -1;
//            for (int k: clusters.get(g)) {
//                double improvement = 0;
//                for (int j = objStart_[k]; j <= objEnd_[k]; j++) {
//                    improvement += (oldIdeal[j] - ideals_[j]);
//                }
//                if (improvement > maxImprovement) {
//                    leaders_[g] = k;
//                    maxImprovement = improvement;
//                }
//            }
//            for (int k: clusters.get(g)) {
//                if (k != leaders_[g])
//                    groups_[k] = leaders_[g];
//            }
//        }

        // 2. 先选leader后分组
        improvements = new double[taskNum_];
        // 2.1 取前k个improvement最大的task作为leader
        for (int k = 0; k < taskNum_; k++){
            for (int j = objStart_[k]; j <= objEnd_[k]; j++)
                improvements[k] += (oldIdeal[j] - ideals_[j]);
            if (improvements[k] < 1e-4)
                improvements[k] = 0;
        }
        if (Arrays.stream(improvements).sum() > 0) {
            int[] idxs = SortingIdx.SortingIdx(improvements, true);
            for (int i = 0; i < leaders_.length; i++)
                leaders_[i] = idxs[i];

            // 2.2 其余task成员leader进组
            UpdateDistances();
            for (int k = 0; k < taskNum_; k++) {
                int finalK = k;
                if (IntStream.of(leaders_).anyMatch(x -> x == finalK))
                    continue;
                double[] finalScore = new double[leaders_.length];
                for (int i = 0; i < finalScore.length; i++) {
//                double factor = scores[k][leaders_[i]] > 0 ? scores[k][leaders_[i]] : 0.5 * Math.pow(2, scores[k][leaders_[i]]);
                    // TODO: 调整计算公式
//                finalScore[i] = factor * Math.exp(1 / (1 + distances[k][leaders_[i]]) - 1);
                    finalScore[i] = 0.5 * scores[k][leaders_[i]] * Math.exp(1 / (1 + distances[k][leaders_[i]]) - 1);
//                finalScore[i] = scores[k][leaders_[i]] * Math.exp(1 / (1 + distances[k][leaders_[i]]) - 1);
                }

                if (Arrays.stream(finalScore).sum() == 0)
                    continue;
                else {
                    groups_[k] = leaders_[Utils.rouletteExceptZero(finalScore)];
//                groups_[k] = leaders_[Utils.roulette(finalScore)];
                }
            }
        }else{
            for (int k = 0; k < taskNum_; k++){
                if (PseudoRandom.randDouble() < 0.1){
                    int idx = k;
                    while (idx == k)
                        idx = PseudoRandom.randInt(0, taskNum_ - 1);
                    groups_[k] = idx;
                }
            }
        }

//        // DEBUG
//        System.out.println(evaluations_ + ": " + Arrays.toString(groups_));
    }

    private void transferConverge(int k2) throws JMException, ClassNotFoundException {
        // Leader先收敛k2代。
        int[] runTimes = new int[taskNum_];
        for (int leader: leaders_)
            if (leader >= 0)
                runTimes[leader] = k2;
        Converge(runTimes);

        // 开始迁移
        for (int k = 0; k < taskNum_; k++) {
            // 如果该任务是leader或者已经收敛完成，则跳过
            if (groups_[k] < 0 || isFinished_[k])
                continue;

//            // 如果该任务之前失败过，则跳过
//            if (skips[k] > 0){
//                skips[k] -= 1;
//                continue;
//            }

            double[] tmpIdeal = ideals_.clone();
            double[] tmpNadir = nadirs_.clone();
            SolutionSet candidates = new SolutionSet(transferNum[k][groups_[k]]);

            int evaluatedCount = 0;
            int[] toTransferIdx = new int[transferNum[k][groups_[k]]];
            randomPermutation(toTransferIdx, toTransferIdx.length, populations_[groups_[k]].size());
            for (int i = 0; i < transferNum[k][groups_[k]] && i <  populations_[k].size(); i++) {
                Solution toTransferSolution;
                if (i < (int) (transferNum[k][groups_[k]] * transferScale[k][groups_[k]])) {
                    // 显式迁移
                    toTransferSolution = new Solution(populations_[groups_[k]].get(i));
                    toTransferSolution.setFlag(1);

                } else {
                    // 隐式迁移
                    Solution[] parents = new Solution[2];
                    parents[0] = populations_[k].get(toTransferIdx[i]);
                    parents[1] = populations_[groups_[k]].get(i);
                    Solution[] children = (Solution[]) crossover_.execute(parents);
                    toTransferSolution = new Solution(children[PseudoRandom.randInt(0, children.length - 1)]);
                    mutation_.execute(toTransferSolution);
                    toTransferSolution.setFlag(2);
                }

                toTransferSolution.setSkillFactor(k);
                toTransferSolution.resetObjective();
                problemSet_.get(k).evaluate(toTransferSolution);
                evaluations_++;

//                // 直接替换
//                int flag = (new DominanceComparator()).compare(toTransferSolution, populations_[k].get(toTransferIdx[i]));
//                if (flag == -1) {
//                    populations_[k].replace(toTransferIdx[i], toTransferSolution);
//                }

                boolean better = false;
                for (int j = objStart_[k]; j <= objEnd_[k]; j++){
                    if (toTransferSolution.getObjective(j) < ideals_[j]){
                        better = true;
                        break;
                    }
                }
                if (better){
                    candidates.add(toTransferSolution);
                }
            }

            SolutionSet union = populations_[k].union(candidates);
            rankSolutionOnTask(union, k, true);
            for (int i = 0; i < populations_[k].size(); i++) {
                union.get(i).setFlag(0);
                populations_[k].replace(i, union.get(i));
            }

            int eSize = 0;
            int iSize = 0;
            for (int i = 0; i < populations_[k].size(); i++){
                int flag = populations_[k].get(i).getFlag();
                switch (flag){
                    case 1:
                        eSize += 1;
                        break;
                    case 2:
                        iSize += 1;
                        break;
                }
            }

//            transferNum[k][groups_[k]] = (int)(0.1 * transferNum[k][groups_[k]] + 0.9 * (eSize + iSize));
//            if (transferNum[k][groups_[k]] > 0)
//                transferScale[k][groups_[k]] = Math.max((double)1 / transferNum[k][groups_[k]], 0.1 * transferScale[k][groups_[k]] + 0.9 * ((double) eSize / (eSize + iSize + 1e-4)));
//            else
//                transferScale[k][groups_[k]] = 0;
//            System.out.println(k + "->" + groups_[k] + " transfer num: " + transferNum[k][groups_[k]]);
//            System.out.println(k + "->" + groups_[k] + " transfer scale: " + transferScale[k][groups_[k]]);
        }

        transferCount ++;
        if (transferCount >= Math.sqrt(maxEvaluations_ / populationSize_ / taskNum_))
            resetTransferParameters();

        updateIdealPoint();
        updateNadirPoint();
    }

    private void UpdateDistances() throws JMException {
        // Wasserstein Distance (psedo)
        for (int i = 0; i < taskNum_ - 1; i++){
            for (int j = i + 1; j < taskNum_; j++){
                double d1 = WassersteinDistance.getWD(populations_[i].getMat(), populations_[j].getMat());
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

    void partlyImplicitTransfer(int src, int trg, int num) throws JMException {
        // produce offspring
        SolutionSet offspring = new SolutionSet(populationSize_);
        for (int i = 0; i < num; i++){
            Solution[] parents = new Solution[2];
            parents[0] = populations_[trg].get(i);
            parents[1] = populations_[src].get(i);
            Solution[] children = (Solution[]) crossover_.execute(parents);
            Solution child = new Solution(children[PseudoRandom.randInt(0, children.length - 1)]);
            mutation_.execute(child);
            child.setSkillFactor(trg);
            problemSet_.get(trg).evaluate(child);
            evaluations_++;
            offspring.add(child);
        }
        // merge parents and children and process environment selection
        SolutionSet union = populations_[trg].union(offspring);
        rankSolutionOnTask(union, trg, true);
        for (int i = 0; i < populations_[trg].size(); i++) {
            populations_[trg].replace(i, union.get(i));
        }
    }

    private void rankSolutionOnTask(SolutionSet pop, int taskId, boolean sorting) {
        for (int i = 0; i < pop.size(); i++)
            pop.get(i).setLocation(Integer.MAX_VALUE);

        boolean selec[] = new boolean[problemSet_.getTotalNumberOfObjs()];

        for (int i = 0; i < selec.length; i++) {
            if (i < objStart_[taskId] || i > objEnd_[taskId])
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

    private SolutionSet[] prepareRepresentatives(int repNum){
        SolutionSet[] rSet = new SolutionSet[taskNum_];
        for (int k = 0; k < taskNum_; k++){
            rankSolutionOnTask(populations_[k], k, false);
            rSet[k] = new SolutionSet(repNum);
            for (int i = 0; i < populations_[k].size() && rSet[k].size() < repNum; i++){
                if (populations_[k].get(i).getLocation() < repNum) {
                    Solution tmp = new Solution(populations_[k].get(i));
                    rSet[k].add(tmp);
                }
            }
        }
        return rSet;
    }

    private List<List<Integer>> Grouping(){
        // 1. Initialize groups
        List<List<Integer>> groups = new ArrayList<>();
        for (int i = 0; i < minGroupNum_; i++)
            groups.add(new ArrayList<Integer>());

        // 2. iterate begin with a random task
        int[] perm = Utils.permutation(taskNum_, taskNum_);
        for (int k: perm){
            // the score list of groups for task k
            List<Double> kthScores = new ArrayList<>();
            for (List<Integer> group: groups){
                double tmpScore = 1;
                // calculate the total score of this group
                for (int i = 0; i < group.size(); i++){
                    tmpScore += scores[k][group.get(i)];
                }
                tmpScore /= Math.max(1, group.size() - taskNum_ / minGroupNum_ + 2);
                kthScores.add(tmpScore);
            }
            double[] scoreList = kthScores.stream().mapToDouble(Double::doubleValue).toArray();
            // if the score of all groups is negative, then create a new group
            boolean needNewGroup = true;
            for (double score: scoreList){
                if (score >= 0) {
                    needNewGroup = false;
                    break;
                }
            }
            int groupIdx;
            if (needNewGroup){
                groupIdx = groups.size();
                groups.add(new ArrayList<>());
            } else {
                groupIdx = Utils.roulette(scoreList);
            }
            groups.get(groupIdx).add(k);
        }

        return groups;
    }

    private void updateIdealPoint() {
        for (int k = 0; k < taskNum_; k++){
            for (int j = objStart_[k]; j <= objEnd_[k]; j++){
                for (int i = 0; i < populations_[k].size(); i++){
                    ideals_[j] = Math.min(ideals_[j], populations_[k].get(i).getObjective(j));
                }
            }
        }
    }

    private void updateNadirPoint() {
        for (int k = 0; k < taskNum_; k++){
            for (int j = objStart_[k]; j <= objEnd_[k]; j++){
                for (int i = 0; i < populations_[k].size(); i++){
                    nadirs_[j] = Math.max(nadirs_[j], populations_[k].get(i).getObjective(j));
                }
            }
        }
    }

    private void Converge(int times) throws JMException, ClassNotFoundException {
        int[] tmpTimes = new int[taskNum_];
        Arrays.fill(tmpTimes, times);
        Converge(tmpTimes);
    }

    private void Converge(int[] times) throws JMException, ClassNotFoundException {
        for (int k = 0; k < taskNum_; k++) {
            for (int t = 0; t < times[k]; t++) {
                if (!isFinished_[k]) {
                    isFinished_[k] = !optimizers_[k].step();
                    evaluations_ += populationSize_;

//                    if (((evaluations_ / 100) * 100) % (20 * taskNum_ * populationSize_) == 0) {
//                        LogPopulation.LogPopulation("MaTBML", populations_, problemSet_, evaluations_);
//                    }
                }
            }
        }
    }

    private void initState() {
        evaluations_ = 0;
        populationSize_ = (Integer) getInputParameter("populationSize");
        maxEvaluations_ = (Integer) getInputParameter("maxEvaluations");
        k1_ = (Integer) getInputParameter("k1");
        k2_ = (Integer) getInputParameter("k2");
        P_ = (double) getInputParameter("P_");
        algoName_ = (String) getInputParameter("algoName");

        crossover_ = operators_.get("crossover");
        mutation_ = operators_.get("mutation");
        selection_ = operators_.get("selection");

        taskNum_ = problemSet_.size();
        objNum_ = problemSet_.getTotalNumberOfObjs();
        minGroupNum_ = (int) Math.sqrt(taskNum_);
        isFinished_ = new boolean[taskNum_];
        ideals_ = new double[objNum_];
        nadirs_ = new double[objNum_];
        objStart_ = new int[taskNum_];
        objEnd_ = new int[taskNum_];

        scores = new double[taskNum_][taskNum_];
        distances = new double[taskNum_][taskNum_];
        transferNum = new int[taskNum_][taskNum_];
        transferScale = new double[taskNum_][taskNum_];

        fails = new int[taskNum_];
        skips = new int[taskNum_];

        leaders_ = new int[minGroupNum_];
        groups_ = new int[taskNum_];

//        Arrays.fill(isFinished_, false);
        Arrays.fill(ideals_, Double.MAX_VALUE);
//        Arrays.fill(nadirs_, 0);

        transferCount = 0;

        for (int k = 0; k < taskNum_; k++){
            objStart_[k] = problemSet_.get(k).getStartObjPos();
            objEnd_[k] = problemSet_.get(k).getEndObjPos();

            Arrays.fill(scores[k], 3);
            Arrays.fill(transferNum[k], populationSize_);
            Arrays.fill(transferScale[k], 0.5);
        }

        // DEBUG
        leaderTimes = new int[taskNum_];
        transferredTimes = new int[taskNum_];
        CBTimes = 0;
        PBTimes = 0;
        NDTimes = 0;
        NDTimes = 0;
        WDTimes = 0;
        savedEvalTimes = 0;
        proceed = new double[taskNum_][maxEvaluations_/populationSize_/taskNum_/k1_];
        proceedIdx = 0;
    }

    private void resetTransferParameters(){
        for (int k = 0; k < taskNum_; k++) {
            for (int kk = 0; kk < taskNum_; kk++)
                transferNum[k][kk] = populations_[kk].size();
            Arrays.fill(transferScale[k], 0.5);
        }
    }

    private void initPopulations() throws JMException, ClassNotFoundException {
        populations_ = new SolutionSet[taskNum_];
        for (int k = 0; k < taskNum_; k++){
            populations_[k] = new SolutionSet(populationSize_);
            for (int i = 0; i < populationSize_; i++){
                Solution individual = new Solution(problemSet_);
                individual.setSkillFactor(k);
                problemSet_.get(k).evaluate(individual);
                evaluations_ ++;
                populations_[k].add(individual);
            }
        }
        updateIdealPoint();
        updateNadirPoint();

//        LogPopulation.LogPopulation("MaTBML", populations_, problemSet_, evaluations_);
    }

    private void initOptimizers() throws JMException, ClassNotFoundException {
        optimizers_ = new MaTAlgorithm[taskNum_];
        if (algoName_.equalsIgnoreCase("MOEAD")){
            Operator crossover;
            Operator mutation;

            HashMap parameters;

            for (int k = 0; k < taskNum_; k++) {
                optimizers_[k] = new MOEAD(problemSet_, populations_[k], k);
                optimizers_[k].setInputParameter("populationSize", 100);
                optimizers_[k].setInputParameter("maxEvaluations", 1000 * 100);

                optimizers_[k].setInputParameter("dataDirectory", "D:\\_r\\EA\\ETMO\\MTO-cec2021-\\resources\\weightVectorFiles\\moead");

                optimizers_[k].setInputParameter("T", 20);
                optimizers_[k].setInputParameter("delta", 0.9);
                optimizers_[k].setInputParameter("nr", 2);

                parameters = new HashMap();
                parameters.put("CR", 1.0);
                parameters.put("F", 0.5);
                crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet_.get(k).getNumberOfVariables());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                optimizers_[k].addOperator("crossover", crossover);
                optimizers_[k].addOperator("mutation", mutation);

                optimizers_[k].initState();
            }
        }
        else if (algoName_.equalsIgnoreCase("MaOEAC")){
            Operator crossover;
            Operator mutation;
            Operator selection;
            HashMap parameters;

            for (int k = 0; k < taskNum_; k++) {
                optimizers_[k] = new MaOEAC(problemSet_, populations_[k], k);

                optimizers_[k].setInputParameter("populationSize",100);
                optimizers_[k].setInputParameter("maxGenerations",1000);

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

                optimizers_[k].addOperator("crossover", crossover);
                optimizers_[k].addOperator("mutation", mutation);
                optimizers_[k].addOperator("selection", selection);

                optimizers_[k].initState();
            }
        }
        else if (algoName_.equalsIgnoreCase("NSGAII")){
            Operator crossover;
            Operator mutation;
            Operator selection;

            HashMap parameters;

            for (int k = 0; k < taskNum_; k++){
                optimizers_[k] = new NSGAII(problemSet_, populations_[k], k);

                optimizers_[k].setInputParameter("populationSize", 100);
                optimizers_[k].setInputParameter("maxEvaluations", 100 * 1000);

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

                optimizers_[k].addOperator("crossover", crossover);
                optimizers_[k].addOperator("mutation", mutation);
                optimizers_[k].addOperator("selection", selection);

                optimizers_[k].initState();
            }
        }
        else{
            throw new JMException("Exception in " + algoName_ + ": Not such algorithm.");
        }
    }
}
