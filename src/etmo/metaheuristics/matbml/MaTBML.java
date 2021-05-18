package etmo.metaheuristics.matbml;

import etmo.core.*;
import etmo.metaheuristics.matbml.libs.*;
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

public class MaTBML extends MtoAlgorithm {
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
    int implicitTransferNum_;
    // implicit transfer probability
    double P_;

    double CBRate;
    double PBRate;

    int[] objStart_;
    int[] objEnd_;
    boolean[] isFinished_;
    double[] ideals_;
    double[] nadirs_;

    double[][] scores;
    double[][] distances;

    int[] groups_;
    int[] leaders_;
    int[] fails;
    int[] skips;

    // DEBUG
    Map<String, Integer> CBD = new HashMap<>();
    Map<String, Integer> PBD = new HashMap<>();
    Map<String, Integer> ND = new HashMap<>();
    Map<String, Integer> WD = new HashMap<>();

    int[] leaderTimes;
    int[] transferredTimes;

    int[][] CBTimes;
    int[][] PBTimes;
    int[][] NDTImes;
    int[][] WDTimes;

    public MaTBML(ProblemSet problemSet){
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

        return populations_;
    }

    private void solelyConverge(int k1) throws JMException, ClassNotFoundException {
        double[] oldIdeal = ideals_.clone();
        Converge(k1);
        updateIdealPoint();
        updateNadirPoint();

        leaders_ = new int[minGroupNum_];

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
        double[] improvements = new double[taskNum_];
        // 2.1 取前k个improvement最大的task作为leader
        for (int k = 0; k < taskNum_; k++){
            for (int j = objStart_[k]; j <= objEnd_[k]; j++)
                improvements[k] += (oldIdeal[j] - ideals_[j]);
        }
        int[] idxs = SortingIdx.SortingIdx(improvements, true);
        for (int i = 0; i < leaders_.length; i++)
            leaders_[i] = idxs[i];

        Arrays.fill(groups_, -1);

        // 2.2 其余task成员leader进组
        UpdateDistances();
        for (int k = 0; k < taskNum_; k++){
            int finalK = k;
            if (IntStream.of(leaders_).anyMatch(x -> x == finalK))
                continue;
            double[] finalScore = new double[leaders_.length];
            for (int i = 0; i < finalScore.length; i++){
                // TODO: 可调整计算公式
                finalScore[i] = scores[k][leaders_[i]] * Math.exp(1 / (1 + distances[k][leaders_[i]]) - 1);
//                // TODO: 可调整最低阈值
//                if (finalScore[i] < 1)
//                    finalScore[i] = 0;
            }
            if (Arrays.stream(finalScore).sum() == 0)
                continue;
            else{
                groups_[k] = leaders_[Utils.roulette(finalScore)];
            }
        }
    }

    private void UpdateDistances() throws JMException {
//        // Wasserstein Distance (psedo)
//        for (int i = 0; i < taskNum_ - 1; i++){
//            for (int j = i + 1; j < taskNum_; j++){
//                double d1 = WassersteinDistance.getWD(populations_[i].getMat(), populations_[j].getMat());
//                distances[i][j] = distances[j][i] = d1;
//            }
//        }

//        // KL Diversity
//        KLD kld = new KLD(problemSet_, populations_);
//        for (int i = 0; i < taskNum_; i++)
//            distances[i] = kld.getKDL(i);
    }

    private void transferConverge(int k2) throws JMException, ClassNotFoundException {
        // Leader先收敛k2代。
        int[] runTimes = new int[taskNum_];
        for (int leader: leaders_)
            if (leader >= 0)
                runTimes[leader] = k2;
        Converge(runTimes);

        for (int k = 0; k < taskNum_; k++){
            // 如果该任务是leader或者已经收敛完成，则跳过
            if (groups_[k] < 0 || isFinished_[k])
                continue;

            // 如果该任务之前失败过，则跳过
            if (skips[k] > 0){
                skips[k] -= 1;
                continue;
            }

            double[] tmpIdeal = ideals_.clone();
            double[] tmpNadir = nadirs_.clone();
            SolutionSet tmpSet = new SolutionSet(populationSize_);
            for (int i = 0; i < populations_[groups_[k]].size(); i++) {
                Solution tmp = new Solution(populations_[groups_[k]].get(i));
                tmp.setSkillFactor(k);
                problemSet_.get(k).evaluate(tmp);
                evaluations_++;
                for (int j = objStart_[k]; j <= objEnd_[k]; j++) {
                    tmpIdeal[j] = Math.min(tmpIdeal[j], tmp.getObjective(j));
                    tmpNadir[j] = Math.max(tmpNadir[j], tmp.getObjective(j));
                }
                tmpSet.add(tmp);
            }

            boolean completelyBetter = true;
            boolean partlyBetter = false;
            boolean isWorse = false;
            for (int j = objStart_[k]; j <= objEnd_[k]; j++) {
                if (tmpIdeal[j] < ideals_[j])
                    // 只要有一个好，就是partly better
                    partlyBetter = true;
                else
                    // 只要有一个不好，就不是complete better
                    completelyBetter = false;
                if (tmpNadir[j] > nadirs_[j])
                    // 只要有一个更差，就是worse
                    isWorse = true;
            }

            String key = groups_[k] + "->" + k;
            if (completelyBetter) {
//                // apply the whole population
//                for (int i = 0; i < tmpSet.size(); i++) {
//                    populations_[k].replace(i, tmpSet.get(i));
//                }

                // Union and selection
                SolutionSet union = populations_[k].union(tmpSet);
                rankSolutionOnTask(union, k, true);
                for (int i = 0; i < populations_[k].size(); i++) {
                    populations_[k].replace(i, union.get(i));
                }

                scores[k][groups_[k]] *= CBRate;
                fails[k] = 0;

                CBTimes[k][groups_[k]] += 1;
                if (CBD.containsKey(key))
                    CBD.put(key, CBD.get(key) + 1);
                else
                    CBD.put(key, 1);
            }
            else if (partlyBetter) {
                // Union and selection
                SolutionSet union = populations_[k].union(tmpSet);
                rankSolutionOnTask(union, k, true);
                for (int i = 0; i < populations_[k].size(); i++) {
                    populations_[k].replace(i, union.get(i));
                }

                scores[k][groups_[k]] *= PBRate;
                fails[k] = 0;

                PBTimes[k][groups_[k]] += 1;
                if (PBD.containsKey(key))
                    PBD.put(key, PBD.get(key) + 1);
                else
                    PBD.put(key, 1);
            }
            else if (!isWorse){
                if (PseudoRandom.randDouble() < P_) {
                    partlyImplicitTransfer(groups_[k], k, implicitTransferNum_);
                }
                fails[k] += 1;
                skips[k] = fails[k] - 1;

                if (ND.containsKey(key))
                    ND.put(key, ND.get(key) + 1);
                else
                    ND.put(key, 1);

                scores[k][groups_[k]] /= CBRate;
            }
            else{
                fails[k] += 1;
                skips[k] = fails[k] - 1;

                if (WD.containsKey(key))
                    WD.put(key, WD.get(key) + 1);
                else
                    WD.put(key, 1);

                scores[k][groups_[k]] = 0;
            }
        }

        updateIdealPoint();
        updateNadirPoint();
    }

    void partlyImplicitTransfer(int src, int trg, int num) throws JMException {
        SolutionSet transferPopulation = new SolutionSet(num);
        // use selector to pick a number of individuals.
        for (int i = 0; i < num; i++) {
            Solution newOne = new Solution((Solution) selection_.execute(populations_[src]));
            transferPopulation.add(newOne);
        }
        // produce offspring
        SolutionSet offspring = new SolutionSet(populationSize_);
        for (int i = 0; i < populations_[src].size(); i++){
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
        implicitTransferNum_ = (Integer) getInputParameter("implicitTransferNum");
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

        fails = new int[taskNum_];
        skips = new int[taskNum_];

        groups_ = new int[taskNum_];

//        Arrays.fill(isFinished_, false);
        Arrays.fill(ideals_, Double.MAX_VALUE);
//        Arrays.fill(nadirs_, 0);

        for (int k = 0; k < taskNum_; k++){
            objStart_[k] = problemSet_.get(k).getStartObjPos();
            objEnd_[k] = problemSet_.get(k).getEndObjPos();

            Arrays.fill(scores[k], 1);
        }

        CBRate = 1.2;
        PBRate = 1.1;

        // DEBUG
        leaderTimes = new int[taskNum_];
        transferredTimes = new int[taskNum_];
        CBTimes = new int[taskNum_][taskNum_];
        PBTimes = new int[taskNum_][taskNum_];
        NDTImes = new int[taskNum_][taskNum_];
        WDTimes = new int[taskNum_][taskNum_];
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
