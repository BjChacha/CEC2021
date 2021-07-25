package etmo.metaheuristics.matbml;

import etmo.core.*;
import etmo.metaheuristics.matbml.libs.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.*;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.LocationComparator;
import etmo.util.logging.LogIGD;
import etmo.util.sorting.SortingIdx;
import org.apache.commons.lang3.ArrayUtils;
import scala.Array;

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

    int[] objStart_;
    int[] objEnd_;
    boolean[] isFinished_;
    double[] ideals_;
    double[] nadirs_;
    double[][] scores;
    double[][] distances;
    double[][] lastBetterRate;
    double[][] implicitTransferP;
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
        for (int k = 0; k < taskNum_; k++){
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
            int checkSize = populations_[groups_[k]].size();
            SolutionSet BSet = new SolutionSet(checkSize);
            SolutionSet NBSet = new SolutionSet(checkSize);

            int evaluatedCount = 0;
            for (int i = 0; i < checkSize; i++) {
//                 计算中断概率，结合历史和当前信息
//                double stopP;
//                if (lastBetterRate[k][groups_[k]] < 0)
//                    stopP = 0;
//                else if (lastBetterRate[k][groups_[k]] == 0)
//                    stopP = (i - BSet.size()) / checkSize;
//                else
//                    stopP = ((i - BSet.size()) - BSet.size()/lastBetterRate[k][groups_[k]]) / checkSize;
                // DEBUG
//                System.out.println("stop P: " + stopP);
                // 如果触发中断概率，则停止评价。
//                if (PseudoRandom.randDouble() < stopP)
//                    break;

                Solution tmp = new Solution(populations_[groups_[k]].get(i));
                tmp.setSkillFactor(k);
                tmp.resetObjective();
                problemSet_.get(k).evaluate(tmp);
                evaluations_++;

                boolean added = false;
                for (int j = objStart_[k]; j <= objEnd_[k]; j++) {
                    if (tmp.getObjective(j) < tmpIdeal[j]){
                        if (!added){
                            tmp.setFlag(1);
                            BSet.add(tmp);
                            added = true;
                        }
                        tmpIdeal[j] = tmp.getObjective(j);
                    }
                    tmpNadir[j] = Math.max(tmpNadir[j], tmp.getObjective(j));
                }
                if (!added) {
                    tmp.setFlag(3);
                    NBSet.add(tmp);
                }

                evaluatedCount ++;
            }
            lastBetterRate[k][groups_[k]] = ((double)BSet.size()) / checkSize;

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

            if (partlyBetter) {
//                boolean isImplicitTransfer = false;
//                if (PseudoRandom.randDouble() < implicitTransferP[k][groups_[k]])
//                    isImplicitTransfer = true;
                boolean isImplicitTransfer = true;

                SolutionSet offspringSet = new SolutionSet(NBSet.size());
                if (isImplicitTransfer) {
                    for (int i = 0; i < NBSet.size(); i++) {
                        int r1 = PseudoRandom.randInt(0, populations_[k].size() - 1);
                        Solution[] parents = new Solution[2];
                        parents[0] = new Solution(populations_[k].get(r1));
                        parents[1] = new Solution(NBSet.get(i));
                        Solution[] offsprings = (Solution[]) crossover_.execute(parents);
                        Solution offspring = offsprings[PseudoRandom.randInt(0, 1)];
                        mutation_.execute(offspring);
                        offspring.setSkillFactor(k);
                        problemSet_.get(k).evaluate(offspring);
                        evaluations_++;
                        offspring.setFlag(2);
                        offspringSet.add(offspring);
                    }
                }

                SolutionSet union = populations_[k].union(BSet);
                union = union.union(NBSet);
                if (isImplicitTransfer)
                    union = union.union(offspringSet);

                rankSolutionOnTask(union, k, true);

                int etbc = 0;
                int itbc = 0;
                int nbc = 0;
                int other = 0;
                for (int i = 0; i < populations_[k].size(); i++) {
                    if (union.get(i).getFlag() == 1)
                        etbc ++;
                    else if (union.get(i).getFlag() == 2)
                        itbc ++;
                    else if (union.get(i).getFlag() == 3)
                        nbc ++;
                    else
                        other ++;
                    union.get(i).setFlag(0);
                    populations_[k].replace(i, union.get(i));
                }
//                System.out.println("显式迁移存活率：" + ((double)etbc / BSet.size()));
                originalIndividualCount.add(other);
                explicitTransferCount.add(etbc);
                implicitTransferCount.add(itbc);
                otherTransferCount.add(nbc);

//                // DEBUG: Explicit & Implicit transfer
//                System.out.println(
//                        evaluations_ + "\t" +
//                        groups_[k] + "\t" +
//                        k + "\t" +
//                        etbc + "\t" +
//                        itbc + "\t" +
//                        nbc + "\t" +
//                        other);

                scores[k][groups_[k]] += (((double)(etbc * 2 + itbc + nbc - other))/ populations_[k].size());
//                if (isImplicitTransfer)
//                    implicitTransferP[k][groups_[k]] = ((double) itbc) / (offspringSet.size());
//                else
//                    implicitTransferP[k][groups_[k]] *= 2;

                fails[k] = 0;

                savedEvalTimes += k2_ * populationSize_;
                PBTimes += 1;
                if (PBD.containsKey(key))
                    PBD.put(key, PBD.get(key) + 1);
                else
                    PBD.put(key, 1);

                proceed[k][proceedIdx] = groups_[k] + 1;
            }

            else if (!isWorse){
//                ArrayList<Integer> groupMembers = new ArrayList<Integer>();
//                groupMembers.add(groups_[k]);
//                for (int i = 0; i < taskNum_; i ++)
//                    if (groups_[i] == groups_[k])
//                        groupMembers.add(i);
//
//                double[] score = new double[groupMembers.size()];
//                for (int i = 0; i < score.length; i++){
//                    double f = groupMembers.get(i) == groups_[k] ? 2 : 1;
//                    score[i] = scores[k][groupMembers.get(i)] * f;
//                }
//                int src = groupMembers.get(Utils.roulette(score));

//                boolean better = false;
//                if (PseudoRandom.randDouble() < P_) {
//                    double[] tmpIdeal2 = ideals_.clone();
//
//                    partlyImplicitTransfer(groups_[k], k, implicitTransferNum_);
//
//                    for (int i = 0; i < populations_[k].size(); i++)
//                        for (int j = objStart_[k]; j <= objEnd_[k]; j++) {
//                            tmpIdeal2[j] = Math.min(tmpIdeal2[j], populations_[k].get(i).getObjective(j));
//                            if (tmpIdeal[j] < ideals_[j]) {
//                                better = true;
//                                break;
//                        }
//                    }
//                }
//
//                if (!better) {
//                    fails[k] += 1;
//                    skips[k] = fails[k] - 1;
//
//                    proceed[k][proceedIdx] = -(groups_[k] + 1.1);
//                    if (NDW.containsKey(key))
//                        NDW.put(key, NDW.get(key) + 1);
//                    else
//                        NDW.put(key, 1);
//
//                }else{
//                    fails[k] = 1;
//
//                    proceed[k][proceedIdx] = groups_[k] + 1.1;
//                    if (NDB.containsKey(key))
//                        NDB.put(key, NDB.get(key) + 1);
//                    else
//                        NDB.put(key, 1);
//                }

//                SolutionSet offspringSet = new SolutionSet(populations_[k].size());
//                for (int i = 0; i < offspringSet.size(); i++){
//                    int r1 = PseudoRandom.randInt(0, populations_[k].size() - 1);
//                    Solution[] parents = new Solution[2];
//                    parents[0] = new Solution(populations_[k].get(r1));
//                    parents[1] = new Solution(populations_[groups_[k]].get(i));
//                    Solution[] offsprings = (Solution[]) crossover_.execute(parents);
//                    Solution offspring = offsprings[PseudoRandom.randInt(0, 1)];
//                    mutation_.execute(offspring);
//                    offspring.setSkillFactor(k);
//                    problemSet_.get(k).evaluate(offspring);
//                    evaluations_ ++;
//                    offspring.setFlag(2);
//                    offspringSet.add(offspring);
//                }
//
//                SolutionSet union = populations_[k].union(offspringSet);
//
//                rankSolutionOnTask(union, k, true);
//
//                int etbc = 0;
//                int itbc = 0;
//                int nbc = 0;
//                int other = 0;
//                for (int i = 0; i < populations_[k].size(); i++) {
//                    if (union.get(i).getFlag() == 1)
//                        etbc ++;
//                    else if (union.get(i).getFlag() == 2)
//                        itbc ++;
//                    else if (union.get(i).getFlag() == 3)
//                        nbc ++;
//                    else
//                        other ++;
//                    union.get(i).setFlag(0);
//                    populations_[k].replace(i, union.get(i));
//                }

                savedEvalTimes -= evaluatedCount;
                NDTimes += 1;
                if (NDB.containsKey(key))
                    NDB.put(key, NDB.get(key) + 1);
                else
                    NDB.put(key, 1);
                scores[k][groups_[k]] = Math.max(0, scores[k][groups_[k]] - 1);
            }
            else{
                fails[k] += 1;
                skips[k] = fails[k] - 1;

                savedEvalTimes -= evaluatedCount;
                WDTimes += 1;
                if (WD.containsKey(key))
                    WD.put(key, WD.get(key) + 1);
                else
                    WD.put(key, 1);

                scores[k][groups_[k]] = 0;

                proceed[k][proceedIdx] = -(groups_[k] + 1);
            }
        }

        proceedIdx ++;

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
//                if (!isFinished_[k]) {
                    isFinished_[k] = !optimizers_[k].step();
                    evaluations_ += populationSize_;

//                    if (((evaluations_ / 100) * 100) % (20 * taskNum_ * populationSize_) == 0) {
//                        LogPopulation.LogPopulation("MaTBML", populations_, problemSet_, evaluations_);
//                    }
//                }
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
        implicitTransferP = new double[taskNum_][taskNum_];
        lastBetterRate = new double[taskNum_][taskNum_];

        fails = new int[taskNum_];
        skips = new int[taskNum_];

        leaders_ = new int[minGroupNum_];
        groups_ = new int[taskNum_];

//        Arrays.fill(isFinished_, false);
        Arrays.fill(ideals_, Double.MAX_VALUE);
//        Arrays.fill(nadirs_, 0);


        for (int k = 0; k < taskNum_; k++){
            objStart_[k] = problemSet_.get(k).getStartObjPos();
            objEnd_[k] = problemSet_.get(k).getEndObjPos();

            Arrays.fill(scores[k], 3);
            Arrays.fill(lastBetterRate[k], -1);
            Arrays.fill(implicitTransferP[k], 1);
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
