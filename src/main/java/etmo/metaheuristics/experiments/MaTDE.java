package etmo.metaheuristics.experiments;

import etmo.core.*;
import etmo.metaheuristics.matbml.Utils;
import etmo.util.*;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.DominanceComparator;
import etmo.util.comparators.LocationComparator;
import etmo.util.sorting.SortingIdx;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

public class MaTDE extends MtoAlgorithm {
    private int populationSize;
    private int archiveSize;
    private SolutionSet[] population;
    private SolutionSet[] offspring;
    private SolutionSet[] archives;

    private int evaluations;
    private int maxEvaluations;

    Operator crossover1;
    Operator crossover2;
    Comparator dominance;

    double alpha;
    double ro;
    double shrinkRate;
    double replaceRate;

    double[][] probability;
    double[][] reward;

    // ------------
    int taskNum;
    int[][] objPos;
    double[] ideals;
    double[] nadirs;

    int leaderNum;
    int[] groups;
    int[] leaders;
    double[][] distances;
    double[][] scores;
    double[][] lastBetterRate;

    int k1 = 1;
    int k2 = 3;

    double P = 0.1;

    // DEBUG
    int[] runTimes;
    int savedTimes;
    int BTimes;
    int NTimes;
    int WTimes;

    Distance distance = new Distance();

    /**
     * Constructor
     *
     * @param problemSet The problem to be solved
     */
    public MaTDE(ProblemSet problemSet) {
        super(problemSet);
    }


    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initState();
        initPopulation();

        savedTimes = 0;
        evaluations = 0;

        while (evaluations < maxEvaluations) {
            solelyConvergence(k1);
            transferConvergence(k2);
            updateArchives();
        }
        return population;
    }

    private void updateArchives() {
        for (int k = 0; k < taskNum; k++) {
            for (int i = 0; i < populationSize; i++) {
                if (PseudoRandom.randDouble() < replaceRate)
                    putArchive(k, population[k].get(i));
            }
        }
    }

    private void initState() {
        populationSize = (Integer) getInputParameter("populationSize");
        archiveSize = (Integer) getInputParameter("archiveSize");
        maxEvaluations = (Integer) getInputParameter("maxEvaluations");

        alpha = (Double) getInputParameter("alpha");
        ro = (Double) getInputParameter("ro");
        shrinkRate = (Double) getInputParameter("shrinkRate");
        replaceRate = (Double) getInputParameter("replaceRate");

        probability = new double[problemSet_.size()][problemSet_.size()];
        reward = new double[problemSet_.size()][problemSet_.size()];

        crossover1 = operators_.get("crossover1");
        crossover2 = operators_.get("crossover2");
        dominance = new DominanceComparator();

        // ------------
        taskNum = problemSet_.size();
        population = new SolutionSet[taskNum];
        offspring = new SolutionSet[taskNum];

        objPos = new int[taskNum][2];

        ideals = new double[problemSet_.getTotalNumberOfObjs()];
        nadirs = new double[problemSet_.getTotalNumberOfObjs()];

        leaderNum = (int) Math.sqrt(taskNum);
        groups = new int[taskNum];
        leaders = new int[leaderNum];
        distances = new double[taskNum][taskNum];
        scores = new double[taskNum][taskNum];
        lastBetterRate = new double[taskNum][taskNum];

        Arrays.fill(ideals, Double.POSITIVE_INFINITY);
        for (int k = 0; k < taskNum; k++) {
            objPos[k][0] = problemSet_.get(k).getStartObjPos();
            objPos[k][1] = problemSet_.get(k).getEndObjPos();
            Arrays.fill(scores[k], 3);
            Arrays.fill(lastBetterRate[k], -1);
        }

        // DEBUG
        runTimes = new int[taskNum];
        savedTimes = 0;
        BTimes = 0;
        NTimes = 0;
        WTimes = 0;
    }

    private void initPopulation() throws JMException, ClassNotFoundException {
        population = new SolutionSet[taskNum];
        offspring = new SolutionSet[taskNum];
        archives = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
            population[k] = new SolutionSet(populationSize);
            offspring[k] = new SolutionSet(populationSize);
            archives[k] = new SolutionSet(archiveSize);
            for (int i = 0; i < populationSize; i++) {
                Solution newSolution = new Solution(problemSet_);
                newSolution.setSkillFactor(k);
                problemSet_.get(k).evaluate(newSolution);
                evaluations++;
                population[k].add(newSolution);
                putArchive(k, newSolution);
            }

            for (int kk = 0; kk < taskNum; kk++) {
                probability[k][kk] = 0.0;
                reward[k][kk] = 1.0;
            }
        }
        updateINPoint();
    }

    private void createOffspring(int task) throws JMException {
        offspring[task].clear();

        for (int i = 0; i < populationSize; i++) {
            Solution off;
            int r1 = i;
            while (r1 == i)
                r1 = PseudoRandom.randInt(0, populationSize - 1);

            Solution[] parents = new Solution[3];
            parents[0] = new Solution(population[task].get(r1));
            parents[1] = new Solution(population[task].get(i));
            parents[2] = new Solution(population[task].get(i));
            off = (Solution) crossover1.execute(new Object[]{population[task].get(i), parents});
            off.setSkillFactor(task);
            problemSet_.get(task).evaluate(off);
            evaluations++;

            int flag = dominance.compare(off, population[task].get(i));
            if (flag < 0)
                population[task].replace(i, off);
            else if (flag == 0)
                offspring[task].add(off);
        }
        unionAndSelection(population[task], offspring[task]);
    }

    private void createTransferOffspring(int task1, int task2) throws JMException {
        int assist = task2;
        offspring[task1].clear();


        double[] pBest = new double[problemSet_.get(task1).getNumberOfObjectives()];
        for (int j = 0; j <= problemSet_.get(task1).getEndObjPos() - problemSet_.get(task1).getStartObjPos(); j++) {
            pBest[j] = Double.MAX_VALUE;
            for (int i = 0; i < populationSize; i++) {
                pBest[j] = Math.min(pBest[j], population[task1].get(i).getObjective(problemSet_.get(task1).getStartObjPos() + j));
            }
        }

        for (int i = 0; i < populationSize; i++) {
            int r1 = PseudoRandom.randInt(0, populationSize - 1);
            Solution[] parents = new Solution[2];
            parents[0] = new Solution(population[assist].get(r1));
            parents[1] = new Solution(population[task1].get(i));
            Solution[] offSprings = (Solution[]) crossover2.execute(parents);
            Solution offSpring = offSprings[PseudoRandom.randInt(0, 1)];
            offSpring.setSkillFactor(task1);
            problemSet_.get(task1).evaluate(offSpring);
            evaluations++;

            int flag = dominance.compare(offSpring, population[task1].get(i));
            if (flag < 0) {
                population[task1].replace(i, offSpring);
            } else if (flag == 0) {
                offspring[task1].add(offSpring);
            }
        }
        unionAndSelection(population[task1], offspring[task1]);

        double[] pBestAfter = new double[problemSet_.get(task1).getNumberOfObjectives()];
        for (int j = 0; j <= problemSet_.get(task1).getEndObjPos() - problemSet_.get(task1).getStartObjPos(); j++) {
            pBestAfter[j] = Double.MAX_VALUE;
            for (int i = 0; i < populationSize; i++) {
                pBestAfter[j] = Math.min(pBestAfter[j], population[task1].get(i).getObjective(problemSet_.get(task1).getStartObjPos() + j));
            }
        }

        boolean isBetter = false;
        for (int i = 0; i < problemSet_.get(task1).getNumberOfObjectives(); i++) {
            if (pBestAfter[i] < pBest[i]) {
                isBetter = true;
                break;
            }
        }

        if (isBetter) {
            reward[task1][assist] /= shrinkRate;
        } else {
            reward[task1][assist] *= shrinkRate;
        }
    }

    void solelyConvergence(int times) throws JMException {
        double[] oldIdeal = ideals.clone();
        for (int t = 0; t < times; t++){
            for (int k = 0; k < taskNum; k++){
                if (runTimes[k] > populationSize * 1000)
                    continue;

                createOffspring(k);
                runTimes[k] += populationSize;
            }
        }
        updateINPoint();
        Arrays.fill(leaders, -1);
        Arrays.fill(groups, -1);
        double[] improvements = new double[taskNum];
        for (int k = 0; k < taskNum; k++){
            for (int j = objPos[k][0]; j <= objPos[k][1]; j++){
                improvements[k] += (oldIdeal[j] - ideals[j]);
            }
        }

        if (Arrays.stream(improvements).sum() > 0){
            int[] idx = SortingIdx.sort(improvements, true);
            for (int i = 0; i < leaderNum; i++)
                leaders[i] = idx[i];

            UpdateDistances();
            for (int k = 0; k < taskNum; k++) {
                int finalK = k;
                if (IntStream.of(leaders).anyMatch(x -> x == finalK))
                    continue;
                double[] finalScore = new double[leaders.length];
                for (int i = 0; i < finalScore.length; i++) {
                    finalScore[i] = scores[k][leaders[i]] > 0 ? ro * probability[k][leaders[i]] + reward[k][leaders[i]] / (1 + Math.log(1 + distances[k][leaders[i]])) : 0;
                }

                if (Arrays.stream(finalScore).sum() == 0)
                    continue;
                else {
                    groups[k] = leaders[Utils.rouletteExceptZero(finalScore)];
//				groups[k] = leaders[Utils.roulette(finalScore)];
                }
            }
        }else{
            for (int k = 0; k < taskNum; k++){
                if (PseudoRandom.randDouble() < P){
                    int idx = k;
                    while (idx == k)
                        idx = PseudoRandom.randInt(0, taskNum - 1);
                    groups[k] = idx;
                }
            }
        }
    }

    void transferConvergence(int times) throws JMException {
        for (int t = 0; t < times; t++) {
            for (int leader : leaders) {
                if (leader < 0 || runTimes[leader] > populationSize * 1000)
                    continue;

                createOffspring(leader);
                runTimes[leader] += populationSize;
            }
        }

        for (int k = 0; k < taskNum; k++){
            if (groups[k] < 0)
                continue;

            // 用leader种群来评价本任务
            double[] tmpIdeal = ideals.clone();
            double[] tmpNadir = nadirs.clone();
            int checkSize = population[groups[k]].size();
            SolutionSet BSet = new SolutionSet(checkSize);
            SolutionSet NBSet = new SolutionSet(checkSize);
            int evaluatedCount = 0;
            int leader = groups[k];

            for (int i = 0; i < population[leader].size(); i++){
                double stopP;
                if (lastBetterRate[k][groups[k]] < 0)
                    stopP = 0;
                else if (lastBetterRate[k][groups[k]] == 0)
                    stopP = (i - BSet.size()) / checkSize;
                else
                    stopP = ((i - BSet.size()) - BSet.size()/lastBetterRate[k][groups[k]]) / checkSize;
                if (PseudoRandom.randDouble() < stopP)
                    break;

                Solution tmp = new Solution(population[leader].get(i));
                tmp.setSkillFactor(k);
                problemSet_.get(k).evaluate(tmp);
                problemSet_.get(k).evaluateConstraints(tmp);
                evaluations ++;

                boolean added = false;
                for (int j = objPos[k][0]; j <= objPos[k][1]; j++) {
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
            lastBetterRate[k][groups[k]] = ((double)BSet.size()) / checkSize;

            // 判断leader种群对本任务是否有帮助
            boolean better = false;
            boolean worse = false;
            for (int j = objPos[k][0]; j <= objPos[k][1]; j++){
                if (tmpIdeal[j] < ideals[j])
                    better = true;
                if (tmpNadir[j] > nadirs[j])
                    worse = true;
            }

            // 根据情况应用leader种群
            if (better) {
                // Union and selection
//				unionAndSelection(population[k], tmpSet);
//                createTransferOffspring(k, groups[k]);

                SolutionSet offspringSet = new SolutionSet(NBSet.size());
                for (int i = 0; i < NBSet.size(); i++){
                    int r1 = PseudoRandom.randInt(0, population[k].size() - 1);
                    Solution[] parents = new Solution[2];
                    parents[0] = new Solution(population[k].get(r1));
                    parents[1] = new Solution(NBSet.get(i));
                    Solution[] offsprings = (Solution[]) crossover2.execute(parents);
                    Solution offspring = offsprings[PseudoRandom.randInt(0, 1)];
                    offspring.setSkillFactor(k);
                    problemSet_.get(k).evaluate(offspring);
                    evaluations ++;
                    offspring.setFlag(2);
                    offspringSet.add(offspring);
                }

                SolutionSet union = population[k].union(BSet);
                union = union.union(offspringSet);
                union = union.union(NBSet);

                assignFitness(union, k);
                union.sort(new LocationComparator());

                int etbc = 0;
                int itbc = 0;
                int nbc = 0;
                int other = 0;
                for (int i = 0; i < population[k].size(); i++) {
                    if (union.get(i).getFlag() == 1)
                        etbc ++;
                    else if (union.get(i).getFlag() == 2)
                        itbc ++;
                    else if (union.get(i).getFlag() == 3)
                        nbc ++;
                    else
                        other ++;
                    union.get(i).setFlag(0);
                    population[k].replace(i, union.get(i));
                }

                scores[k][groups[k]] += (((double)(etbc * 2 + itbc + nbc - other))/ population[k].size());

                BTimes ++;
                savedTimes += k2 * populationSize;
            }
            else if (worse){
                scores[k][leader] = 0;

                WTimes ++;
                savedTimes -= populationSize;
            }
            else{
                scores[k][leader] -= 1;

                NTimes ++;
                savedTimes -= populationSize;
            }
        }

        updateINPoint();
    }

    void unionAndSelection(SolutionSet P, SolutionSet Q){
        int taskId = P.get(0).getSkillFactor();
        SolutionSet union = P.union(Q);
        assignFitness(union, taskId);
        union.sort(new LocationComparator());

        P.clear();
        for (int i = 0; i < populationSize; i++)
            P.add(union.get(i));
    }

    void assignFitness(SolutionSet pop, int taskId) {
        for (int i = 0; i < pop.size(); i++)
            pop.get(i).setLocation(Integer.MAX_VALUE);

        rankSolutionOnTask(pop, taskId);
    }

    void rankSolutionOnTask(SolutionSet pop, int taskId) {
        boolean selec[] = new boolean[problemSet_.getTotalNumberOfObjs()];

        for (int i = 0; i < selec.length; i++) {
            if (i < objPos[taskId][0] || i > objPos[taskId][1])
                selec[i] = false;
            else
                selec[i] = true;
        }

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
    }

//    private void createOffspringPopulation() throws JMException {
//        for (int k = 0; k < problemSet_.size(); k++){
//            List<Solution> offspringList = new ArrayList<>();
//            double p = PseudoRandom.randDouble();
//            if (p > alpha){
//                for (int i = 0; i < populationSize; i++){
//                    Solution offSpring;
//                    int r1 = i;
//                    while (r1 == i)
//                        r1 = PseudoRandom.randInt(0, populationSize - 1);
//
//                    Solution[] parents = new Solution[3];
//                    parents[0] = new Solution(population[k].get(r1));
//                    parents[1] = new Solution(population[k].get(i));
//                    parents[2] = new Solution(population[k].get(i));
//                    offSpring = (Solution) crossover1.execute(new Object[] {population[k].get(i), parents});
//
//                    problemSet_.get(k).evaluate(offSpring);
//                    evaluations ++;
//
//                    int flag = dominance.compare(offSpring, population[k].get(i));
//                    if (flag < 0) {
//                        population[k].replace(i, offSpring);
//                    }
//                    else if (flag == 0){
//                        offspringList.add(offSpring);
//                    }
//                }
//            }
//            else{
//                int assist = findAssistTask(k);
//
//                double[] pBest = new double[problemSet_.get(k).getNumberOfObjectives()];
//                for (int j = 0; j <= problemSet_.get(k).getEndObjPos() - problemSet_.get(k).getStartObjPos(); j++) {
//                    pBest[j] = Double.MAX_VALUE;
//                    for (int i = 0; i < populationSize; i++) {
//                        pBest[j] = Math.min(pBest[j], population[k].get(i).getObjective(problemSet_.get(k).getStartObjPos()+j));
//                    }
//                }
//
//                for (int i = 0; i < populationSize; i++){
//                    int r1 = PseudoRandom.randInt(0, populationSize - 1);
//                    Solution[] parents = new Solution[2];
//                    parents[0] = new Solution(population[assist].get(r1));
//                    parents[1] = new Solution(population[k].get(i));
//                    Solution[] offSprings = (Solution[]) crossover2.execute(parents);
//                    Solution offSpring = offSprings[PseudoRandom.randInt(0,1)];
//
//                    problemSet_.get(k).evaluate(offSpring);
//                    evaluations ++;
//
//                    int flag = dominance.compare(offSpring, population[k].get(i));
//                    if (flag < 0) {
//                        population[k].replace(i, offSpring);
//                    }
//                    else if (flag == 0){
//                        offspringList.add(offSpring);
//                    }
//                }
//
//                double[] pBestAfter = new double[problemSet_.get(k).getNumberOfObjectives()];
//                for (int j = 0; j <= problemSet_.get(k).getEndObjPos() - problemSet_.get(k).getStartObjPos(); j++) {
//                    pBestAfter[j] = Double.MAX_VALUE;
//                    for (int i = 0; i < populationSize; i++) {
//                        pBestAfter[j] = Math.min(pBestAfter[j], population[k].get(i).getObjective(problemSet_.get(k).getStartObjPos()+j));
//                    }
//                }
//
//                boolean isBetter = false;
//                for (int i = 0; i < problemSet_.get(k).getNumberOfObjectives(); i++){
//                    if (pBestAfter[i] < pBest[i]){
//                        isBetter = true;
//                        break;
//                    }
//                }
//
//                if (isBetter){
//                    reward[k][assist] /= shrinkRate;
//                }
//                else{
//                    reward[k][assist] *= shrinkRate;
//                }
//            }
//
//            // 未淘汰父代的子种群
//            SolutionSet offspringPopulation = new SolutionSet(offspringList);
//            // 与原种群合并
//            SolutionSet union = population[k].union(offspringPopulation);
//
//            // 最终选择原种群大小那么多的个体
//            int remain = populationSize;
//            // pf层级
//            int idx = 0;
//            // 设置个体目标值掩码
//            boolean[] chosen = new boolean[problemSet_.getTotalNumberOfObjs()];
//            for (int i = problemSet_.get(k).getStartObjPos(); i <= problemSet_.get(k).getEndObjPos(); i++){
//                chosen[i] = true;
//            }
//            // 非支配排序
//            Ranking ranking = new Ranking(union);
//            SolutionSet front;
//            // 原种群清空，其个体已被保留到合并种群union中
//            population[k].clear();
//            Distance distance = new Distance();
//            while ((remain > 0)) {
//                // 计算拥挤度
//                front = ranking.getSubfront(idx);
//                if (remain >= front.size()) {
//                    distance.crowdingDistanceAssignment(front, problemSet_.get(k).getNumberOfObjectives(), chosen);
//                    for (int i = 0; i < front.size(); i++) {
//                        population[k].add(front.get(i));
//                    }
//                    remain -= front.size();
//                }
//                else {
//                    distance.crowdingDistanceAssignment(front, problemSet_.get(k).getNumberOfObjectives(), chosen);
//                    front.sort(new CrowdingComparator());
//                    for (int i = 0; i < remain; i++) {
//                        population[k].add(front.get(i));
//                    }
//                    break;
//                }
//                idx ++;
//            }
//        }
//    }

    private int findAssistTask(int task) throws JMException {
        KLD kldCalculator = new KLD(problemSet_, archives);
        double[] kld = kldCalculator.getKDL(task);
        double sum = 0;
        for (int k = 0; k < taskNum; k++){
            if (k == task)
                continue;
            probability[task][k] = ro * probability[task][k] + reward[task][k] / (1 + Math.log(1 + kld[k]));
            sum += probability[task][k];
        }
        double s = 0;
        double p = PseudoRandom.randDouble();
        int idx;
        // 轮盘赌算法
        for (idx = 0; idx < taskNum; idx++) {
            if (idx == task)
                continue;
            s += probability[task][idx] / sum;
            if (s >= p)
                break;
        }
        if (idx >= taskNum)
            idx = taskNum - 1;
        return idx;
    }

    private void putArchive(int task, Solution p){
        if (archives[task].size() < archiveSize){
            archives[task].add(p);
        }
        else{
            int idx = PseudoRandom.randInt(0, archiveSize - 1);
            archives[task].replace(idx, p);
        }
    }

    void updateINPoint(){
        for (int k = 0; k < taskNum; k++){
            for (int j = objPos[k][0]; j <= objPos[k][1]; j++){
                for (int i = 0; i < population[k].size(); i++){
                    ideals[j] = Math.min(ideals[j], population[k].get(i).getObjective(j));
                    nadirs[j] = Math.max(nadirs[j], population[k].get(i).getObjective(j));
                }
            }
        }
    }

    void UpdateDistances() throws JMException {
//        // Wasserstein Distance (psedo)
//        for (int i = 0; i < taskNum - 1; i++){
//            for (int j = i + 1; j < taskNum; j++){
//                double d1 = WassersteinDistance.getWD(population[i].getMat(), population[j].getMat());
//                distances[i][j] = distances[j][i] = d1;
//            }
//        }

        // KL Diversity
        KLD kld = new KLD(problemSet_, population);
        for (int i = 0; i < taskNum; i++)
            distances[i] = kld.getKDL(i);

//        // random
//        for (int i = 0; i < taskNum_; i++)
//            Arrays.fill(distances[i], 0);
    }
}
