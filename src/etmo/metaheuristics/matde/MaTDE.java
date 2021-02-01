package etmo.metaheuristics.matde;

import etmo.core.*;
import etmo.util.*;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.DominanceComparator;
import etmo.util.logging.LogPopulation;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class MaTDE extends MtoAlgorithm {
    // Population Size per task
    private int populationSize;
    // Archive Size per task
    private int archiveSize;
    // Population
    private SolutionSet[] population;
    private SolutionSet[] archives;

    private int evaluations;
    private int maxEvaluations;

    // 不迁移时用的差分交叉
    Operator crossover1;
    // 迁移时用的SBX交叉
    Operator crossover2;
    Operator mutation;
    Operator selector;
    Comparator dominance;

    double alpha;
    double ro;
    double shrinkRate;
    double replaceRate;

    double transferThreshold;

    double[][] probability;
    double[][] reward;


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
        populationSize = (Integer) getInputParameter("populationSize");
        archiveSize = (Integer) getInputParameter("archiveSize");
        maxEvaluations = (Integer) getInputParameter("maxEvaluations");

        alpha = (Double) getInputParameter("alpha");
        ro = (Double) getInputParameter("ro");
        shrinkRate = (Double) getInputParameter("shrinkRate");
        replaceRate = (Double) getInputParameter("replaceRate");

        transferThreshold = (Double) getInputParameter("transferThreshold");

        probability = new double[problemSet_.size()][problemSet_.size()];
        reward = new double[problemSet_.size()][problemSet_.size()];

        crossover1 = operators_.get("crossover1");
        crossover2 = operators_.get("crossover2");
        mutation = operators_.get("mutation");
        selector = operators_.get("operator");
        dominance = new DominanceComparator();

        evaluations = 0;

        initPopulation();
        LogPopulation.LogPopulation("MaTDE",population, problemSet_,evaluations);
        while (evaluations < maxEvaluations){
            createOffspringPopulation();
            updateArchives();
//            for (int k = 0; k < problemSet_.size(); k++)
//                System.out.println("Archive " + k + ": " + archives[k].size());
            // 一共1000代
            if (evaluations % (problemSet_.size() * populationSize * 20) == 0){
                LogPopulation.LogPopulation("MaTDE",population, problemSet_,evaluations);
            }
        }
        return population;
    }

    private void updateArchives() {
        for (int k = 0; k < problemSet_.size(); k++) {
            for (int i = 0; i < populationSize; i++){
                if (PseudoRandom.randDouble() < replaceRate)
                    putArchive(k, population[k].get(i));
            }
        }
    }


    private void initPopulation() throws JMException, ClassNotFoundException {
        population = new SolutionSet[problemSet_.size()];
        archives = new SolutionSet[problemSet_.size()];
        for (int k = 0; k < problemSet_.size(); k++){
            population[k] = new SolutionSet(populationSize);
            archives[k] = new SolutionSet(archiveSize);
            for (int i = 0; i < populationSize; i++){
                Solution newSolution = new Solution(problemSet_);
                problemSet_.get(k).evaluate(newSolution);
                evaluations ++;
                population[k].add(newSolution);
                putArchive(k, newSolution);
            }

            for (int kk = 0; kk < problemSet_.size(); kk++){
                probability[k][kk] = 0.0;
                reward[k][kk] = 1.0;
            }
        }
    }

    private void createOffspringPopulation() throws JMException, ClassNotFoundException {
        // 每个任务单独进行
        for (int k = 0; k < problemSet_.size(); k++){
            List<Solution> offspringList = new ArrayList<Solution>();
            double p = PseudoRandom.randDouble();
            int betterCount = 0;
            if (p > alpha || evaluations <= Math.floor(maxEvaluations * transferThreshold)){
                // 不迁移。子种群内部进行交叉变异。
                for (int i = 0; i < populationSize; i++){
                    Solution offSpring;
                    int r1 = i;
                    while (r1 == i)
                        r1 = PseudoRandom.randInt(0, populationSize - 1);

                    // 差分交叉
                    Solution[] parents = new Solution[3];
                    parents[0] = new Solution(population[k].get(r1));
                    parents[1] = new Solution(population[k].get(i));
                    parents[2] = new Solution(population[k].get(i));
                    offSpring = (Solution) crossover1.execute(
                            new Object[] {population[k].get(i), parents});

                    // 新增：变异
//                    mutation.execute(offSpring);
                    problemSet_.get(k).evaluate(offSpring);
                    evaluations ++;
                    // 由于原算法是单目标，这里多目标比较用支配关系。
                    int flag = dominance.compare(offSpring, population[k].get(i));
                    if (flag == -1) {
                        population[k].replace(i, offSpring);
                    }
                    else if (flag == 0){
                        offspringList.add(offSpring);
                    }
                }
//                System.out.println("No Transfer) Better count: "+betterCount);
            }
            else{
                // 迁移。找一个最合适的辅助种群迁移到目标种群。
                int assist = findAssistTask(k);
//                int betterCount = 0;

                // 方案2：先记录迁移前的种群最优值。
                double[] pBest = new double[problemSet_.get(k).getNumberOfObjectives()];
                for (int j = 0; j <= problemSet_.get(k).getEndObjPos() - problemSet_.get(k).getStartObjPos(); j++) {
                    pBest[j] = Double.MAX_VALUE;
                    for (int i = 0; i < populationSize; i++) {
                        pBest[j] = Math.min(pBest[j], population[k].get(i).getObjective(problemSet_.get(k).getStartObjPos()+j));
                    }
                }

                for (int i = 0; i < populationSize; i++){
                    int r1 = PseudoRandom.randInt(0, populationSize - 1);
                    Solution[] parents = new Solution[2];
                    parents[0] = new Solution(population[assist].get(r1));
                    parents[1] = new Solution(population[k].get(i));
                    Solution[] offSprings = (Solution[]) crossover2.execute(parents);
                    Solution offSpring = offSprings[PseudoRandom.randInt(0,1)];

                    // 新增：变异
//                    mutation.execute(offSpring);

                    // 子代生成完毕。评价。
                    problemSet_.get(k).evaluate(offSpring);
                    evaluations ++;

                    int flag = dominance.compare(offSpring, population[k].get(i));
                    if (flag == -1) {
//                        betterCount ++;
                        if (k == 0){
                            Solution tmp = population[k].get(i);
                        }
                        population[k].replace(i, offSpring);
                    }
                    else if (flag == 0){
                        offspringList.add(offSpring);
                    }
                }
//                System.out.println("Task "+k+": transfer Better count: "+betterCount);
                // 原算法：若pBest更新，则奖励增加。
                // 这里是多目标，不好判断pBest是否更新。

                // 方案1：用子代淘汰父代次数来代替。
                // 这里可以做自适应：与不迁移时的淘汰频率比较，等。
//                if (betterCount >= Math.round(0.35 * populationSize)){
//                    reward[k][assist] /= shrink_rate;
//                }
//                else{
//                    reward[k][assist] += shrink_rate;
//                }

                // 方案2：用目标函数总和来判断。
                double[] pBestAfter = new double[problemSet_.get(k).getNumberOfObjectives()];
                for (int j = 0; j <= problemSet_.get(k).getEndObjPos() - problemSet_.get(k).getStartObjPos(); j++) {
                    pBestAfter[j] = Double.MAX_VALUE;
                    for (int i = 0; i < populationSize; i++) {
                        pBestAfter[j] = Math.min(pBestAfter[j], population[k].get(i).getObjective(problemSet_.get(k).getStartObjPos()+j));
                    }
                }

                boolean isBetter = false;
                for (int i = 0; i < problemSet_.get(k).getNumberOfObjectives(); i++){
                    if (pBestAfter[i] < pBest[i]){
                        isBetter = true;
                        break;
                    }
                }

                if (isBetter){
                    reward[k][assist] /= shrinkRate;
                }
                else{
                    reward[k][assist] *= shrinkRate;
                }
            }

            // 未淘汰父代的子种群
            SolutionSet offspringPopulation = new SolutionSet(offspringList);
            // 与原种群合并
            SolutionSet union = population[k].union(offspringPopulation);

            // 最终选择原种群大小那么多的个体
            int remain = populationSize;
            // pf层级
            int idx = 0;
            // 设置个体目标值掩码
            boolean[] chosen = new boolean[problemSet_.getTotalNumberOfObjs()];
            for (int i = problemSet_.get(k).getStartObjPos(); i <= problemSet_.get(k).getEndObjPos(); i++){
                chosen[i] = true;
            }
            // 非支配排序
            Ranking ranking = new Ranking(union);
            SolutionSet front = null;
            // 原种群清空，其个体已被保留到合并种群union中
            population[k].clear();
            Distance distance = new Distance();
            while ((remain > 0)) {
                // 计算拥挤度
                front = ranking.getSubfront(idx);
                if (remain >= front.size()) {
                    distance.crowdingDistanceAssignment(front, problemSet_.get(k).getNumberOfObjectives(), chosen);
                    for (int i = 0; i < front.size(); i++) {
                        population[k].add(front.get(i));
                    }
                    remain -= front.size();
                }
                else {
                    distance.crowdingDistanceAssignment(front, problemSet_.get(k).getNumberOfObjectives(), chosen);
                    front.sort(new CrowdingComparator());
                    for (int i = 0; i < remain; i++) {
                        population[k].add(front.get(i));
                    }
                    break;
                }
                idx ++;
            }
        }
//        System.out.println("Better Count: " + betterCount);
    }

    private int findAssistTask(int task) throws JMException {
        KLD kldCalculator = new KLD(problemSet_, archives);
        double[] kld = kldCalculator.getKDL(task);
        double sum = 0;
        for (int k = 0; k < problemSet_.size(); k++){
            if (k == task)
                continue;
            probability[task][k] = ro * probability[task][k] + reward[task][k] / (1 + Math.log(1 + kld[k]));
            sum += probability[task][k];
        }
        double s = 0;
        double p = PseudoRandom.randDouble();
        int idx = 0;
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

    private void putArchive(int task, Solution p){
        if (archives[task].size() < archiveSize){
            archives[task].add(p);
        }
        else{
            int idx = PseudoRandom.randInt(0, archiveSize - 1);
            archives[task].replace(idx, p);
        }
    }
}
