package etmo.metaheuristics.matde;

import etmo.core.*;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.KLD;

public class MaTDE extends MtoAlgorithm {
    // Population Size per task
    private int populationSize;
    // Archive Size per task
    private int archiveSize;
    // Population
    private SolutionSet[] population;
    private SolutionSet[] archives;

    int evaluations;
    int maxEvaluations;

    Operator crossover;

    double alpha;
    double ro;
    double shrinkRate;
    double replaceRate;

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

        probability = new double[problemSet_.size()][problemSet_.size()];
        reward = new double[problemSet_.size()][problemSet_.size()];

        crossover = operators_.get("crossover");

        evaluations = 0;

        initPopulation();
        while (evaluations < maxEvaluations){
//            long startTime = System.currentTimeMillis();
            createOffspringPopulation();
//            long endTime = System.currentTimeMillis();
//            System.out.println("eval: " + evaluations + "\tcreateOffspringPopulation: " + (endTime-startTime) + " ms.");

//            startTime = System.currentTimeMillis();
            updateArchives();
//            long endTime = System.currentTimeMillis();
//            System.out.println("eval: " + evaluations + "\ttime: " + (endTime - startTime) + " ms.");
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
            double p = PseudoRandom.randDouble();
            if (p > alpha){
                // 不迁移。子种群内部进行交叉变异。
                int betterCount = 0;
                for (int i = 0; i < populationSize; i++){
                    Solution offSpring;
                    int r1 = i;
                    while (r1 == i)
                        r1 = PseudoRandom.randInt(0, populationSize - 1);

                    // 差分变异。这里按照原算法是不完全差分（仅用了两个父体）。
                    // TODO 可以用完整版DE，或者写一个参数随机波动的差分进化。
                    Solution[] parents = new Solution[3];
                    parents[0] = population[k].get(r1);
                    parents[1] = population[k].get(i);
                    parents[2] = population[k].get(i);
                    // TODO 这里调试看一下有没有正常执行
                    offSpring = (Solution) crossover.execute(
                            new Object[] {population[k].get(i), parents});

                    problemSet_.get(k).evaluate(offSpring);
                    evaluations ++;
                    // 由于原算法是单目标，这里多目标比较用支配关系。
                    // TODO 这里可以用非支配排序，若只有一层PF，则随机选一个。
                    boolean dominated = true;
                    for (int j = 0; j < problemSet_.get(k).getNumberOfObjectives(); j++) {
                        if (offSpring.getObjective(j) > population[k].get(i).getObjective(j)) {
                            dominated = false;
                            break;
                        }
                    }
                    // 若子代比父代好，则直接替换掉父代。
                    if (dominated) {
                        population[k].replace(i, offSpring);
                        betterCount ++;
                    }
                }
//                System.out.println("No Transfer) Better count: "+betterCount);
            }
            else{
                // 迁移。找一个最合适的辅助种群迁移到目标种群。
                int assist = findAssistTask(k);
                int betterCount = 0;

                // 方案2：先记录迁移前的种群最优值。
                double pBest = Double.MAX_VALUE;
                for (int i = 0; i < populationSize; i++){
                    pBest = Math.min(pBest, population[k].get(i).getObjectiveWeightedSum());
                }

                for (int i = 0; i < populationSize; i++){
                    int r1 = PseudoRandom.randInt(0, populationSize - 1);
                    // TODO 好像没有均匀交叉。
                    // 手动均匀交叉。
                    Variable[] offVar = population[k].get(i).getDecisionVariables();
                    // 选择迁移个体基因的概率。
                    double CR = PseudoRandom.randDouble(0.1, 0.9);
                    // 确保至少有一个基因被选择。
                    int l = PseudoRandom.randInt(0, offVar.length - 1);
                    for (int j = 0; j < offVar.length; j++){
                        if (j == l || PseudoRandom.randDouble() < CR) {
                            offVar[j] = population[assist].get(r1).getDecisionVariables()[j];
                        }
                        else{
                            offVar[j] = population[k].get(i).getDecisionVariables()[j];
                        }
                    }
                    Solution offSpring = new Solution(problemSet_);
                    offSpring.setDecisionVariables(offVar);
                    // 子代生成完毕。评价。
                    problemSet_.get(k).evaluate(offSpring);
                    evaluations ++;

                    boolean dominated = true;
                    for (int j = 0; j < problemSet_.get(k).getNumberOfObjectives(); j++) {
                        if (offSpring.getObjective(j) > population[k].get(i).getObjective(j)) {
                            dominated = false;
                            break;
                        }
                    }
                    // 每次生成子代都会更新种群
                    // 若子代比父代好，则直接替换掉父代。
                    if (dominated) {
                        population[k].replace(i, offSpring);
                        betterCount ++;
                    }
                }
//                System.out.println("Transfer) Better count: "+betterCount);
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
                double pBestAfter = Double.MAX_VALUE;
                for (int i = 0; i < populationSize; i++){
                    pBestAfter = Math.min(pBest, population[k].get(i).getObjectiveWeightedSum());
                }

                if (pBestAfter < pBest){
                    reward[k][assist] /= shrinkRate;
                }
                else{
                    reward[k][assist] *= shrinkRate;
                }
            }
        }
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
            int replace = PseudoRandom.randInt(0, archiveSize - 1);
            archives[task].replace(replace, p);
        }
    }
}
