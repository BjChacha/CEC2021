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
    private int populationSize;
    private int archiveSize;
    private SolutionSet[] population;
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

        crossover1 = operators_.get("crossover1");
        crossover2 = operators_.get("crossover2");
        dominance = new DominanceComparator();

        evaluations = 0;

        initPopulation();
        LogPopulation.LogPopulation("MaTDE",population, problemSet_,evaluations);
        while (evaluations < maxEvaluations){
            createOffspringPopulation();
            updateArchives();
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
                newSolution.setSkillFactor(k);
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

    private void createOffspringPopulation() throws JMException {
        for (int k = 0; k < problemSet_.size(); k++){
            List<Solution> offspringList = new ArrayList<>();
            double p = PseudoRandom.randDouble();
            if (p > alpha){
                for (int i = 0; i < populationSize; i++){
                    Solution offSpring;
                    int r1 = i;
                    while (r1 == i)
                        r1 = PseudoRandom.randInt(0, populationSize - 1);

                    Solution[] parents = new Solution[3];
                    parents[0] = new Solution(population[k].get(r1));
                    parents[1] = new Solution(population[k].get(i));
                    parents[2] = new Solution(population[k].get(i));
                    offSpring = (Solution) crossover1.execute(new Object[] {population[k].get(i), parents});
                    offSpring.setSkillFactor(k);
                    problemSet_.get(k).evaluate(offSpring);
                    evaluations ++;

                    int flag = dominance.compare(offSpring, population[k].get(i));
                    if (flag < 0) {
                        population[k].replace(i, offSpring);
                    }
                    else if (flag == 0){
                        offspringList.add(offSpring);
                    }
                }
            }
            else{
                int assist = findAssistTask(k);

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
                    offSpring.setSkillFactor(k);
                    problemSet_.get(k).evaluate(offSpring);
                    evaluations ++;

                    int flag = dominance.compare(offSpring, population[k].get(i));
                    if (flag < 0) {
                        population[k].replace(i, offSpring);
                    }
                    else if (flag == 0){
                        offspringList.add(offSpring);
                    }
                }

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
            SolutionSet front;
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
        int idx;
        // 轮盘赌算法
        for (idx = 0; idx < problemSet_.size() - 1; idx++) {
            if (idx == task)
                continue;
            s += probability[task][idx] / sum;
            if (s >= p)
                break;
        }
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
