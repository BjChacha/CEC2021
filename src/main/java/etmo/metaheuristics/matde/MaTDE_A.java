package etmo.metaheuristics.matde;

import etmo.core.*;
import etmo.metaheuristics.matbml.Utils;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.*;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.DominanceComparator;
import etmo.util.logging.LogIGD;
import etmo.util.sorting.SortingIdx;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.IntStream;

public class MaTDE_A extends MtoAlgorithm {
    private int populationSize;
    private int archiveSize;
    private SolutionSet[] population;
    private SolutionSet[] offspring;
    private SolutionSet[] archives;

    private int taskNum;
    private int evaluations;
    private int maxEvaluations;

    Operator crossover1;
    Operator crossover2;
    Comparator dominance;

    double alpha;
    double ro;
    double shrinkRate;
    double replaceRate;

    int convergeStep;
    double stepShrinkRate;
    double transferConvergeStep;

    double[][] probability;
    double[][] reward;


    // IGD
    ArrayList<Double>[] igds;
    String[] pf;
    QualityIndicator[] indicators;

    /**
     * Constructor
     *
     * @param problemSet The problem to be solved
     */
    public MaTDE_A(ProblemSet problemSet) {
        super(problemSet);
    }


    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initState();
        initPopulation();
        saveIGD();

        // original
//        while (evaluations < maxEvaluations){
//            createOffspringPopulation();
//            updateArchives();
//        }

        // A1
        while (evaluations < maxEvaluations){
            // 单独收敛5次，并收集各个任务的提升度
            for (int k = 0; k < taskNum; k++){
                for (int i = 0; i < convergeStep; i++){
                    normalReproduce(k);
                    unionAndRankSelection(k);
                    updateArchives(k);
                    saveIGD(k);
                }
            }
            // 迁移，迁移源任务额外收敛5次
            for (int k = 0; k < taskNum; k++){
                int assist = findAssistTask(k);
                for (int i = 0; i < transferConvergeStep; i++){
                    normalReproduce(assist);
                    unionAndRankSelection(assist);
                    updateArchives(assist);
                    saveIGD(assist);
                }
                transferReproduce(k, assist);
                unionAndRankSelection(k);
                updateArchives(k);
                saveIGD(k);
            }
            transferConvergeStep *= stepShrinkRate;
        }

//        // A2
//        while (evaluations < maxEvaluations){
//            solelyConverge(convergeStep);
//            for (int k = 0; k < taskNum; k++){
//                int assist = findAssistTask(k);
//                transferReproduce(k, assist);
//                unionAndRankSelection(k);
//                updateArchives(k);
//                saveIGD(k);
//            }
//        }

//        double[][] IGD = new double[taskNum][];
//        for (int k = 0; k < taskNum; k++){
//            IGD[k] = igds[k].stream().mapToDouble(Double::doubleValue).toArray();
//        }
//        LogIGD.MarkLog("MaTDE_A_", problemSet_.get(0).getName(), IGD);

        return population;
    }

    private void solelyConverge(int times) throws JMException, ClassNotFoundException {
        double[] improvements = new double[taskNum];
        for (int k = 0; k < taskNum; k++){
            double[] oldIdeal = getIdealPoint(k);
            for (int t = 0; t < times; t++){
                normalReproduce(k);
                unionAndRankSelection(k);
                updateArchives(k);
                saveIGD(k);
            }
            double[] newIdeal = getIdealPoint(k);
            for (int i = 0; i < newIdeal.length; i++){
                improvements[k] += (oldIdeal[i] - newIdeal[i]);
            }
        }

        int[] idxs = SortingIdx.SortingIdx(improvements, true);
        for (int i = 0; i < (int) Math.sqrt(taskNum); i++){
            for (int t = 0; t < times; t++){
                normalReproduce(idxs[i]);
                unionAndRankSelection(idxs[i]);
                updateArchives(idxs[i]);
                saveIGD(idxs[i]);
            }
        }
    }


    private void updateArchives() {
        for (int k = 0; k < problemSet_.size(); k++) {
            updateArchives(k);
        }
    }

    private void updateArchives(int task) {
        for (int i = 0; i < populationSize; i++){
            if (PseudoRandom.randDouble() < replaceRate)
                putArchive(task, population[task].get(i));
        }
    }

    private void initState(){
        populationSize = (Integer) getInputParameter("populationSize");
        archiveSize = (Integer) getInputParameter("archiveSize");
        maxEvaluations = (Integer) getInputParameter("maxEvaluations");

        alpha = (Double) getInputParameter("alpha");
        ro = (Double) getInputParameter("ro");
        shrinkRate = (Double) getInputParameter("shrinkRate");
        replaceRate = (Double) getInputParameter("replaceRate");
        convergeStep = (Integer) getInputParameter("convergeStep");

        stepShrinkRate = 0.99;
        transferConvergeStep = 2 * convergeStep;

        probability = new double[problemSet_.size()][problemSet_.size()];
        reward = new double[problemSet_.size()][problemSet_.size()];

        crossover1 = operators_.get("crossover1");
        crossover2 = operators_.get("crossover2");
        dominance = new DominanceComparator();

        evaluations = 0;
        taskNum = problemSet_.size();

        // IGD
        pf = new String[taskNum];
        igds = new ArrayList[taskNum];
        indicators = new QualityIndicator[taskNum];
        for (int k = 0; k < pf.length; k++){
            pf[k] = "PF/StaticPF/" + problemSet_.get(k).getHType() + "_" + problemSet_.get(k).getNumberOfObjectives() + "D.pf";
            igds[k] = new ArrayList<>();
            indicators[k] = new QualityIndicator(problemSet_.get(k), pf[k]);
        }
    }

    private void initPopulation() throws JMException, ClassNotFoundException {
        population = new SolutionSet[problemSet_.size()];
        offspring = new SolutionSet[problemSet_.size()];
        archives = new SolutionSet[problemSet_.size()];
        for (int k = 0; k < problemSet_.size(); k++){
            population[k] = new SolutionSet(populationSize);
            offspring[k] = new SolutionSet(populationSize);
            archives[k] = new SolutionSet(archiveSize);
            for (int i = 0; i < populationSize; i++){
                Solution newSolution = new Solution(problemSet_);
                newSolution.setSkillFactor(k);
                newSolution.resetObjective();
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

    private void normalReproduce(int task) throws JMException {
        offspring[task].clear();
        for (int i = 0; i < populationSize; i++){
            Solution offSpring;
            int r1 = i;
            while (r1 == i)
                r1 = PseudoRandom.randInt(0, populationSize - 1);

            Solution[] parents = new Solution[3];
            parents[0] = new Solution(population[task].get(r1));
            parents[1] = new Solution(population[task].get(i));
            parents[2] = new Solution(population[task].get(i));
            offSpring = (Solution) crossover1.execute(new Object[] {population[task].get(i), parents});
            offSpring.setSkillFactor(task);
            offSpring.resetObjective();
            problemSet_.get(task).evaluate(offSpring);
            evaluations ++;

            int flag = dominance.compare(offSpring, population[task].get(i));
            if (flag < 0) {
                population[task].replace(i, offSpring);
            }
            else if (flag == 0){
                offspring[task].add(offSpring);
            }
        }
    }

    private void transferReproduce(int targetTask, int sourceTask) throws JMException {
        offspring[targetTask].clear();
        for (int i = 0; i < populationSize; i++) {
            int r1 = PseudoRandom.randInt(0, populationSize - 1);
            Solution[] parents = new Solution[2];
            parents[0] = new Solution(population[sourceTask].get(r1));
            parents[1] = new Solution(population[targetTask].get(i));
            Solution[] offSprings = (Solution[]) crossover2.execute(parents);
            Solution offSpring = offSprings[PseudoRandom.randInt(0, 1)];
            offSpring.setSkillFactor(targetTask);
            offSpring.resetObjective();
            problemSet_.get(targetTask).evaluate(offSpring);
            evaluations++;

            int flag = dominance.compare(offSpring, population[targetTask].get(i));
            if (flag < 0) {
                population[targetTask].replace(i, offSpring);
            } else if (flag == 0) {
                offspring[targetTask].add(offSpring);
            }
        }
    }

    private void createOffspringPopulation() throws JMException {
        for (int k = 0; k < problemSet_.size(); k++){
            double p = PseudoRandom.randDouble();
            if (p > alpha){
                normalReproduce(k);
            } else{
                int assist = findAssistTask(k);
                double[] pBest = getIdealPoint(k);
                transferReproduce(k, assist);
                double[] pBestAfter = getIdealPoint(k);

                boolean isBetter = false;
                for (int i = 0; i < pBest.length; i++){
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
            unionAndRankSelection(k);
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

    private double[] getIdealPoint(int task){
        double[] res = new double[problemSet_.get(task).getNumberOfObjectives()];
        for (int j = 0; j <= problemSet_.get(task).getEndObjPos() - problemSet_.get(task).getStartObjPos(); j++) {
            res[j] = Double.MAX_VALUE;
            for (int i = 0; i < populationSize; i++) {
                res[j] = Math.min(res[j], population[task].get(i).getObjective(problemSet_.get(task).getStartObjPos()+j));
            }
        }
        return res;
    }

    private void unionAndRankSelection(int task){
        SolutionSet union = population[task].union(offspring[task]);
        // 最终选择原种群大小那么多的个体
        int remain = populationSize;
        // pf层级
        int idx = 0;
        // 设置个体目标值掩码
        boolean[] chosen = new boolean[problemSet_.getTotalNumberOfObjs()];
        for (int i = problemSet_.get(task).getStartObjPos(); i <= problemSet_.get(task).getEndObjPos(); i++){
            chosen[i] = true;
        }
        // 非支配排序
        Ranking ranking = new Ranking(union);
        SolutionSet front;
        // 原种群清空，其个体已被保留到合并种群union中
        population[task].clear();
        Distance distance = new Distance();
        while ((remain > 0)) {
            // 计算拥挤度
            front = ranking.getSubfront(idx);
            if (remain >= front.size()) {
                distance.crowdingDistanceAssignment(front, problemSet_.get(task).getNumberOfObjectives(), chosen);
                for (int i = 0; i < front.size(); i++) {
                    population[task].add(front.get(i));
                }
                remain -= front.size();
            }
            else {
                distance.crowdingDistanceAssignment(front, problemSet_.get(task).getNumberOfObjectives(), chosen);
                front.sort(new CrowdingComparator());
                for (int i = 0; i < remain; i++) {
                    population[task].add(front.get(i));
                }
                break;
            }
            idx ++;
        }
    }

    private SolutionSet[] splitPopulation() {
        SolutionSet[] resPopulation = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
            resPopulation[k] = splitPopulation(k);
        }
        return resPopulation;
    }

    private SolutionSet splitPopulation(int task) {
        SolutionSet resPopulation = new SolutionSet();
        for (int i = 0; i < population[task].size(); i++) {
            Solution sol = population[task].get(i);
            int start = problemSet_.get(task).getStartObjPos();
            int end = problemSet_.get(task).getEndObjPos();
            Solution newSolution = new Solution(end - start + 1);
            for (int kk = start; kk <= end; kk++)
                newSolution.setObjective(kk - start, sol.getObjective(kk));
            resPopulation.add(newSolution);
        }
        return resPopulation;
    }

    private void saveIGD(){
        SolutionSet[] resPopulation = splitPopulation();
        for (int k = 0; k < taskNum; k++){
            igds[k].add(indicators[k].getIGD(resPopulation[k]));
        }
    }

    private void saveIGD(int task){
        SolutionSet resPopulation = splitPopulation(task);
        igds[task].add(indicators[task].getIGD(resPopulation));
    }
}
