package etmo.metaheuristics.matmy;

import etmo.core.*;
import etmo.util.*;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.DominanceComparator;
import etmo.util.logging.LogPopulation;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class MaTMY extends MtoAlgorithm {
    // Population Size per task
    private int populationSize;
    // Archive Size per task
    private int archiveSize;
    // Population
    private SolutionSet[] population;
    private SolutionSet[] archives;

    private int evaluations;
    private int maxEvaluations;

    Operator crossover1;
    Operator crossover2;
    Operator mutation;
    Operator selector;
    Comparator dominance;

    double alpha;
    double tentative_alpha;
    double selective_alpha;

    double ro;
    double shrinkRate;
    double replaceRate;
    double transferThreshold;

    double[][] probability;
    double[][] reward;
    // 计算收敛幅度
    double[] omega;

    int[] transferSource;

    /**
     * Constructor
     *
     * @param problemSet The problem to be solved
     */
    public MaTMY(ProblemSet problemSet) {
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initialState();

        initPopulation();              // evaluations: T x N x 1
        LogPopulation.LogPopulation("MaTMY",population, problemSet_,evaluations);

        while (evaluations < maxEvaluations){
            for (int t = 0; t < 10; t++) {
                soloConverge();            // evaluations: T x N x times
                if (evaluations % (problemSet_.size() * populationSize * 20) == 0)
                    LogPopulation.LogPopulation("MaTMY",population, problemSet_,evaluations);
            }

            tentativeTransfer();           // evaluations: T x N x 1
            if (evaluations % (problemSet_.size() * populationSize * 20) == 0)
                LogPopulation.LogPopulation("MaTMY",population, problemSet_,evaluations);

            for (int t = 0; t < 10; t++) {
                selectiveTransfer();       // evaluations: T x N x times
                if (evaluations % (problemSet_.size() * populationSize * 20) == 0)
                    LogPopulation.LogPopulation("MaTMY",population, problemSet_,evaluations);
            }
//            if (evaluations % (problemSet_.size() * populationSize * 20) == 0)
//                LogPopulation.LogPopulation("MaTMY",population, problemSet_,evaluations);
        }
        return population;
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

    // 收敛阶段
    private void soloConverge() throws JMException {
        for (int k = 0; k < problemSet_.size(); k++) {
            List<Solution> offspringList = new ArrayList<Solution>();
            for (int i = 0; i < populationSize; i++) {
                Solution offSpring;
                int r1, r2;
                r1 = r2 = i;
                while (r1 == i || r2 == i || r1 == r2) {
                    r1 = PseudoRandom.randInt(0, populationSize - 1);
                    r2 = PseudoRandom.randInt(0, populationSize - 1);
                }
                // 差分交叉
                Solution[] parents = new Solution[3];
                parents[0] = new Solution(population[k].get(r1));
                parents[1] = new Solution(population[k].get(r2));
                parents[2] = new Solution(population[k].get(i));
                offSpring = (Solution) crossover1.execute(new Object[]{population[k].get(i), parents});
                // 变异
                mutation.execute(offSpring);

                problemSet_.get(k).evaluate(offSpring);
                evaluations++;

                int flag = dominance.compare(offSpring, population[k].get(i));
                if (flag == -1)
                    population[k].replace(i, offSpring);
                else if (flag == 0)
                    offspringList.add(offSpring);
            }
            // 更新种群
            EliteSelect(k, offspringList);
        }
    }

    private void tentativeTransfer() throws JMException {
        updateArchives();
        for (int k = 0; k < problemSet_.size(); k++) {
            List<Solution> offspringList = new ArrayList<Solution>();

            int assistK = findAssistTask(k);

            double[] pBest = new double[problemSet_.get(k).getNumberOfObjectives()];
            for (int j = 0; j <= problemSet_.get(k).getEndObjPos() - problemSet_.get(k).getStartObjPos(); j++) {
                pBest[j] = Double.MAX_VALUE;
                for (int i = 0; i < populationSize; i++)
                    pBest[j] = Math.min(pBest[j], population[k].get(i).getObjective(problemSet_.get(k).getStartObjPos() + j));
            }

            for (int i = 0; i < populationSize; i++) {
                int r1 = PseudoRandom.randInt(0, populationSize - 1);
                Solution[] parents = new Solution[2];
                parents[0] = new Solution(population[assistK].get(r1));
                parents[1] = new Solution((Solution) selector.execute(population[k]));
                Solution[] offSprings = (Solution[]) crossover2.execute(parents);
                Solution offSpring = offSprings[PseudoRandom.randInt(0, 1)];

                problemSet_.get(k).evaluate(offSpring);
                evaluations++;

                int flag = dominance.compare(offSpring, population[k].get(i));
                if (flag == -1)
                    population[k].replace(i, offSpring);
                else if (flag == 0)
                    offspringList.add(offSpring);
            }
            EliteSelect(k, offspringList);

            double[] pBestAfter = new double[problemSet_.get(k).getNumberOfObjectives()];
            for (int j = 0; j <= problemSet_.get(k).getEndObjPos() - problemSet_.get(k).getStartObjPos(); j++) {
                pBestAfter[j] = Double.MAX_VALUE;
                for (int i = 0; i < populationSize; i++)
                    pBestAfter[j] = Math.min(pBestAfter[j], population[k].get(i).getObjective(problemSet_.get(k).getStartObjPos() + j));
            }

            boolean isBetter = false;
            for (int i = 0; i < problemSet_.get(k).getNumberOfObjectives(); i++) {
                if (pBestAfter[i] < pBest[i]) {
                    isBetter = true;
                    break;
                }
            }

            if (isBetter) {
                reward[k][assistK] /= shrinkRate;
                transferSource[k] = assistK;
            }
            else
                reward[k][assistK] = 0;

        }
    }

    private void selectiveTransfer() throws JMException {
        for (int k = 0; k < problemSet_.size(); k++){
            if (transferSource[k] < 0)
                continue;

            List<Solution> offspringList = new ArrayList<Solution>();

            int assistK = transferSource[k];

            double[] pBest = new double[problemSet_.get(k).getNumberOfObjectives()];
            for (int j = 0; j <= problemSet_.get(k).getEndObjPos() - problemSet_.get(k).getStartObjPos(); j++) {
                pBest[j] = Double.MAX_VALUE;
                for (int i = 0; i < populationSize; i++)
                    pBest[j] = Math.min(pBest[j], population[k].get(i).getObjective(problemSet_.get(k).getStartObjPos() + j));
            }

            for (int i = 0; i < populationSize; i++) {
                int r1 = PseudoRandom.randInt(0, populationSize - 1);
                Solution[] parents = new Solution[2];
                parents[0] = new Solution(population[assistK].get(r1));
                parents[1] = new Solution((Solution) selector.execute(population[k]));
                Solution[] offSprings = (Solution[]) crossover2.execute(parents);
                Solution offSpring = offSprings[PseudoRandom.randInt(0, 1)];

                problemSet_.get(k).evaluate(offSpring);
                evaluations++;

                int flag = dominance.compare(offSpring, population[k].get(i));
                if (flag == -1)
                    population[k].replace(i, offSpring);
                else if (flag == 0) {
                    offspringList.add(offSpring);
                }
            }

            EliteSelect(k, offspringList);

            double[] pBestAfter = new double[problemSet_.get(k).getNumberOfObjectives()];
            for (int j = 0; j <= problemSet_.get(k).getEndObjPos() - problemSet_.get(k).getStartObjPos(); j++) {
                pBestAfter[j] = Double.MAX_VALUE;
                for (int i = 0; i < populationSize; i++)
                    pBestAfter[j] = Math.min(pBestAfter[j], population[k].get(i).getObjective(problemSet_.get(k).getStartObjPos() + j));
            }

            boolean isBetter = false;
            for (int i = 0; i < problemSet_.get(k).getNumberOfObjectives(); i++) {
                if (pBestAfter[i] < pBest[i]) {
                    isBetter = true;
                    break;
                }
            }

            if (isBetter)
                reward[k][assistK] /= shrinkRate;
            else
                reward[k][assistK] *= shrinkRate;

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


    private void updateArchives() {
        for (int k = 0; k < problemSet_.size(); k++) {
            for (int i = 0; i < populationSize; i++){
                if (PseudoRandom.randDouble() < replaceRate)
                    putArchive(k, population[k].get(i));
            }
        }
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

    private void EliteSelect(int task, List<Solution> offspringList){
        SolutionSet offspringPopulation = new SolutionSet(offspringList);
        SolutionSet union = population[task].union(offspringPopulation);

        int remain = populationSize;
        int idx = 0;
        boolean[] chosen = new boolean[problemSet_.getTotalNumberOfObjs()];
        for (int i = problemSet_.get(task).getStartObjPos(); i <= problemSet_.get(task).getEndObjPos(); i++)
            chosen[i] = true;

        Ranking ranking = new Ranking(union);
        Distance distance = new Distance();
        population[task].clear();
        SolutionSet front = null;

        while ((remain > 0)) {
            front = ranking.getSubfront(idx);
            if (remain >= front.size()) {
                distance.crowdingDistanceAssignment(front, problemSet_.get(task).getNumberOfObjectives(), chosen);
                for (int i = 0; i < front.size(); i++)
                    population[task].add(front.get(i));
                remain -= front.size();
            } else {
                distance.crowdingDistanceAssignment(front, problemSet_.get(task).getNumberOfObjectives(), chosen);
                front.sort(new CrowdingComparator());
                for (int i = 0; i < remain; i++)
                    population[task].add(front.get(i));
                break;
            }
            idx++;
        }
    }

    private void initialState(){
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

        omega = new double[problemSet_.size()];
        transferSource = new int[problemSet_.size()];

        crossover1 = operators_.get("crossover1");
        crossover2 = operators_.get("crossover2");
        mutation = operators_.get("mutation");
        selector = operators_.get("operator");
        dominance = new DominanceComparator();

        evaluations = 0;
    }
}
