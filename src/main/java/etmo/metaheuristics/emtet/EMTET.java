package etmo.metaheuristics.emtet;

import etmo.core.*;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PORanking;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.LocationComparator;

import java.util.Map;
import java.util.TreeMap;

public class EMTET extends MtoAlgorithm {


    private int populationSize;

    private SolutionSet[] population;
    private SolutionSet offspringPopulation;
//    private SolutionSet transferPopulation;
//    private SolutionSet union;

    int evaluations;
    int maxEvaluations;

    Operator crossover;
    Operator mutation;
    Operator selection;

    int G;

    Distance distance = new Distance();

    /**
     * Constructor
     *
     * @param problemSet The problem to be solved
     */
    public EMTET(ProblemSet problemSet) {
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        populationSize = ((Integer) getInputParameter("populationSize")).intValue();
        maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();
        G = ((Integer) getInputParameter("transferNum")).intValue();

        crossover = operators_.get("crossover");
        mutation = operators_.get("mutation");
        selection = operators_.get("selection");

        evaluations = 0;
        initPopulation();

        while (evaluations < maxEvaluations) {
            createTransferPopulation();
            createOffspringPopulation();
            getNextPopulation();
        }

        return population;
    }

//    未进行迁移解微调：* λ ∼ U(0, 2), p = 0.5,
    private void createTransferPopulation() throws JMException {
        for (int i = 0; i < population.length; i++){
//            存活个体即迁移解
            SolutionSet survivalPopulation = new SolutionSet(G);

//            目标任务中标记为1且位于0层的个体为存活个体
            for (int j = 0; j < G && population[i].get(j).getRank() == 0; j++){
                Solution survival = population[i].get(j);
                if(survival.getSkillFactor() == 1){
                    survivalPopulation.add(survival);
                }
            }

//            目标任务标记清零
            for (int j = 0; j < populationSize; j++){
                population[i].get(j).setSkillFactor(0);
            }

            for (int j = 0; j < population.length; j++){
                if (i == j) continue;
                //            若无存活个体则直接从源任务里面选排名前G个并作标记
                SolutionSet rankHignG = new SolutionSet(G);
                if (survivalPopulation.size() == 0){
                    for (int k = 0; k < G; k++){
                        population[j].get(k).setSkillFactor(1);
                        problemSet_.get(i).evaluate(population[j].get(k));
                        evaluations++;
                        rankHignG.add(population[j].get(k));
                    }
                    population[i] = population[i].union(rankHignG);
                }
//                存在存活个体，转换到源任务j寻找最近的邻居作为新的迁移解
                else {
                    int cycle = G / survivalPopulation.size();
                    int more = G % survivalPopulation.size();
                    for (int k = 0; k < survivalPopulation.size(); k++){
//                        线性变换不影响直接在统一域内计算距离
                           Variable[] x = survivalPopulation.get(k).getDecisionVariables();
                           int num = cycle;
                           if (k < more) num += 1;
                           population[i] = population[i].union(transferNeighbor(x,j,i,num));
                    }

                }

            }
        }
    }

    private SolutionSet transferNeighbor(Variable[] x, int source, int target, int num) throws JMException {
//        记得要做标记
        SolutionSet closestNeighbor = new SolutionSet(num);
        Map<Double,Integer> dis = new TreeMap<Double,Integer>();
        for (int i = 0; i < populationSize; i++){
            Variable[] y = population[source].get(i).getDecisionVariables();
            double sum = 0;
            for(int j = 0 ; j < x.length; j++){
                sum += Math.pow(x[j].getValue() - y[j].getValue(),2);
            }
            dis.put(Math.sqrt(sum), i);
        }
        int cnt = 0;
        for (double key : dis.keySet()){
            if (cnt > num)  break;
            Solution t = population[source].get(dis.get(key));
            t.setSkillFactor(1);
            problemSet_.get(target).evaluate(t);
            evaluations++;
            closestNeighbor.add(t);
            cnt++;
            if (cnt >= closestNeighbor.size())
                break;
        }
        return closestNeighbor;
    }


    private void createOffspringPopulation() throws JMException {
        for (int i = 0; i < population.length; i++){
            int numOffspring = populationSize * 2 - population[i].size();
            offspringPopulation = new SolutionSet(numOffspring);
            Solution[] parents = new Solution[2];
            for (int j = 0; j < numOffspring; j++){
                parents[0] = (Solution) selection.execute(population[i]);
                parents[1] = (Solution) selection.execute(population[i]);

                Solution[] offspring = (Solution[]) crossover.execute(parents);

                mutation.execute(offspring[0]);

                offspring[0].setSkillFactor(0);
                resetObjectives(offspring[0]);
                problemSet_.get(i).evaluate(offspring[0]);
                problemSet_.get(i).evaluateConstraints(offspring[0]);

                evaluations += 1;
                offspringPopulation.add(offspring[0]);
            }
            population[i] = population[i].union(offspringPopulation);
        }

    }

    void resetObjectives(Solution sol) {
        for (int i = 0; i < sol.getNumberOfObjectives(); i++)
            sol.setObjective(i, Double.POSITIVE_INFINITY);
    }

    private void getNextPopulation() {
        for (int i = 0; i < population.length; i++){
            assignFitness(population[i]);
            population[i].sort(new LocationComparator());
            SolutionSet nextPop = new SolutionSet();
            nextPop = nextPop.union(population[i]);
            population[i].clear();
            for (int j = 0; j < populationSize; j++)
                population[i].add(nextPop.get(j));
        }
    }

    private void initPopulation() throws ClassNotFoundException, JMException {

        population = new SolutionSet[problemSet_.size()];
        for (int i = 0; i < problemSet_.size(); i++) {
            population[i] = new SolutionSet(populationSize);
            for (int j = 0; j < populationSize; j++) {
                Solution newSolution = new Solution(problemSet_);
                problemSet_.get(i).evaluate(newSolution);
                problemSet_.get(i).evaluateConstraints(newSolution);
                evaluations++;
//                这里的skillFactor为迁移存活标志
                newSolution.setSkillFactor(0);
                population[i].add(newSolution);
            } // for
            assignFitness(population[i]);
        }

    }

    void assignFitness(SolutionSet pop) {
        for (int i = 0; i < pop.size(); i++)
            pop.get(i).setLocation(Integer.MAX_VALUE);
        for (int i = 0; i < problemSet_.size(); i++)
            rankSolutionOnTask(pop, i);
    }

    private void rankSolutionOnTask(SolutionSet pop, int taskId) {
        int start = problemSet_.get(taskId).getStartObjPos();
        int end = problemSet_.get(taskId).getEndObjPos();

        boolean selec[] = new boolean[problemSet_.getTotalNumberOfObjs()];

        for (int i = 0; i < selec.length; i++) {
            if (i < start || i > end)
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
}
