package etmo.metaheuristics.matde;

import etmo.core.*;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import sun.font.TrueTypeFont;

import java.beans.PropertyChangeSupport;

public class MaTDE extends MtoAlgorithm {
    // Population Size
    private int populationSize;
    private int archiveSize;
    // Population
    private SolutionSet[] population;
    private SolutionSet[] archive;
    private SolutionSet offspringPopulation;

    int evaluations;
    int maxEvaluations;

    Operator crossover;
    Operator mutation;

    double alpha;

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

        probability = new double[problemSet_.size()][problemSet_.size()];
        reward = new double[problemSet_.size()][problemSet_.size()];

        crossover = operators_.get("crossover");
        mutation = operators_.get("mutation");



        evaluations = 0;
        initPopulation();
        while (evaluations < maxEvaluations){
            createOffspringPopulation();

        }
        return population;
    }

    private void initPopulation() throws JMException, ClassNotFoundException {
        population = new SolutionSet[problemSet_.size()];
        archive = new SolutionSet[problemSet_.size()];
        for (int k = 0; k < problemSet_.size(); k++){
            population[k] = new SolutionSet(populationSize);
            archive[k] = new SolutionSet(archiveSize);
            for (int i = 0; i < populationSize; i++){
                Solution newSolution = new Solution(problemSet_);
                problemSet_.get(k).evaluate(newSolution);
                problemSet_.get(k).evaluateConstraints(newSolution);
                population[k].add(newSolution);
                putArchive(k, newSolution);
            }

            for (int kk = 0; kk < problemSet_.size(); kk++){
                probability[k][kk] = 0.0;
                reward[k][kk] = 1.0;
            }
//            TODO 不知道需不需要这个
//            assignFitness(population[k]);
        }
    }

//    private void assignFitness(SolutionSet pop) {
//        for (int i = 0; i < pop.size(); i++)
//            pop.get(i).setLocation(Integer.MAX_VALUE);
//    }


    private void createOffspringPopulation() {
        // 每个任务单独进行
        for (int k = 0; k < problemSet_.size(); k++){
            offspringPopulation = new SolutionSet(populationSize);
            p = PseudoRandom.randDouble();
            if (p > alpha){
                // 不迁移。子种群内部进行交叉变异。
                for (int i = 0; i < populationSize; i++){
                    Solution offSpring;
                    int r1 = i;
                    while (r1 == i)
                        r1 = PseudoRandom.randInt(0, populationSize);

                    // 差分变异。这里按照原算法是不完全差分（仅用了两个父体）。
                    Solution[] parents = new Solution[3];
                    parents[0] = population[k].get(r1);
                    parents[1] = population[k].get(i);
                    parents[2] = population[k].get(i);
                    // TODO 这里调试看一下有没有正常执行
                    offSpring = (Solution) crossover.execute(
                            new Object[] {population[k].get(i), parents});

                    problemSet_.get(k).evaluate(offSpring);
                    // 由于原算法是单目标，这里多目标比较用支配关系。
                    boolean dominated = true;
                    for (int j = 0; j < problemSet_.get(k).getNumberOfObjectives(); j++) {
                        if (offSpring.getObjective(j) > population[k].get(i).getObjective(j)) {
                            dominated = false;
                            break
                        }
                    }
                    // 若子代比父代好，则直接替换掉父代。
                    if (dominated)
                        population[k].replace(i, offSpring);
                }
            }
            else{
                // 迁移。找一个最合适的辅助种群迁移到目标种群。
                int assist = findAssistTask(k);

                for (int i = 0; i < populationSize; i++){
                    Solution offSpring;
                    r1 = PseudoRandom.randInt(0, populationSize);

                    // TODO 好像没有均匀交叉。
                }
            }
        }
    }

    private void findAssistTask(int task){

    }

    private void putArchive(int task, Solution p){
        if (archive[task].size() < archiveSize){
            archive[task].add(p);
        }
        else{
            int replace = PseudoRandom.randInt(0, archiveSize) % archiveSize;
            archive[task].replace(replace, p);
        }
    }
}
