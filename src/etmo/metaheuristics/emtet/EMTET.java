package etmo.metaheuristics.emtet;

import etmo.core.*;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PORanking;
import etmo.util.PseudoRandom;
import etmo.util.comparators.CrowdingComparator;

public class EMTET extends MtoAlgorithm {


    private int populationSize;

    private SolutionSet[] population;
//    private SolutionSet offspringPopulation;


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

        initPopulation();
        evaluations = 0;

        while (evaluations < maxEvaluations) {
            
            createOffspringPopulation();
            getNextPopulation();
        }

        return population;
    }

    private void createOffspringPopulation() {
    }

    private void getNextPopulation() {
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
                population[i].add(newSolution);
            } // for
            assignFitness(population);
        }


    }

    private void assignFitness(SolutionSet[] population) {
        for (int j = 0; j < population.length; j++){
            for (int i = 0; i < population[j].size(); i++){
                population[j].get(i).setLocation(Integer.MAX_VALUE);
            }
            rankSolutionOnTask(population[j], j);
        }

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
