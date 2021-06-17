package etmo.metaheuristics.mteaad;

import etmo.core.*;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.Ranking;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.PdfComparator;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.Covariance;

import java.util.Arrays;

public class MTEA_AD extends Algorithm {
    public MTEA_AD(ProblemSet problemSet) {
        super(problemSet);
        //	System.out.println("sup: " + problemSet_.get(0).getHType());
    } // NSGAII

    public SolutionSet execute() throws JMException, ClassNotFoundException {
        int populationSize;
        int maxEvaluations;
        int evaluations;
        int taskNum;


        SolutionSet population;
        SolutionSet offspringPopulation;
        SolutionSet union;

        Operator mutationOperator;
        Operator crossoverOperator;
        Operator selectionOperator;

        Distance distance = new Distance();

        // Read the parameters
        populationSize = ((Integer) getInputParameter("populationSize")).intValue();
        maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();

        // Initialize the variables
        population = new SolutionSet(populationSize);
        taskNum = problemSet_.size();
        evaluations = 0;


        // Read the operators
        mutationOperator = operators_.get("mutation");
        crossoverOperator = operators_.get("crossover");
        selectionOperator = operators_.get("selection");

        // Create the initial solutionSet
        Solution newSolution;

        for (int i = 0; i < populationSize; i++) {
                newSolution = new Solution(problemSet_);
                problemSet_.get(0).evaluate(newSolution);
                problemSet_.get(0).evaluateConstraints(newSolution);
                evaluations++;
                population.add(newSolution);
        } // for



        // Generations
        while (evaluations < maxEvaluations) {

            // Create the offSpring solutionSet
            offspringPopulation = new SolutionSet(populationSize);
            Solution[] parents = new Solution[2];
            for (int i = 0; i < (populationSize / 2); i++) {
                if (evaluations < maxEvaluations) {
                    // obtain parents
                    parents[0] = (Solution) selectionOperator.execute(population);
                    parents[1] = (Solution) selectionOperator.execute(population);
                    Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
                    mutationOperator.execute(offSpring[0]);
                    mutationOperator.execute(offSpring[1]);
                    problemSet_.get(0).evaluate(offSpring[0]);
                    problemSet_.get(0).evaluateConstraints(offSpring[0]);
                    problemSet_.get(0).evaluate(offSpring[1]);
                    problemSet_.get(0).evaluateConstraints(offSpring[1]);
                    offspringPopulation.add(offSpring[0]);
                    offspringPopulation.add(offSpring[1]);
                    evaluations += 2;
                } // if
            } // for

            // Create the solutionSet union of solutionSet and offSpring
            union = ((SolutionSet) population).union(offspringPopulation);

            // Ranking the union
            Ranking ranking = new Ranking(union);

            int remain = populationSize;
            int index = 0;
            SolutionSet front = null;
            population.clear();

            // Obtain the next front
            front = ranking.getSubfront(index);

            while ((remain > 0) && (remain >= front.size())) {
                // Assign crowding distance to individuals
                distance.crowdingDistanceAssignment(front, problemSet_.get(0).getNumberOfObjectives());
                // Add the individuals of this front
                for (int k = 0; k < front.size(); k++) {
                    population.add(front.get(k));
                } // for

                // Decrement remain
                remain = remain - front.size();

                // Obtain the next front
                index++;
                if (remain > 0) {
                    front = ranking.getSubfront(index);
                } // if
            } // while

            // Remain is less than front(index).size, insert only the best one
            if (remain > 0) { // front contains individuals to insert
                distance.crowdingDistanceAssignment(front, problemSet_.get(0).getNumberOfObjectives());
                front.sort(new CrowdingComparator());
                for (int k = 0; k < remain; k++) {
                    population.add(front.get(k));
                } // for

                remain = 0;
            } // if

        } // while

        Ranking ranking = new Ranking(population);

        return ranking.getSubfront(0);

    } // execute

    public SolutionSet[] executeMultiTask() throws JMException, ClassNotFoundException {
        int populationSize;
        int maxEvaluations;
        int evaluations;
        int taskNum;
        double[] epsilon;


        SolutionSet[] population;
        SolutionSet offspringPopulation;
        SolutionSet union;

        Operator mutationOperator;
        Operator crossoverOperator;
        Operator selectionOperator;

        Distance distance = new Distance();

        // Read the parameters
        populationSize = ((Integer) getInputParameter("populationSize")).intValue();
        maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();

        // Initialize the variables
        taskNum = problemSet_.size();
        evaluations = 0;
        population = new SolutionSet[taskNum];
        epsilon = new double[taskNum];
        Arrays.fill(epsilon,0.2);




        // Read the operators
        mutationOperator = operators_.get("mutation");
        crossoverOperator = operators_.get("crossover");
        selectionOperator = operators_.get("selection");

        // Create the initial solutionSet
        Solution newSolution;

        for (int t = 0; t < taskNum; t++){
            population[t] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++) {
                newSolution = new Solution(problemSet_);
                problemSet_.get(t).evaluate(newSolution);
                problemSet_.get(t).evaluateConstraints(newSolution);
                evaluations++;
                population[t].add(newSolution);
            } // for
        }


        // Generations
        while (evaluations < maxEvaluations) {

            for (int t = 0; t < taskNum; t++){

                // Create the offSpring solutionSet
                offspringPopulation = new SolutionSet(populationSize);
                Solution[] parents = new Solution[2];
                for (int i = 0; i < (populationSize / 2); i++) {
                    if (evaluations < maxEvaluations) {
                        // obtain parents
                        parents[0] = (Solution) selectionOperator.execute(population[t]);
                        parents[1] = (Solution) selectionOperator.execute(population[t]);
                        Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
                        mutationOperator.execute(offSpring[0]);
                        mutationOperator.execute(offSpring[1]);
                        problemSet_.get(t).evaluate(offSpring[0]);
                        problemSet_.get(t).evaluateConstraints(offSpring[0]);
                        problemSet_.get(t).evaluate(offSpring[1]);
                        problemSet_.get(t).evaluateConstraints(offSpring[1]);
                        offspringPopulation.add(offSpring[0]);
                        offspringPopulation.add(offSpring[1]);
                        evaluations += 2;
                    } // if
                } // for


                double transRand = PseudoRandom.randDouble();
                if (transRand <= 0.1){
                    // init ad model
                    int vNum = problemSet_.get(t).getNumberOfVariables();
                    double[][] v = new double[populationSize][vNum];
                    for (int i = 0; i < populationSize; i++){
                        for (int j = 0; j < vNum; j++){
                            v[i][j] = population[t].get(i).getDecisionVariables()[j].getValue() + PseudoRandom.randDouble(-1e-5, 1e-5);
                        }
                    }
                    RealMatrix mx = MatrixUtils.createRealMatrix(v);
                    RealMatrix cov = new Covariance(mx).getCovarianceMatrix();
                    double[] vMean = new double[vNum];
                    for (int i = 0; i < vNum; i++){
                        vMean[i] = StatUtils.mean(mx.getColumn(i));
                    }

                    // how to debug singular?
                    MultivariateNormalDistribution adModeal = null;
                    try {
                        adModeal = new MultivariateNormalDistribution(vMean, cov.getData());
                    } catch (SingularMatrixException e) {
                        e.printStackTrace();
                    } catch (DimensionMismatchException e) {
                        e.printStackTrace();
                    } catch (NonPositiveDefiniteMatrixException e) {
                        e.printStackTrace();
                    }

                    // generate tranfer solutions
                    SolutionSet transfer = new SolutionSet();
                    for (int i = 0; i < taskNum; i++){
                        if (i == t) continue;
                        for (int j = 0; j < populationSize; j++){
                            int maxD = adModeal.getDimension();
                            int num = problemSet_.get(i).getNumberOfVariables();
                            double[] indiV;
                            if (maxD <= num){
                                indiV = new double[maxD];
                                for (int k = 0; k < maxD; k++){
                                    indiV[k] = population[i].get(j).getDecisionVariables()[k].getValue();
                                }
                            }
                            else {
                                indiV = new double[maxD];
                                for (int k = 0; k < maxD; k++){
                                    if (k >= num) indiV[k] = PseudoRandom.randDouble();
                                    else indiV[k] = population[i].get(j).getDecisionVariables()[k].getValue();
                                }
                            }
                            Solution tmp = new Solution(population[i].get(j));
                            tmp.setPdf(adModeal.density(indiV));
                            transfer.add(tmp);
                        }
                    }
                    transfer.sort(new PdfComparator());
                    SolutionSet transferRes = new SolutionSet();
                    double transNum = epsilon[t] * populationSize * (taskNum - 1);
                    for (int i = 0; i < (int)transNum; i++){
//                        System.out.println(transNum + " : " + epsilon[t]);
                        Solution so = new Solution(transfer.get(i));
                        problemSet_.get(t).evaluate(so);
                        so.setIsTran(1);
                        evaluations++;
                        transferRes.add(so);
                    }

                    union = population[t].union(offspringPopulation).union(transferRes);

                    // generate new pop
                    {
                        // Ranking the union
                        Ranking ranking = new Ranking(union);

                        int remain = populationSize;
                        int index = 0;
                        SolutionSet front = null;
                        population[t].clear();

                        // Obtain the next front
                        front = ranking.getSubfront(index);

                        while ((remain > 0) && (remain >= front.size())) {
                            // Assign crowding distance to individuals
                            distance.crowdingDistanceAssignment(front, problemSet_.getTotalNumberOfObjs(), TaskObjSelection(t));
                            // Add the individuals of this front
                            for (int k = 0; k < front.size(); k++) {
                                population[t].add(front.get(k));
                            } // for

                            // Decrement remain
                            remain = remain - front.size();

                            // Obtain the next front
                            index++;
                            if (remain > 0) {
                                front = ranking.getSubfront(index);
                            } // if
                        } // while

                        // Remain is less than front(index).size, insert only the best one
                        if (remain > 0) { // front contains individuals to insert
                            distance.crowdingDistanceAssignment(front, problemSet_.getTotalNumberOfObjs(), TaskObjSelection(t));
                            front.sort(new CrowdingComparator());
                            for (int k = 0; k < remain; k++) {
                                population[t].add(front.get(k));
                            } // for

                            remain = 0;
                        } // if
                    }

                    // update epsilon
//                    double successNum = 0;
//                    for (int i = 0; i < populationSize; i++){
//                        if (population[t].get(i).getIsTran() == 1) successNum += 1;
//                    }
////                    if (transNum != 0)  epsilon[t] = successNum / transNum;
//                    if (transNum != 0)  epsilon[t] = successNum / populationSize;
//                    else epsilon[t] = 0;
//                    System.out.println("successNum: " + successNum + " transNum: " + transNum);
//                    System.out.println("epsilon " + t + ":" + epsilon[t]);

                    for (int i = 0; i < populationSize; i++){
                        population[t].get(i).setIsTran(0);
                    }

                }// if a <= 0.1

                else {
                    union = population[t].union(offspringPopulation);
                    // Ranking the union
                    Ranking ranking = new Ranking(union);

                    int remain = populationSize;
                    int index = 0;
                    SolutionSet front = null;
                    population[t].clear();

                    // Obtain the next front
                    front = ranking.getSubfront(index);

                    while ((remain > 0) && (remain >= front.size())) {
                        // Assign crowding distance to individuals
                        distance.crowdingDistanceAssignment(front, problemSet_.getTotalNumberOfObjs(), TaskObjSelection(t));
                        // Add the individuals of this front
                        for (int k = 0; k < front.size(); k++) {
                            population[t].add(front.get(k));
                        } // for

                        // Decrement remain
                        remain = remain - front.size();

                        // Obtain the next front
                        index++;
                        if (remain > 0) {
                            front = ranking.getSubfront(index);
                        } // if
                    } // while

                    // Remain is less than front(index).size, insert only the best one
                    if (remain > 0) { // front contains individuals to insert
                        distance.crowdingDistanceAssignment(front, problemSet_.getTotalNumberOfObjs(), TaskObjSelection(t));
                        front.sort(new CrowdingComparator());
                        for (int k = 0; k < remain; k++) {
                            population[t].add(front.get(k));
                        } // for

                        remain = 0;
                    } // if


                }

            } // for each task



        } // while

        for (int i = 0; i < taskNum; i++){
            Ranking ranking = new Ranking(population[i]);
            population[i] = ranking.getSubfront(0);
        }

        return population;


    } // executeMutiTask

    public boolean[] TaskObjSelection(int taskId){
        int start = problemSet_.get(taskId).getStartObjPos();
        int end = problemSet_.get(taskId).getEndObjPos();

        boolean selec[] = new boolean[problemSet_.getTotalNumberOfObjs()];

        for (int i = 0; i < selec.length; i++) {
            if (i < start || i > end)
                selec[i] = false;
            else
                selec[i] = true;
        }
        return selec;
    }
}
