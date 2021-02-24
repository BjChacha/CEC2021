package etmo.metaheuristics.matmy2;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.util.JMException;
import etmo.metaheuristics.matmy2.libs.*;
import etmo.util.PseudoRandom;
import etmo.util.logging.LogPopulation;

import java.util.Arrays;
import java.util.HashMap;

public class MaTMY2 extends MtoAlgorithm {
    MaTAlgorithm[] optimizers;
    SolutionSet[] populations;

    private int populationSize;
    private int maxEvaluations;
    private int evaluations;
    int taskNum;

    int[] transferSourceIndexes;
    int transferVolume;

    boolean[] inProcess;

    private String algoName;

    public MaTMY2(ProblemSet problemSet) {
        super(problemSet);
    }

    private void initOptimizers() throws JMException, ClassNotFoundException {
        optimizers = new MaTAlgorithm[taskNum];
        initPopulations();

        if (algoName.equalsIgnoreCase("MOEAD")){
            Operator crossover;
            Operator mutation;

            HashMap parameters;

            for (int k = 0; k < taskNum; k++) {
                ProblemSet pS = problemSet_.getTask(k);
                optimizers[k] = new MOEAD(pS, populations[k]);
                optimizers[k].setInputParameter("populationSize", 100);
                optimizers[k].setInputParameter("maxEvaluations", 1000 * 100);

                optimizers[k].setInputParameter("dataDirectory", "D:\\_r\\EA\\ETMO\\MTO-cec2021-\\resources\\weightVectorFiles\\moead");

                optimizers[k].setInputParameter("T", 20);
                optimizers[k].setInputParameter("delta", 0.9);
                optimizers[k].setInputParameter("nr", 2);

                parameters = new HashMap();
                parameters.put("CR", 1.0);
                parameters.put("F", 0.5);
                crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / pS.get(0).getNumberOfVariables());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                optimizers[k].addOperator("crossover", crossover);
                optimizers[k].addOperator("mutation", mutation);

                ((MOEAD) optimizers[k]).initState();
            }
        }
        else if (algoName.equalsIgnoreCase("MaOEAC")){
            Operator crossover;
            Operator mutation;
            Operator selection;

            HashMap parameters;

            for (int k = 0; k < taskNum; k++) {
                ProblemSet pS = problemSet_.getTask(k);
                optimizers[k] = new MaOEAC(pS, populations[k]);

                optimizers[k].setInputParameter("populationSize",100);
                optimizers[k].setInputParameter("maxGenerations",1000);

                parameters = new HashMap();
                parameters.put("probability", 1.0);
                parameters.put("distributionIndex", 30.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / pS.get(0).getNumberOfVariables());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                parameters = null;
                selection = SelectionFactory.getSelectionOperator("RandomSelection", parameters);

                optimizers[k].addOperator("crossover", crossover);
                optimizers[k].addOperator("mutation", mutation);
                optimizers[k].addOperator("selection", selection);

                ((MaOEAC) optimizers[k]).initState();
            }
        }
        else if (algoName.equalsIgnoreCase("NSGAII")){
            Operator crossover;
            Operator mutation;
            Operator selection;

            HashMap parameters;

            for (int k = 0; k < taskNum; k++){
                ProblemSet pS = problemSet_.getTask(k);
                optimizers[k] = new NSGAII(pS, populations[k]);

                optimizers[k].setInputParameter("populationSize", 100);
                optimizers[k].setInputParameter("maxEvaluations", 100 * 1000);

                parameters = new HashMap();
                parameters.put("probability", 0.9);
                parameters.put("distributionIndex", 20.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / pS.getMaxDimension());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                parameters = null;
                selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters);

                optimizers[k].addOperator("crossover", crossover);
                optimizers[k].addOperator("mutation", mutation);
                optimizers[k].addOperator("selection", selection);

                optimizers[k].initState();
            }
        }
        else{
            throw new JMException("Exception in " + algoName + ": Not such algorithm.");
        }
    }

    private void initPopulations() throws ClassNotFoundException, JMException {
        populations = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++){
            populations[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++){
                Solution newSolution = new Solution(problemSet_.getTask(k));
                problemSet_.get(k).evaluate(newSolution);
                populations[k].add(newSolution);
            }
        }
        evaluations += (taskNum * populationSize);
        LogPopulation.LogPopulation("MaTMY2", populations, problemSet_, evaluations);
    }

    private void initState(){
        taskNum = problemSet_.size();

        populationSize = (Integer) getInputParameter("populationSize");
        maxEvaluations = (Integer) getInputParameter("maxEvaluations");
        transferVolume = (Integer) getInputParameter("transferVolume");

        algoName = (String) getInputParameter("algoName");

        transferSourceIndexes = new int[taskNum];
        inProcess = new boolean[taskNum];
        Arrays.fill(inProcess, true);
    }

    private void Converge(int times) throws JMException, ClassNotFoundException {
        for (int t = 0; t < times; t++) {
            for (int k = 0; k < taskNum; k++) {
                if (inProcess[k]) {
                    inProcess[k] = optimizers[k].step();
                }
            }
            evaluations += (taskNum * populationSize);
            if (evaluations % (20 * taskNum * populationSize) == 0) {
//                System.out.println(totalCount+"x"+taskNum+"x"+populationSize+"="+totalCount * taskNum * populationSize);
                LogPopulation.LogPopulation("MaTMY2", populations, problemSet_, evaluations);
            }
        }
    }

    private void tentativeTransfer() throws JMException, ClassNotFoundException {
        SolutionSet[] subPop = prepareTransferIndividuals();

        int[] srcTask = getTransferSourceIndexes();

        transferIndividuals(subPop, srcTask);

        double[][] bestBefore = new double[taskNum][];
        for (int k = 0; k < taskNum; k++){
            bestBefore[k] = populations[k].getBestObjectiveVector(0);
        }

        Converge(1);

        Arrays.fill(transferSourceIndexes, -1);
        double[][] bestAfter = new double[taskNum][];
        for (int k = 0; k < taskNum; k++){
            bestAfter[k] = populations[k].getBestObjectiveVector(0);

            // 宽松更优
            boolean isBetter = false;
            for (int i = 0; i < problemSet_.get(k).getNumberOfObjectives(); i++) {
                if (bestAfter[k][i] < bestBefore[k][i]) {
                    isBetter = true;
                    break;
                }
            }
//            // 严格更优
//            boolean isBetter = true;
//            for (int i = 0; i < problemSet_.get(k).getNumberOfObjectives(); i++) {
//                if (bestAfter[k][i] > bestBefore[k][i]) {
//                    isBetter = false;
//                    break;
//                }
//            }
            if (isBetter)
                transferSourceIndexes[k] = srcTask[k];
        }
    }

    private void selectiveTransfer(int times) throws JMException, ClassNotFoundException {
        for (int t = 0; t < times; t++) {
            SolutionSet[] subPop = prepareTransferIndividuals();
            transferIndividuals(subPop, transferSourceIndexes);
            Converge(1);
        }
    }

    private int[] getTransferSourceIndexes(){
        int[] srcTask = new int[taskNum];
        for (int k = 0; k < taskNum; k++) {
            srcTask[k] = k;
            while (srcTask[k] == k)
                srcTask[k] = PseudoRandom.randInt(0, taskNum - 1);
        }
        return srcTask;
    }

    private SolutionSet[] prepareTransferIndividuals() {
        SolutionSet[] subPop = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++){
            int[] permutation = new int[transferVolume];
            Utils.randomPermutation(permutation, transferVolume, populations[k].size());

            subPop[k] = new SolutionSet(transferVolume);
            for (int i = 0; i < transferVolume; i++){
                Solution indiv = new Solution(populations[k].get(permutation[i]));
                subPop[k].add(indiv);
            }
        }
        return subPop;
    }

    private void transferIndividuals(SolutionSet[] subPop, int[] srcTaskList){
        for (int k = 0; k < taskNum; k++){
            if (srcTaskList[k] >= 0) {
                int[] permutation = new int[transferVolume];
                Utils.randomPermutation(permutation, transferVolume, populations[k].size());
                for (int i = 0; i < transferVolume; i++) {
                    Solution newIndiv = subPop[srcTaskList[k]].get(i);
                    newIndiv.setProblemSet_(problemSet_.getTask(k));
                    populations[k].replace(permutation[i], newIndiv);
                }
            }
        }
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        evaluations = 0;
        initState();
        initOptimizers();
//        Converge(1);
        while (evaluations < maxEvaluations) {
            tentativeTransfer();
            selectiveTransfer(9);
        }
        return populations;
    }
}