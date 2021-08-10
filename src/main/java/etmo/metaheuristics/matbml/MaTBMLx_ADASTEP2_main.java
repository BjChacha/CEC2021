package etmo.metaheuristics.matbml;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.comparators.LocationComparator;
import etmo.util.logging.LogIGD;

import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class MaTBMLx_ADASTEP2_main {
    // 0: random; 1: WD; 2: KL;
    static int DISTANCE_TPYE = 1;
    static int INIT_SCORE = 3;
    static double BETTER_THRESHOLD = 0.1;
    // 1: ns; 2: dr
    static int ENVIRONMENT_SELECTION_TYPE = 1;
    static double TRANSFER_SCALE = 1;
    static int BASE_RUN_TIME = 5;

    static final String ALGO_NAME = "MaOEAC";
    static final boolean LOG_IGD = true;

    public static void main(String[] args) throws JMException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
        ProblemSet problemSet;
        MtoAlgorithm algorithm;
        Operator crossover;
        HashMap parameters;

        int problemStart = 25;
        int problemEnd = 32;
        int times = 10;

        DecimalFormat form = new DecimalFormat("#.####E0");
        System.out.println("Algo: MaTBML.");

        String benchmark_name;

        if (problemEnd < problemStart)
            problemEnd = problemStart;

        for (int pCase = problemStart; pCase <= problemEnd; pCase++) {
                // CEC2021
                benchmark_name = "CEC2021";
                problemSet = (ProblemSet) Class
                        .forName("etmo.problems.benchmarks_CEC2021.ETMOF" + pCase)
                        .getMethod("getProblem")
                        .invoke(null, null);

                // // WCCI 2020
                // benchmark_name = "WCCI2020";
                // problemSet = (ProblemSet) Class
                //         .forName("etmo.problems.benchmarks_WCCI2020.MATP" + pCase)
                //         .getMethod("getProblem")
                //         .invoke(null, null);


                int taskNum = problemSet.size();

                String[] pf = new String[taskNum];
                for (int k = 0; k < taskNum; k++) {
                    pf[k] = "resources/PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
                }

                String pSName = problemSet.get(0).getName();
                System.out.println(pSName + "\ttaskNum = " + taskNum + "\tfor " + times + " times.");

                algorithm = new MaTBMLx_ADASTEP2(problemSet);
                algorithm.setInputParameter("populationSize", 100);
                algorithm.setInputParameter("maxEvaluations", 1000 * 100 * taskNum);
                algorithm.setInputParameter("transferScale", TRANSFER_SCALE);
                algorithm.setInputParameter("algoName", ALGO_NAME);
                algorithm.setInputParameter("distanceType", DISTANCE_TPYE);
                algorithm.setInputParameter("initScore", INIT_SCORE);
                algorithm.setInputParameter("betterThreshold", BETTER_THRESHOLD);
                algorithm.setInputParameter("environmentSelectionType", ENVIRONMENT_SELECTION_TYPE);
                algorithm.setInputParameter("baseRunTime", BASE_RUN_TIME);

//            // Randomly DE
//            parameters = new HashMap();
//            parameters.put("CR_LB", 0.1);
//            parameters.put("CR_UB", 0.9);
//            parameters.put("F_LB", 0.1);
//            parameters.put("F_UB", 2.0);
//            crossover = CrossoverFactory.getCrossoverOperator("RandomDECrossover",parameters);

//                // Transfer DE
//                parameters = new HashMap();
//                crossover = CrossoverFactory.getCrossoverOperator("TransferDECrossover", parameters);

//                // Random UF
//                parameters = new HashMap();
//                parameters.put("CR_LB", 0.1);
//                parameters.put("CR_UB", 0.9);
//                crossover = CrossoverFactory.getCrossoverOperator("RandomUniformCrossover", parameters);

                // SBX
                parameters = new HashMap();
                parameters.put("probability", 0.9);
                parameters.put("distributionIndex", 20.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                algorithm.addOperator("crossover", crossover);

                double[][] igds = new double[taskNum][times];
                for (int t = 0; t < times; t++) {
                    SolutionSet[] population = algorithm.execute();

                    SolutionSet[] resPopulation = new SolutionSet[taskNum];
                    for (int k = 0; k < taskNum; k++) {
                        resPopulation[k] = new SolutionSet();
                        for (int i = 0; i < population[k].size(); i++) {
                            Solution sol = population[k].get(i);
                            int start = problemSet.get(k).getStartObjPos();
                            int end = problemSet.get(k).getEndObjPos();
                            Solution newSolution = new Solution(end - start + 1);
                            for (int kk = start; kk <= end; kk++)
                                newSolution.setObjective(kk - start, sol.getObjective(kk));
                            resPopulation[k].add(newSolution);
                        }
                    }
                    double igd;
                    for (int k = 0; k < taskNum; k++) {
                        QualityIndicator indicator = new QualityIndicator(problemSet.get(k), pf[k]);
                        if (population[k].size() == 0)
                            continue;

                        igd = indicator.getIGD(resPopulation[k]);
                        igds[k][t] = igd;
                    }
                }

                if (LOG_IGD) {
                    LogIGD.LogIGD("MaTBMLx(MaOEAC_IEMIX_PUNISH_SBX_NoLEADER_WD_NewRA5)" + "_" + benchmark_name, pCase, igds);
                }
            }
    }
}
