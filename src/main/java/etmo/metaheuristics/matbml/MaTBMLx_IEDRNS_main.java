package etmo.metaheuristics.matbml;

import etmo.core.*;
import etmo.metaheuristics.matbml.MaTBML;
import etmo.metaheuristics.matbml.MaTBMLx;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.comparators.LocationComparator;
import etmo.util.logging.LogIGD;

import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;

public class MaTBMLx_IEDRNS_main {
    static int SOLELY_CONVERGE_TIMES = 5;
    static int TRANSFER_CONVERGE_TIMES = 5;
    static int IMPLICIT_TRANSFER_NUM = 10;
    static int EXPLICIT_TRANSFER_NUM = 10;

    static boolean IS_IMPLICIT = true;
    static boolean IS_EXPLICIT = true;
    static boolean IS_NS = true;

    static final String ALGO_NAME = "MaOEAC";
    static final boolean LOG_IGD = false;

    public static void main(String[] args) throws JMException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
        ProblemSet problemSet;
        MtoAlgorithm algorithm;
        Operator crossover;
        Operator mutation;
        Operator selection;
        HashMap parameters;

        int problemStart = 25;
        int problemEnd = 32;
        int times = 10;

        DecimalFormat form = new DecimalFormat("#.####E0");
        System.out.println("Algo: MaTBML.");

        String benchmark_name;

        if (problemEnd < problemStart)
            problemEnd = problemStart;

        for (int flag_ns = 1; flag_ns < 2; flag_ns++) {
            for (int flag_explicit = 0; flag_explicit < 2; flag_explicit++) {
                for (int flag_implicit = 0; flag_implicit < 2; flag_implicit++) {
                    if (flag_implicit + flag_explicit <= 0)
                        continue;
                    IS_IMPLICIT = flag_implicit != 0;
                    IS_EXPLICIT = flag_explicit != 0;
                    IS_NS = flag_ns != 0;

                    for (int pCase = problemStart; pCase <= problemEnd; pCase++) {
                        // CEC2021
                        benchmark_name = "CEC2021";
                        problemSet = (ProblemSet) Class
                                .forName("etmo.problems.benchmarks_ETMO.ETMOF" + pCase)
                                .getMethod("getProblem")
                                .invoke(null, null);

                        // // WCCI 2020
                        // benchmark_name = "WCCI2020";
                        // problemSet = (ProblemSet) Class
                        //         .forName("etmo.problems.benchmarks_WCCI2020.MATP" + pCase)
                        //         .getMethod("getProblem")
                        //         .invoke(null, null);


                        int taskNum = problemSet.size();
                        double[] ave = new double[taskNum];

                        String[] pf = new String[taskNum];
                        for (int k = 0; k < taskNum; k++) {
                            pf[k] = "PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
                        }

                        String pSName = problemSet.get(0).getName();
                        System.out.println(pSName + "\ttaskNum = " + taskNum + "\tfor " + times + " times.");

                        algorithm = new MaTBMLx(problemSet);
                        algorithm.setInputParameter("populationSize", 100);
                        algorithm.setInputParameter("maxEvaluations", 1000 * 100 * taskNum);
                        algorithm.setInputParameter("solelyConvergeTimes", SOLELY_CONVERGE_TIMES);
                        algorithm.setInputParameter("transferConvergeTimes", TRANSFER_CONVERGE_TIMES);
                        algorithm.setInputParameter("implicitTransferNum", IMPLICIT_TRANSFER_NUM);
                        algorithm.setInputParameter("explicitTransferNum", EXPLICIT_TRANSFER_NUM);
                        algorithm.setInputParameter("algoName", ALGO_NAME);
                        algorithm.setInputParameter("isImplicit", IS_IMPLICIT);
                        algorithm.setInputParameter("isExplicit", IS_EXPLICIT);
                        algorithm.setInputParameter("isNS", IS_NS);
                        algorithm.setInputParameter("distanceType", 1);

//            // Randomly DE
//            parameters = new HashMap();
//            parameters.put("CR_LB", 0.1);
//            parameters.put("CR_UB", 0.9);
//            parameters.put("F_LB", 0.1);
//            parameters.put("F_UB", 2.0);
//            crossover = CrossoverFactory.getCrossoverOperator("RandomDECrossover",parameters);

                        // SBX
                        parameters = new HashMap();
                        parameters.put("probability", 0.9);
                        parameters.put("distributionIndex", 20.0);
                        crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
                        // Mutation operator
                        parameters = new HashMap();
                        parameters.put("probability", 1.0 / problemSet.getMaxDimension());
                        parameters.put("distributionIndex", 20.0);
                        mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                        // Selection Operator
                        parameters = new HashMap();
                        parameters.put("comparator", new LocationComparator());
                        selection = SelectionFactory.getSelectionOperator("BinaryTournament",
                                parameters);

                        algorithm.addOperator("crossover", crossover);
                        algorithm.addOperator("mutation", mutation);
                        algorithm.addOperator("selection", selection);

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
                        String signImplicit = IS_IMPLICIT?"I":"";
                        String signExplicit = IS_EXPLICIT?"E":"";
                        String signNS = IS_NS?"NS":"DR";
                        LogIGD.LogIGD("MaTBMLx("+signImplicit+signExplicit+signNS + ")_" + benchmark_name, pCase, igds);
            //            for(int i=0;i<taskNum;i++) {
            //                double[] tmp = new double[times];
            //                for (int t = 0; t < times; t++){
            //                    tmp[t] = igds[i][t];
            //                }
            //                double mean = 0;
            //                mean = Arrays.stream(tmp).sum() / times;
            //                System.out.println(form.format(mean));
            //            }
            //            System.out.println();
                    }
                }
            }
        }
    }
}