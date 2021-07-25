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
import java.util.Arrays;
import java.util.HashMap;

public class MaTBMLI_main {
    public static void main(String[] args) throws JMException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
        int K1 = 5;
        int K2 = 5;
        String ALGO_NAME = "MaOEAC";

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
        System.out.println("Algo: MaTBML_exp.");

        String benchmark_name;

        if (problemEnd < problemStart)
            problemEnd = problemStart;
        for (int pCase = problemStart; pCase <= problemEnd; pCase++){
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
            for (int k = 0; k < taskNum; k++){
                pf[k] = "PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
            }

            String pSName = problemSet.get(0).getName();
            System.out.println(pSName + "\ttaskNum = "+taskNum+"\tfor "+times+" times.");

            algorithm = new MaTBML_I(problemSet);
            algorithm.setInputParameter("populationSize", 100);
            algorithm.setInputParameter("maxEvaluations", 1000 * 100 * taskNum);
            algorithm.setInputParameter("k1", K1);
            algorithm.setInputParameter("k2", K2);
            algorithm.setInputParameter("P_", 0.5);
            algorithm.setInputParameter("implicitTransferNum", 50);
            algorithm.setInputParameter("algoName", ALGO_NAME);

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
            parameters = new HashMap() ;
            parameters.put("comparator", new LocationComparator());
            selection = SelectionFactory.getSelectionOperator("BinaryTournament",
                    parameters);

            algorithm.addOperator("crossover", crossover);
            algorithm.addOperator("mutation", mutation);
            algorithm.addOperator("selection", selection);

            // DEBUG
            double[][] igds = new double[taskNum][times];
            for (int t = 0; t < times; t++) {
                SolutionSet population[] = algorithm.execute();
                SolutionSet resPopulation[] = new SolutionSet[taskNum];
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
//                    resPopulation[k].printObjectivesToFile("MaTBML_" + problemSet.get(k).getNumberOfObjectives() + "Obj_" +
//                            problemSet.get(k).getName() + "_" + problemSet.get(k).getNumberOfVariables() + "D_run_" + t + ".txt");

                }
                double igd;
                for (int k = 0; k < taskNum; k++){
                    QualityIndicator indicator = new QualityIndicator(problemSet.get(k), pf[k]);
                    if (population[k].size() == 0)
                        continue;

                    igd = indicator.getIGD(resPopulation[k]);
//                    ave[k] += igd;
                    // DEBUG
                    igds[k][t] = igd;
                }
//                // DEBUG
            //    LogIGD.LogIGD("MaTBML(MaOEAC-10-1)_" + problemSet.get(0).getName() + "D_run_" + t + ".txt", igds[t]);
            }
            LogIGD.LogIGD("MaTBML("+ALGO_NAME+"-"+K1+"-"+K2+")_IONLY_KL" + "_" + benchmark_name, pCase, igds);
            for(int i=0;i<taskNum;i++) {
                double[] tmp = new double[times];
                for (int t = 0; t < times; t++){
                    tmp[t] = igds[i][t];
                }
//                Arrays.sort(tmp);
                double best, worst, mean, median, std = 0;
//                best = tmp[0];
//                worst = tmp[times-1];
                mean = Arrays.stream(tmp).sum() / times;
//                median = tmp[times/2];
//                for (double e: tmp){
//                    std += Math.pow(e - mean, 2);
//                }
//                std = Math.sqrt(std / （times - 1));
//                System.out.println(best + "\t" + worst + "\t" + mean + "\t" + median + "\t" + std);
                System.out.println(form.format(mean));
            }
            System.out.println();
        }
    }
}
