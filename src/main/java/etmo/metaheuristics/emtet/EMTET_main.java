package etmo.metaheuristics.emtet;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;

import etmo.core.*;
import etmo.problems.CEC2017.*;
import etmo.util.comparators.LocationComparator;

import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;

import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.logging.LogIGD;

public class EMTET_main {
    public static void main(String[] args) throws IOException, JMException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
        ProblemSet problemSet; // The problem to solve
        MtoAlgorithm algorithm; // The algorithm to use
        Operator crossover; // Crossover operator
        Operator mutation; // Mutation operator
        Operator selection;

        HashMap parameters; // Operator parameters

        int problemStart = 1;
        int problemEnd = 9;

        int times = 20;

        DecimalFormat form = new DecimalFormat("#.####E0");

        System.out.println("Algo: EMTET.");

        for (int pCase = problemStart; pCase <= problemEnd; pCase++) {
//            problemSet = (ProblemSet) Class
//                    .forName("etmo.problems.benchmarks_CEC2021.ETMOF" + pCase)
//                    .getMethod("getProblem")
//                    .invoke(null, null);
            // CEC2017
            ProblemSet[] cec2017 = {
                    CIHS.getProblem(),
                    CIMS.getProblem(),
                    CILS.getProblem(),
                    PIHS.getProblem(),
                    PIMS.getProblem(),
                    PILS.getProblem(),
                    NIHS.getProblem(),
                    NIMS.getProblem(),
                    NILS.getProblem()
            };
            problemSet = cec2017[pCase-1];


            int taskNum = problemSet.size();
            double[] ave = new double[taskNum];

            String[] pf = new String[taskNum];
            for (int i = 0; i < pf.length; i++) {
                pf[i] = "resources/PF/StaticPF/" + problemSet.get(i).getHType() + "_" + problemSet.get(i).getNumberOfObjectives() + "D.pf";
            }

            String pSName = problemSet.get(0).getName();
            pSName = pSName.substring(0, pSName.length() - 2);
            System.out.println(pSName + "\ttaskNum = "+taskNum+"\tfor "+times+" times.");

            algorithm = new EMTET(problemSet);

            algorithm.setInputParameter("populationSize", 100);
            algorithm.setInputParameter("maxEvaluations", 100 * taskNum * 1000);
            algorithm.setInputParameter("transferNum", 8);

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

            // Add the operators to the algorithm
            algorithm.addOperator("crossover", crossover);
            algorithm.addOperator("mutation", mutation);
            algorithm.addOperator("selection", selection);

            double[][] IGDs = new double[taskNum][times];
            for (int t = 1; t <= times; t++) {
                long startTime = System.currentTimeMillis();

                SolutionSet[] res = algorithm.execute();

                long endTime = System.currentTimeMillis();
                System.out.println("epoch: " + t + "\trunning: " + (endTime-startTime)/1000 + " s.");

                SolutionSet[] resPopulation = new SolutionSet[taskNum];
                for (int k = 0; k < taskNum; k++){
                    resPopulation[k] = new SolutionSet();
                    for (int i = 0; i < res[k].size(); i++) {
                        Solution sol = res[k].get(i);

                        int start = problemSet.get(k).getStartObjPos();
                        int end = problemSet.get(k).getEndObjPos();

                        Solution newSolution = new Solution(end - start + 1);

                        for (int kk = start; kk <= end; kk++)
                            newSolution.setObjective(kk - start, sol.getObjective(kk));

                        resPopulation[k].add(newSolution);
                    }
                    resPopulation[k].printObjectivesToFile("EMTET_"+problemSet.get(k).getNumberOfObjectives()+"Obj_"+
                            problemSet.get(k).getName()+ "_" + problemSet.get(k).getNumberOfVariables() + "D_run_"+t+".txt");
                }
                double igd;
                for (int k = 0; k < taskNum; k++) {
                    QualityIndicator indicator = new QualityIndicator(problemSet.get(k), pf[k]);
                    if (resPopulation[k].size() == 0)
                        continue;

                    igd = indicator.getIGD(resPopulation[k]);
                    IGDs[k][t-1] = igd;
                    ave[k] += igd;
                }
            }
            LogIGD.LogIGD("EMTET_CEC2017_x" + times, pCase, IGDs);
            for (int i = 0; i < taskNum; i++)
                System.out.println("Average IGD for " + problemSet.get(i).getName() + ": " + form.format(ave[i] / times));
            System.out.println();
        }
    }
}
