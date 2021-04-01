package etmo.metaheuristics.matmy2;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.problems.benchmarks.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class MaTMY2_main {
    public static void main(String[] args) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, JMException, IOException {
        ProblemSet problemSet;
        MtoAlgorithm algorithm;
        Operator crossover;
        Operator selection;

        int problemStart = 25;
        int problemEnd = 25;

        int times = 1;

        DecimalFormat form = new DecimalFormat("#.####E0");

        System.out.println("Algo: MaTMY2.");

        for (int pCase = problemStart; pCase <= problemEnd; pCase++){
            // CEC2021
            problemSet = (ProblemSet) Class
                    .forName("etmo.problems.benchmarks_ETMO.ETMOF" + pCase)
                    .getMethod("getProblem")
                    .invoke(null, null);

            // CEC 2017
//            switch (pCase){
//                case 2:
//                    problemSet = CIMS.getProblem();
//                    break;
//                case 3:
//                    problemSet = CILS.getProblem();
//                    break;
//                case 4:
//                    problemSet = PIHS.getProblem();
//                    break;
//                case 5:
//                    problemSet = PIMS.getProblem();
//                    break;
//                case 6:
//                    problemSet = PILS.getProblem();
//                    break;
//                case 7:
//                    problemSet = NIHS.getProblem();
//                    break;
//                case 8:
//                    problemSet = NIMS.getProblem();
//                    break;
//                case 9:
//                    problemSet = NILS.getProblem();
//                    break;
//                default:
//                    problemSet = CIHS.getProblem();
//            }


            int taskNum = problemSet.size();
            double[] ave = new double[taskNum];

            String[] pf = new String[taskNum];
            for (int k = 0; k < pf.length; k++){
                pf[k] = "PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
            }

            String pSName = problemSet.get(0).getName();
//            pSName = pSName.substring(0, pSName.length()-2);
            System.out.println(pSName + "\ttaskNum = "+taskNum+"\tfor "+times+" times.");

            algorithm = new MaTMY2(problemSet);

            algorithm.setInputParameter("populationSize", 100);
            algorithm.setInputParameter("maxEvaluations", 1000 * 100 * taskNum);
            algorithm.setInputParameter("transferVolume", 10);
            algorithm.setInputParameter("baseRunTime", 3);

            algorithm.setInputParameter("forceTransferRate", 0.2);

            algorithm.setInputParameter("scoreIncrement", 0.5);
            algorithm.setInputParameter("scoreDecreaseRate", 0.2);
            algorithm.setInputParameter("isDRA", false);
            algorithm.setInputParameter("algoName", "MOEAD");

            HashMap parameters;

            parameters = new HashMap();
            parameters.put("probability", 0.9);
            parameters.put("distributionIndex", 20.0);
            crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

            parameters = null;
            selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters);

            algorithm.addOperator("crossover",crossover);
            algorithm.addOperator("selection", selection);

            for (int t = 0; t < times; t++){
//                long startTime = System.currentTimeMillis();
                SolutionSet population[] = algorithm.execute();
//                long endTime = System.currentTimeMillis();
//                System.out.println("epoch: " + t + "\trunning: " + (endTime-startTime)/1000 + " s.");

                // 各个任务的目标数不一，所以评价时需要根据任务来重新设置种群规模。
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
                    resPopulation[k].printObjectivesToFile("MaTMY2_"+problemSet.get(k).getNumberOfObjectives()+"Obj_"+
                            problemSet.get(k).getName()+ "_" + problemSet.get(k).getNumberOfVariables() + "D_run_"+t+".txt");
                }
                double igd;
                for (int k = 0; k < taskNum; k++){
                    QualityIndicator indicator = new QualityIndicator(problemSet.get(k), pf[k]);
                    if (population[k].size() == 0)
                        continue;

                    igd = indicator.getIGD(resPopulation[k]);
                    ave[k] += igd;
                }
            }
            for(int i=0;i<taskNum;i++)
                System.out.println(form.format(ave[i] / times));
            System.out.println();
        }
    }
}
