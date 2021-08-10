package etmo.metaheuristics.matde;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;

import etmo.core.MtoAlgorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.operators.crossover.CrossoverFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.problems.CEC2017.*;
import etmo.util.logging.LogIGD;

public class MaTDE_main {
    public static void main(String[] args) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, JMException, IOException {
        ProblemSet problemSet;
        MtoAlgorithm algorithm;
        Operator crossover1;
        Operator crossover2;

        HashMap parameters;

        int problemStart = 1;
        int problemEnd = 9;

        int times = 20;

        String benchmark_name;
        DecimalFormat form = new DecimalFormat("#.####E0");

        System.out.println("Algo: MaTDE.");

        for (int pCase = problemStart; pCase <= problemEnd; pCase++){
//           // CEC 2021
//           benchmark_name = "CEC2021";
//           problemSet = (ProblemSet) Class
//                   .forName("etmo.problems.benchmarks_CEC2021.ETMOF" + pCase)
//                   .getMethod("getProblem")
//                   .invoke(null, null);

//             // WCCI 2020
//             benchmark_name = "WCCI2020";
//             problemSet = (ProblemSet) Class
//                     .forName("etmo.problems.benchmarks_WCCI2020.MATP" + pCase)
//                     .getMethod("getProblem")
//                     .invoke(null, null);

             // CEC2017
             benchmark_name = "CEC2017";
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
            for (int k = 0; k < pf.length; k++){
                pf[k] = "resources/PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
            }

            String pSName = problemSet.get(0).getName();
            pSName = pSName.substring(0, pSName.length()-2);
            System.out.println(pSName + "\ttaskNum = "+taskNum+"\tfor "+times+" times.");

            algorithm = new MaTDE(problemSet);

            algorithm.setInputParameter("populationSize", 100);
            algorithm.setInputParameter("archiveSize", 300);
            algorithm.setInputParameter("maxEvaluations", 1000 * taskNum * 100);

            // 迁移交叉概率
            algorithm.setInputParameter("alpha", 0.1);
            // 衰减系数
            algorithm.setInputParameter("ro", 0.8);
            // 衰减速率
            algorithm.setInputParameter("shrinkRate", 0.8);
            // Archive更新概率
            algorithm.setInputParameter("replaceRate", 0.2);

            parameters = new HashMap();
            // 原论文：CR = (0.1, 0.9)
            parameters.put("CR_LB", 0.1);
            parameters.put("CR_UB", 0.9);
            // 原论文：F = (0.1, 2)
            parameters.put("F_LB", 0.1);
            parameters.put("F_UB", 2.0);
            crossover1 = CrossoverFactory.getCrossoverOperator("RandomDECrossover",parameters);
            algorithm.addOperator("crossover1", crossover1);

            parameters = new HashMap();
            parameters.put("CR_LB", 0.1);
            parameters.put("CR_UB", 0.9);
            crossover2 = CrossoverFactory.getCrossoverOperator("RandomUniformCrossover", parameters);
            algorithm.addOperator("crossover2", crossover2);

            double[][] IGDs = new double[taskNum][times];
            for (int t = 0; t < times; t++){
                long startTime = System.currentTimeMillis();
                SolutionSet[] population = algorithm.execute();
                long endTime = System.currentTimeMillis();
                System.out.println("epoch: " + t + "\trunning: " + (endTime-startTime)/1000 + " s.");

                // 各个任务的目标数不一，所以评价时需要根据任务来重新设置种群规模。
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
//                    resPopulation[k].printObjectivesToFile("MaTDE_"+problemSet.get(k).getNumberOfObjectives()+"Obj_"+
//                            problemSet.get(k).getName()+ "_" + problemSet.get(k).getNumberOfVariables() + "D_run_"+t+".txt");
                }
                double igd;
                for (int k = 0; k < taskNum; k++){
                    QualityIndicator indicator = new QualityIndicator(problemSet.get(k), pf[k]);
                    if (population[k].size() == 0)
                        continue;

                    igd = indicator.getIGD(resPopulation[k]);
                    IGDs[k][t] = igd;
                    ave[k] += igd;
                }
            }
            LogIGD.LogIGD("MaTDE_" + "x" + times + "_" + benchmark_name, pCase, IGDs);
            for(int i=0;i<taskNum;i++)
                System.out.println(form.format(ave[i] / times));
            System.out.println();
        }
    }
}
