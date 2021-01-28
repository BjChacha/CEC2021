package etmo.metaheuristics.matde;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;

import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class MaTDE_main {
    public static void main(String[] args) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, JMException {
        ProblemSet problemSet;
        MtoAlgorithm algorithm;
        Operator crossover;

        HashMap parameters;

        int problemStart = 25;
        int problemEnd = 32;

        int times = 21;

        DecimalFormat form = new DecimalFormat("#.####E0");

        System.out.println("Algo: MaTDE.");

        for (int pCase = problemStart; pCase <= problemEnd; pCase++){
            problemSet = (ProblemSet) Class
                    .forName("etmo.problems.benchmarks_ETMO.ETMOF" + pCase)
                    .getMethod("getProblem")
                    .invoke(null, null);

            int taskNum = problemSet.size();
            double[] ave = new double[taskNum];

            String[] pf = new String[taskNum];
            for (int k = 0; k < pf.length; k++){
                pf[k] = "PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
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
            // 子代个体加入到Archive中的概率
            algorithm.setInputParameter("replaceRate", 0.2);

            parameters = new HashMap();
            // 原论文：CR = (0.1, 0.9)
            parameters.put("CR", 0.55);
            // 原论文：F = (0.1, 2)
            parameters.put("F", 1.0);
            crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover",parameters);

            algorithm.addOperator("crossover", crossover);

            for (int t = 0; t < times; t++){
                long startTime = System.currentTimeMillis();

                SolutionSet population[] = algorithm.execute();

                long endTime = System.currentTimeMillis();
                System.out.println("epoch: " + t + "\trunning: " + (endTime-startTime)/1000 + " s.");

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
                System.out.println("T" + (i+1) + "\t" + form.format(ave[i] / times));
            System.out.println();
        }
    }
}