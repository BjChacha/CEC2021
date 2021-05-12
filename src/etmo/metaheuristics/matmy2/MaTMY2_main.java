package etmo.metaheuristics.matmy2;

import ch.qos.logback.core.joran.spi.ElementSelector;
import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.problems.benchmarks.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.logging.LogIGD;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;

public class MaTMY2_main {
    public static void main(String[] args) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, JMException, IOException {
        ProblemSet problemSet;
        MtoAlgorithm algorithm;
        Operator crossover;
        Operator selection;

        int problemStart = 17;
        int problemEnd = 24;

        int times = 21;

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

//            // DEBUG
//            ProblemSet tmp = new ProblemSet(problemSet.size());
//            for (int i = 0; i < problemSet.size(); i ++){
//                if (i == 0)
//                    tmp.add(problemSet.get(19));
//                else if (i == 19)
//                    tmp.add(problemSet.get(0));
//                else
//                    tmp.add(problemSet.get(i));
//            }
//            problemSet = tmp;

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

//            parameters = new HashMap();
//            parameters.put("CR_LB", 0.1);
//            parameters.put("CR_UB", 0.9);
//            parameters.put("F_LB", 0.1);
//            parameters.put("F_UB", 2.0);
//            crossover = CrossoverFactory.getCrossoverOperator("RandomDECrossover",parameters);
//            algorithm.addOperator("crossover", crossover);

            parameters = null;
            selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters);

            algorithm.addOperator("crossover",crossover);
            algorithm.addOperator("selection", selection);

            double[][] igds = new double[times][taskNum];
            for (int t = 0; t < times; t++){
//                long startTime = System.currentTimeMillis();
                SolutionSet population[] = algorithm.execute();
//                long endTime = System.currentTimeMillis();
//                System.out.println("epoch: " + t + "\trunning: " + (endTime-startTime)/1000 + " s.");

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

                    // DEBUG
                    igds[t][k] = igd;
                }
                // DEBUG
                LogIGD.LogIGD("MaTMY2_" + problemSet.get(0).getName() + "D_run_" + t + ".txt", igds[t]);
//                System.out.println("Times: " + t + " finished.");
            }
            for(int i=0;i<taskNum;i++) {
                double[] tmp = new double[times];
                for (int t = 0; t < times; t++){
                    tmp[t] = igds[t][i];
                }
                Arrays.sort(tmp);
                double best, worst, mean, median, std = 0;
                best = tmp[0];
                worst = tmp[times-1];
                mean = Arrays.stream(tmp).sum() / times;
                median = tmp[times/2];
                for (double e: tmp){
                    std += Math.pow(e - mean, 2);
                }
                std = Math.sqrt(std / (times-1));
                System.out.println(best + "\t" + worst + "\t" + mean + "\t" + median + "\t" + std);
//                System.out.println(form.format(ave[i] / times));
            }
            System.out.println();
        }
    }
}
