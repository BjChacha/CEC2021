package etmo.metaheuristics.maoeac;

import etmo.core.Algorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.SolutionSet;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.logging.LogIGD;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class MaOEAC_main {
    public static void main(String args[]) throws JMException, ClassNotFoundException, SecurityException, IOException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
        ProblemSet problemSet; // The problem to solve
        Algorithm algorithm; // The algorithm to use
        Operator crossover; // Crossover operator
        Operator mutation; // Mutation operator
        Operator selection; //Selection operator

        HashMap parameters;

        int taskStart = 1;
        int taskEnd = 10;

        int times = 10;

        DecimalFormat form = new DecimalFormat("#.####E0");
        String benchmark_name;
        for (int pCase = taskStart; pCase <= taskEnd; pCase++){
//            benchmark_name = "CEC2021";
//            problemSet = (ProblemSet) Class
//                    .forName("etmo.problems.benchmarks_CEC2021.ETMOF" + pCase)
//                    .getMethod("getProblem")
//                    .invoke(null, null);

             // WCCI 2020
             benchmark_name = "WCCI2020";
             problemSet = (ProblemSet) Class
                     .forName("etmo.problems.benchmarks_WCCI2020.MATP" + pCase)
                     .getMethod("getProblem")
                     .invoke(null, null);

            int taskNum = problemSet.size();

            String[] pf = new String[taskNum];
            double[] ave = new double[taskNum];

            String pSName = problemSet.get(0).getName();
            pSName = pSName.substring(0, pSName.length()-2);
            System.out.println(pSName + "\ttaskNum = "+taskNum+"\tfor "+times+" times.");

            double[][] igds = new double[taskNum][times];
            for (int tsk = 0; tsk < taskNum; tsk++){
                ProblemSet pS = problemSet.getTask(tsk);
//                System.out.println("RunID\t" + "IGD for " + problemSet.get(tsk).getName() + " for " + times + " times.");
                pf[tsk] = "resources/PF/StaticPF/" + problemSet.get(tsk).getHType() + "_" + problemSet.get(tsk).getNumberOfObjectives() + "D.pf";

                for (int t = 1; t <= times; t++){
                    algorithm = new MaOEAC(pS);

                    algorithm.setInputParameter("populationSize",100);
                    algorithm.setInputParameter("maxGenerations",1000);

                    parameters = new HashMap();
                    parameters.put("probability", 1.0);
                    parameters.put("distributionIndex", 30.0);
                    crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
                            parameters);

                    parameters = new HashMap();
                    parameters.put("probability", 1.0 / pS.get(0).getNumberOfVariables());
                    parameters.put("distributionIndex", 20.0);
                    mutation = MutationFactory.getMutationOperator("PolynomialMutation",
                            parameters);

                    parameters = null;
                    selection = SelectionFactory.getSelectionOperator("RandomSelection",
                            parameters);

                    algorithm.addOperator("crossover", crossover);
                    algorithm.addOperator("mutation", mutation);
                    algorithm.addOperator("selection", selection);

                    QualityIndicator indicator = new QualityIndicator(problemSet.get(tsk), pf[tsk]);

                    SolutionSet population = algorithm.execute();
                    population.printObjectivesToFile("MaOEAC_"+problemSet.get(tsk).getNumberOfObjectives()+"Obj_"+
                            problemSet.get(tsk).getName()+ "_" + problemSet.get(tsk).getNumberOfVariables() + "D_run"+t+".txt");
                    double igd =  indicator.getIGD(population);
                    igds[tsk][t-1] = igd;
                    ave[tsk] += igd;
                }
                // System.out.println("T" + (tsk+1) + "\t" + form.format(ave[tsk] / times));
                System.out.println(form.format(ave[tsk] / times));
            }
            LogIGD.LogIGD("MaOEAC" + "_" + benchmark_name, pCase, igds);
            System.out.println();
            // for briefly summarization
//			for(int i=0;i<taskNum;i++) {
//				System.out.println("Average IGD for " + problemSet.get(i).getName() + ": " + form.format(ave[i] / times));
//			}
        }//for-fun
    }//main
}
