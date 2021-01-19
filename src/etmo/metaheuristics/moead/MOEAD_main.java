package etmo.metaheuristics.moead;

import etmo.core.Algorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.SolutionSet;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.problems.benchmarks_ETMO.ETMOF1;
import etmo.problems.benchmarks_ETMO.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.Ranking;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class MOEAD_main {
    public static void main(String args[]) throws IOException, JMException, ClassNotFoundException {
        ProblemSet problemSet1; // The problem to solve
        ProblemSet problemSet2;

        Algorithm algorithm; // The algorithm to use
        Operator crossover; // Crossover operator
        Operator mutation; // Mutation operator
        Operator selection;

        HashMap parameters; // Operator parameters

        for (int pCase = 9; pCase <= 9; pCase++ ){
            switch (pCase){
                case 1:
                    problemSet1 = ETMOF1.getProblem();
                    break;
                case 2:
                    problemSet1 = ETMOF2.getProblem();
                    break;
                case 3:
                    problemSet1 = ETMOF3.getProblem();
                    break;
                case 4:
                    problemSet1 = ETMOF4.getProblem();
                    break;
                case 5:
                    problemSet1 = ETMOF5.getProblem();
                    break;
                case 6:
                    problemSet1 = ETMOF6.getProblem();
                    break;
                case 7:
                    problemSet1 = ETMOF7.getProblem();
                    break;
                case 8:
                    problemSet1 = ETMOF8.getProblem();
                    break;
                case 9:
                    problemSet1 = ETMOF9.getProblem();
                    break;
                default:
                    problemSet1 = ETMOF1.getProblem();
            }

            int taskNumber = problemSet1.size();
            System.out.println("taskNumber = "+taskNumber);
            for (int tsk=0; tsk < taskNumber; tsk++) {

                problemSet2 = problemSet1.getTask(tsk);
                algorithm = new MOEAD(problemSet2);

                String pf = "PF/StaticPF/" + problemSet2.get(0).getHType() + "_" + problemSet2.get(0).getNumberOfObjectives() + "D.pf";

                algorithm.setInputParameter("populationSize", 495);
                algorithm.setInputParameter("maxEvaluations", 495 * 1000);

                algorithm.setInputParameter("dataDirectory", "D:\\Workspace\\EMTO2021\\myRes\\MTO-cec2021-\\resources\\weightVectorFiles\\moead");


                algorithm.setInputParameter("T", 20);
                algorithm.setInputParameter("delta", 0.9);
                algorithm.setInputParameter("nr", 2);

                parameters = new HashMap();
                parameters.put("CR", 1.0);
                parameters.put("F", 0.5);
                crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover",parameters);

                // Mutation operator
                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet2.get(0).getNumberOfVariables());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);


                algorithm.addOperator("crossover", crossover);
                algorithm.addOperator("mutation", mutation);

                System.out.println("RunID\t" + "IGD for " + problemSet2.get(0).getName());
                DecimalFormat form = new DecimalFormat("#.####E0");
                QualityIndicator indicator = new QualityIndicator(problemSet2.get(0), pf);
                int times = 21;
                double aveIGD = 0;
                for (int i = 1; i <= times; i++) {
                    SolutionSet population = algorithm.execute();
//                population.printObjectivesToFile("MOEAD_"+problemSet.get(0).getNumberOfObjectives()+"Obj_"+
//                        problemSet.get(0).getName()+ "_" + problemSet.get(0).getNumberOfVariables() + "D_run"+i+".txt");
                    double igd = indicator.getIGD(population);
                    aveIGD += igd;
//                    System.out.println(i + "\t" + form.format(igd));
                }
                System.out.println("Average IGD for " + problemSet2.get(0).getName() + ": " + form.format(aveIGD / times));
                System.out.println();




            }

        }





        }
}
