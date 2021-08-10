package etmo.metaheuristics.mteaad;

import etmo.core.*;
import etmo.operators.crossover.Crossover;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class MTEA_AD_main {
    public static void main(String[] args) throws IOException, JMException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
        ProblemSet problemSet; // The problem to solve
        Algorithm algorithm; // The algorithm to use
        Crossover crossover; // Crossover operator
        Operator mutation; // Mutation operator
        Operator selection;

        HashMap parameters; // Operator parameters

        for (int pCase = 1; pCase <= 10; pCase++ ){
            // switch (pCase){
            //     case 1:
            //         problemSet = ETMOF1.getProblem();
            //         break;
            //     case 2:
            //         problemSet = ETMOF2.getProblem();
            //         break;
            //     case 3:
            //         problemSet = ETMOF3.getProblem();
            //         break;
            //     case 4:
            //         problemSet = ETMOF4.getProblem();
            //         break;
            //     case 5:
            //         problemSet = ETMOF5.getProblem();
            //         break;
            //     case 6:
            //         problemSet = ETMOF6.getProblem();
            //         break;
            //     case 7:
            //         problemSet = ETMOF7.getProblem();
            //         break;
            //     case 8:
            //         problemSet = ETMOF8.getProblem();
            //         break;
            //     case 9:
            //         problemSet = ETMOF9.getProblem();
            //         break;
            //     case 10:
            //         problemSet = ETMOF10.getProblem();
            //         break;
            //     case 11:
            //         problemSet = ETMOF11.getProblem();
            //         break;
            //     case 12:
            //         problemSet = ETMOF12.getProblem();
            //         break;
            //     case 13:
            //         problemSet = ETMOF13.getProblem();
            //         break;
            //     case 14:
            //         problemSet = ETMOF14.getProblem();
            //         break;
            //     case 15:
            //         problemSet = ETMOF15.getProblem();
            //         break;
            //     case 16:
            //         problemSet = ETMOF16.getProblem();
            //         break;
            //     default:
            //         problemSet = ETMOF1.getProblem();
            // }

//            // CEC2021
//            problemSet = (ProblemSet) Class
//                    .forName("etmo.problems.benchmarks_CEC2021.ETMOF" + pCase)
//                    .getMethod("getProblem")
//                    .invoke(null, null);

            // WCCI 2020
            problemSet = (ProblemSet) Class
                    .forName("etmo.problems.benchmarks_WCCI2020.MATP" + pCase)
                    .getMethod("getProblem")
                    .invoke(null, null);


            int taskNumber = problemSet.size();
            System.out.println("taskNumber = "+taskNumber);
            System.out.println("TaskID\t" + "IGD for " + problemSet.get(0).getName()+" to " +problemSet.get(taskNumber-1).getName());
            algorithm = new MTEA_AD(problemSet);

            String[] pf = new String[taskNumber];
            for (int i = 0; i < pf.length; i++){
                pf[i] = "resources/PF/StaticPF/" + problemSet.get(i).getHType() + "_" + problemSet.get(i).getNumberOfObjectives() + "D.pf";
            }//String pf = "PF/StaticPF/" + "convex.pf";
            //System.out.println(pf);
            algorithm.setInputParameter("populationSize", 100);
            algorithm.setInputParameter("maxEvaluations",1000 * 100 * taskNumber);

            parameters = new HashMap();
            parameters.put("probability", 0.9);
            parameters.put("distributionIndex", 20.0);
//                crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);
            crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

            //Crossover Operators
//                parameters = new HashMap();
//                parameters.put("distributionIndex", 30.0);
//                parameters.put("Lr", 0.2);
//                parameters.put("Dr", 0.2);
//                parameters.put("Er", 0.7);
//                crossover = CrossoverFactory.getCrossoverOperator("EGG",parameters);

            // Mutation operator
            parameters = new HashMap();
            parameters.put("probability", 1.0 / problemSet.getMaxDimension());
            parameters.put("distributionIndex", 20.0);
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

            // Selection Operator
            parameters = null ;
            selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters) ;

            // Add the operators to the algorithm
            algorithm.addOperator("crossover", crossover);
            algorithm.addOperator("mutation", mutation);
            algorithm.addOperator("selection", selection);

//                System.out.println("RunID\t" + "IGD for " + problemSet2.get(0).getName());
            DecimalFormat form = new DecimalFormat("#.####E0");


            int times = 1;
            double ave[] = new double[taskNumber];
            double cpIGD[][] = new double[taskNumber][times];

            for (int t = 1; t <= times; t++) {
                SolutionSet[] resPopulation = ((MTEA_AD) algorithm).executeMultiTask();


                SolutionSet[] res = new SolutionSet[taskNumber];
                for (int i = 0; i < taskNumber; i++){
                    res[i] = new SolutionSet();
                    for (int k = 0; k < resPopulation[i].size(); k++){
                        Solution sol = resPopulation[i].get(k);
                        int start = problemSet.get(i).getStartObjPos();
                        int end = problemSet.get(i).getEndObjPos();
                        Solution newSolution = new Solution(end - start + 1);

                        for (int l = start; l <= end; l++)
                            newSolution.setObjective(l - start, sol.getObjective(l));
                        res[i].add(newSolution);
                    }
                }

                double igd;
//				System.out.print(t + "\t");
                for(int i = 0; i < taskNumber; i++){
                    QualityIndicator indicator = new QualityIndicator(problemSet.get(i), pf[i]);
                    if(res[i].size()==0)
                        continue;
                    igd =  indicator.getIGD(res[i]);
//					System.out.print(form.format(igd) + "\t" );
                    ave[i] += igd;
                    cpIGD[i][t-1] = igd;
                }
//                    System.out.println(i + "\t" + form.format(igd));
            }

            for(int i=0;i<taskNumber;i++)
                System.out.println(form.format(ave[i] / times));

//            String path = "MaOEA_ACT_2021F9-16.txt";
//            printIGD.printIGDtoText(path, cpIGD, taskNumber, times);


        }



    }


}
