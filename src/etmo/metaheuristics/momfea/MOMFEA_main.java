package etmo.metaheuristics.momfea;

import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;

import etmo.util.comparators.LocationComparator;
import etmo.core.Algorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.problems.benchmarks_ETMO.*;

import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;

import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;

public class MOMFEA_main {
	public static void main(String args[]) throws IOException, JMException, ClassNotFoundException {
		ProblemSet problemSet; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection;

		HashMap parameters; // Operator parameters


		for (int pCase = 9; pCase <= 16; pCase++ ){
			switch (pCase){
				case 1:
					problemSet = ETMOF1.getProblem();
					break;
				case 2:
					problemSet = ETMOF2.getProblem();
					break;
				case 3:
					problemSet = ETMOF3.getProblem();
					break;
				case 4:
					problemSet = ETMOF4.getProblem();
					break;
				case 5:
					problemSet = ETMOF5.getProblem();
					break;
				case 6:
					problemSet = ETMOF6.getProblem();
					break;
				case 7:
					problemSet = ETMOF7.getProblem();
					break;
				case 8:
					problemSet = ETMOF8.getProblem();
					break;
				case 9:
					problemSet = ETMOF9.getProblem();
					break;
				case 10:
					problemSet = ETMOF10.getProblem();
					break;
				case 11:
					problemSet = ETMOF11.getProblem();
					break;
				case 12:
					problemSet = ETMOF12.getProblem();
					break;
				case 13:
					problemSet = ETMOF13.getProblem();
					break;
				case 14:
					problemSet = ETMOF14.getProblem();
					break;
				case 15:
					problemSet = ETMOF15.getProblem();
					break;
				case 16:
					problemSet = ETMOF16.getProblem();
					break;
				default:
					problemSet = ETMOF1.getProblem();
			}

			int taskNumber = problemSet.size();
			System.out.println("taskNumber = "+taskNumber);


			String[] pf = new String[problemSet.size()];
			for (int i = 0; i < pf.length; i++){
				pf[i] = "PF/StaticPF/" + problemSet.get(i).getHType() + "_" + problemSet.get(i).getNumberOfObjectives() + "D.pf";
			}

			algorithm = new MOMFEA(problemSet);

			algorithm.setInputParameter("populationSize",100*taskNumber);
			algorithm.setInputParameter("maxEvaluations",100*taskNumber * 1000);
			algorithm.setInputParameter("rmp", 0.9);

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

			// Add the operators to the algorithm
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);
			algorithm.addOperator("selection", selection);

			DecimalFormat form = new DecimalFormat("#.####E0");

			System.out.println("RunID\t" + "IGD for "+problemSet.get(0).getName()+" to "+problemSet.get(taskNumber-1).getName());

			int times = 21;

			double ave[] = new double[taskNumber];
			for (int t = 1; t <= times; t++) {
				SolutionSet population = algorithm.execute();

				SolutionSet[] resPopulation = new SolutionSet[problemSet.size()];
				for (int i = 0; i < problemSet.size(); i++)
					resPopulation[i] = new SolutionSet();

				for (int i = 0; i < population.size(); i++) {
					Solution sol = population.get(i);

					int pid = sol.getSkillFactor();

					int start = problemSet.get(pid).getStartObjPos();
					int end = problemSet.get(pid).getEndObjPos();

					Solution newSolution = new Solution(end - start + 1);

					for (int k = start; k <= end; k++)
						newSolution.setObjective(k - start, sol.getObjective(k));

					resPopulation[pid].add(newSolution);
				}

				//测试种群中不同任务的个体数量
//			for (int i = 0; i < problemSet.size(); i++){
//				System.out.print(i+":"+resPopulation[i].size()+"\t");
//			}
//			System.out.println("");

				double igd;
				System.out.print(t + "\t");
				for(int i = 0; i < taskNumber; i++){
					QualityIndicator indicator = new QualityIndicator(problemSet.get(i), pf[i]);
					if(resPopulation[i].size()==0)
						continue;
//				getTask中用到add影响problem起始和结束值
//				resPopulation[i].printObjectivesToFile("MOMFEA_"+problemSet.getTask(i).get(0).getNumberOfObjectives()+"Obj_"+
//						problemSet.getTask(i).get(0).getName()+ "_" + problemSet.getTask(i).get(0).getNumberOfVariables() + "D_run"+t+".txt");
//					resPopulation[i].printObjectivesToFile("MOMFEA_"+problemSet.get(i).getNumberOfObjectives()+"Obj_"+
//							problemSet.get(i).getName()+ "_" + problemSet.get(i).getNumberOfVariables() + "D_run"+t+".txt");

					igd =  indicator.getIGD(resPopulation[i]);
					System.out.print(form.format(igd) + "\t" );
					ave[i] += igd;
				}
				System.out.println("");
			}

			System.out.println();
			for(int i=0;i<taskNumber;i++)
				System.out.println("Average IGD for " + problemSet.get(i).getName()+ ": " + form.format(ave[i] / times));

		}



	}
}
