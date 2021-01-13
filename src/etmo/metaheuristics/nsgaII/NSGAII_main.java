package etmo.metaheuristics.nsgaII;

import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.problems.benchmarks_ETMO.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.Ranking;

public class NSGAII_main {
	public static void main(String args[]) throws IOException, JMException, ClassNotFoundException {
		ProblemSet problemSet1; // The problem to solve
		ProblemSet problemSet2;
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection;

		HashMap parameters; // Operator parameters

		problemSet1 = ETMOF3.getProblem();
		int taskNumber = problemSet1.size();
		System.out.println("taskNumber = "+taskNumber);
		for (int tsk=0;tsk<taskNumber;tsk++) {
			 
			problemSet2 = problemSet1.getTask(tsk);
			algorithm = new NSGAII(problemSet2);
	
			String pf = "PF/StaticPF/" + problemSet2.get(0).getHType() + "_" + problemSet2.get(0).getNumberOfObjectives() + "D.pf";
			//String pf = "PF/StaticPF/" + "convex.pf";
		    //System.out.println(pf);	
			algorithm.setInputParameter("populationSize", 100);
			algorithm.setInputParameter("maxEvaluations", 100 * 1000);
	
			parameters = new HashMap();
			parameters.put("probability", 0.9);
			parameters.put("distributionIndex", 20.0);
			crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
	
			// Mutation operator
			parameters = new HashMap();
			parameters.put("probability", 1.0 / problemSet2.getMaxDimension());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
	
			// Selection Operator
		    parameters = null ;
		    selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters) ;  
		    
			// Add the operators to the algorithm
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);
			algorithm.addOperator("selection", selection);
	
			System.out.println("RunID\t" + "IGD for " + problemSet2.get(0).getName());
			DecimalFormat form = new DecimalFormat("#.####E0");
			QualityIndicator indicator = new QualityIndicator(problemSet2.get(0), pf);
			int times = 1;
			double aveIGD = 0;
				for (int i = 1; i <= times; i++) {
					SolutionSet population = algorithm.execute();
					Ranking ranking = new Ranking(population);
					population = ranking.getSubfront(0);
					population.printObjectivesToFile("NSGAII_"+problemSet2.get(0).getNumberOfObjectives()+"Obj_"+
					                                 problemSet2.get(0).getName()+ "_" + problemSet2.get(0).getNumberOfVariables() + "D_run"+i+".txt");
					double igd = indicator.getIGD(population);
					aveIGD += igd;
					System.out.println(i + "\t" + form.format(igd));
				} 					
				System.out.println("Average IGD for " + problemSet2.get(0).getName() + ": " + form.format(aveIGD / times));
				System.out.println();
			}
	}	

}
