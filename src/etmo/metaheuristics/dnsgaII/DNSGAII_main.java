package etmo.metaheuristics.dnsgaII;

import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.problems.benchmarks_ETMO.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.Ranking;

public class DNSGAII_main {
	public static void main(String args[]) throws IOException, JMException, ClassNotFoundException {
		ProblemSet problemSet1; // The problem to solve
		ProblemSet problemSet2;
		DynamicAlgorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection;
		HashMap parameters; // Operator parameters

		/*
		 * fc indicates the frequency of change
		 * sc indicates the severity of change
		 * nc number of changes in each run
		*/
		int fc = 20;
		int sc = 10;
		int nc = 30;

		problemSet1 = ETMOF34.getProblem(fc,sc);
		int taskNumber = problemSet1.size();
		System.out.println("taskNumber = "+taskNumber);
		for (int tsk=0;tsk<taskNumber;tsk++) {
			 
			problemSet2 = problemSet1.getTask(tsk);
			algorithm = new DNSGAII(problemSet2);
			algorithm.setInputParameter("populationSize", 100);
		
			algorithm.setInputParameter("maxEvaluations", 100*fc*(nc+1));
	
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
	
			System.out.println("RunID\t" + "MIGD for " + problemSet2.get(0).getName());
			DecimalFormat form = new DecimalFormat("#.####E0");
			
			String[] pf = new String[nc+1];
			for(int i=0;i<nc+1;i++){
				pf[i] = "PF/DynamicPF/" + problemSet2.get(0).getName() + "/POF_Tt=" + i + ".txt";
			}
			
			QualityIndicator indicator;
			int times = 1;
			double aveIGD = 0;
				for (int i = 1; i <= times; i++) {
					double MIGD = 0;
					List<SolutionSet> dynamicPopulations = algorithm.execute();
					for(int p=0;p<dynamicPopulations.size();p++){
						indicator = new QualityIndicator(problemSet2.get(0), pf[p]);
						dynamicPopulations.get(p).printObjectivesToFile("DNSGAII_"+problemSet2.get(0).getNumberOfObjectives()+"Obj_"+
                                problemSet2.get(0).getName()+ "_" + problemSet2.get(0).getNumberOfVariables() + "D_run"+i+"_"+p+".txt");
						MIGD += indicator.getIGD(dynamicPopulations.get(p));
					}
					MIGD = MIGD/dynamicPopulations.size();
					aveIGD += MIGD;
					System.out.println(i + "\t" + form.format(MIGD));
				} 					
				System.out.println("Average MIGD for " + problemSet2.get(0).getName() + ": " + form.format(aveIGD / times));
				System.out.println();
			}
	}	

}
