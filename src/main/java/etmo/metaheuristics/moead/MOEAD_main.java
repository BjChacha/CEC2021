package etmo.metaheuristics.moead;

import etmo.core.Algorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.SolutionSet;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.logging.LogIGD;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;

public class MOEAD_main {
	public static void main(String[] args) throws JMException, SecurityException, IOException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
		ProblemSet problemSet; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator

		HashMap parameters; // Operator parameters

		int taskStart = 1;
		int taskEnd = 10;

		int times = 10;

		DecimalFormat form = new DecimalFormat("#.####E0");
		String benchmark_name;
		for (int pCase = taskStart; pCase <= taskEnd; pCase++){
			// benchmark_name = "CEC2021";
			// problemSet = (ProblemSet) Class
			// 		.forName("etmo.problems.CEC2021.ETMOF" + pCase)
			// 		.getMethod("getProblem")
			// 		.invoke(null, null);

			// WCCI 2020
			benchmark_name = "WCCI2020";
			problemSet = (ProblemSet) Class
					.forName("etmo.problems.WCCI2020.MATP" + pCase)
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

//				System.out.println("RunID\t" + "IGD for " + problemSet.get(tsk).getName() + " for " + times + " times.");

				pf[tsk] = "resources/PF/StaticPF/" + problemSet.get(tsk).getHType() + "_" + problemSet.get(tsk).getNumberOfObjectives() + "D.pf";

				for (int t = 1; t <= times; t++){
					algorithm = new MOEAD(pS);

					algorithm.setInputParameter("populationSize", 100);
					algorithm.setInputParameter("maxEvaluations", 1000 * 100);

					algorithm.setInputParameter("dataDirectory", "D:\\_r\\EA\\ETMO\\MTO-cec2021-\\resources\\weightVectorFiles\\moead");

					algorithm.setInputParameter("T", 20);
					algorithm.setInputParameter("delta", 0.9);
					algorithm.setInputParameter("nr", 2);

					parameters = new HashMap();
					parameters.put("CR", 1.0);
					parameters.put("F", 0.5);
					crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover",parameters);

					// Mutation operator
					parameters = new HashMap();
					parameters.put("probability", 1.0 / pS.get(0).getNumberOfVariables());
					parameters.put("distributionIndex", 20.0);
					mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

					algorithm.addOperator("crossover", crossover);
					algorithm.addOperator("mutation", mutation);

					QualityIndicator indicator = new QualityIndicator(problemSet.get(tsk), pf[tsk]);

					SolutionSet population = algorithm.execute();
//					population.printObjectivesToFile("MOEAD_"+problemSet.get(tsk).getNumberOfObjectives()+"Obj_"+
//							problemSet.get(tsk).getName()+ "_" + problemSet.get(tsk).getNumberOfVariables() + "D_run"+t+".txt");
					double igd =  indicator.getIGD(population);
					ave[tsk] += igd;
					igds[tsk][t-1] = igd;
				}
				System.out.println("T" + (tsk+1) + "\t" + form.format(ave[tsk] / times));
//				System.out.println("Average IGD for " + problemSet.get(tsk).getName() + ": " + form.format(ave[tsk] / times));
			}
			LogIGD.LogIGD("MOEAD" + "_" + benchmark_name, pCase, igds);
			System.out.println();
			// for briefly summarization
//			for(int i=0;i<taskNum;i++) {
//				System.out.println("Average IGD for " + problemSet.get(i).getName() + ": " + form.format(ave[i] / times));
//			}
		}
	} // main
} // MOEAD_main
