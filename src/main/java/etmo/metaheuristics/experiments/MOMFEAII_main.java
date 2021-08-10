package etmo.metaheuristics.experiments;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.comparators.LocationComparator;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;

public class MOMFEAII_main {
	public static void main(String args[]) throws IOException, JMException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
		ProblemSet problemSet; // The problem to solve
		MtoAlgorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection;

		HashMap parameters; // Operator parameters

		int taskStart = 25;
		int taskEnd = 25;

		int times = 1;

		DecimalFormat form = new DecimalFormat("#.####E0");

		System.out.println("Algo: MOMFEAII.");

		for (int pCase = taskStart; pCase <= taskEnd; pCase++ ){
			problemSet = (ProblemSet) Class
					.forName("etmo.problems.CEC2021.ETMOF" + pCase)
					.getMethod("getProblem")
					.invoke(null, null);

			int taskNum = problemSet.size();
			double ave[] = new double[taskNum];

			String[] pf = new String[taskNum];
			for (int i = 0; i < pf.length; i++){
				pf[i] = "resources/PF/StaticPF/" + problemSet.get(i).getHType() + "_" + problemSet.get(i).getNumberOfObjectives() + "D.pf";
			}

			String pSName = problemSet.get(0).getName();
			pSName = pSName.substring(0, pSName.length()-2);
			System.out.println(pSName + "\ttaskNum = "+taskNum+"\tfor "+times+" times.");

			algorithm = new MOMFEAII(problemSet);

			algorithm.setInputParameter("populationSize",100);
			algorithm.setInputParameter("maxEvaluations",1000 * taskNum * 100);
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

			double[][] igds = new double[times][taskNum];
			for (int t = 0; t < times; t++) {
				SolutionSet population[] = algorithm.execute();
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
					resPopulation[k].printObjectivesToFile("MOMFEA_" + problemSet.get(k).getNumberOfObjectives() + "Obj_" +
							problemSet.get(k).getName() + "_" + problemSet.get(k).getNumberOfVariables() + "D_run_" + t + ".txt");

				}
				double igd;
				for (int k = 0; k < taskNum; k++){
					QualityIndicator indicator = new QualityIndicator(problemSet.get(k), pf[k]);
					if (population[k].size() == 0)
						continue;

					igd = indicator.getIGD(resPopulation[k]);
//                    ave[k] += igd;
					// DEBUG
					igds[t][k] = igd;
				}
//                // DEBUG
//                LogIGD.LogIGD("MaTBML_" + problemSet.get(0).getName() + "D_run_" + t + ".txt", igds[t]);
			}
			for(int i = 0; i < taskNum; i++) {
				double[] tmp = new double[times];
				for (int t = 0; t < times; t++){
					tmp[t] = igds[t][i];
				}
//                Arrays.sort(tmp);
				double best, worst, mean, median, std = 0;
//                best = tmp[0];
//                worst = tmp[times-1];
				mean = Arrays.stream(tmp).sum() / times;
//                median = tmp[times/2];
//                for (double e: tmp){
//                    std += Math.pow(e - mean, 2);
//                }
//                std = Math.sqrt(std / （times - 1));
//                System.out.println(best + "\t" + worst + "\t" + mean + "\t" + median + "\t" + std);
				System.out.println(form.format(mean));
			}
			System.out.println();
		}
	}
}
