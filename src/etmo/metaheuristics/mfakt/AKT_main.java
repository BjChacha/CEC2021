package etmo.metaheuristics.mfakt;

import etmo.core.*;
import etmo.operators.adaCrossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.comparators.LocationComparator;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;


public class AKT_main {
	public static void main(String args[]) throws IOException, JMException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
		ProblemSet problemSet; // The problem to solve
		Algorithm algorithm; // The algorithm to use

		AdaOperator crossover; // Crossover operator
		AdaOperator crossover1;
		Operator mutation; // Mutation operator
		Operator selection;

		HashMap parameters; // Operator parameters	

		int problemStart = 25;
		int problemEnd = 32;

		int times = 21;

		DecimalFormat form = new DecimalFormat("#.####E0");

		System.out.println("Algo: MFEA_AKT.");

		for (int pCase = problemStart; pCase <= problemEnd; pCase++) {
			problemSet = (ProblemSet) Class
					.forName("etmo.problems.benchmarks_ETMO.ETMOF" + pCase)
					.getMethod("getProblem")
					.invoke(null, null);

			int taskNum = problemSet.size();
			double[] ave = new double[taskNum];

			String[] pf = new String[taskNum];
			for (int k = 0; k < pf.length; k++)
				pf[k] = "PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";

			String pSName = problemSet.get(0).getName();
			pSName = pSName.substring(0, pSName.length() - 2);
			System.out.println(pSName + "\ttaskNum = " + taskNum + "\tfor " + times + " times.");

			algorithm = new MFEAAKT(problemSet);

			algorithm.setInputParameter("populationSize", 100 * taskNum);
			algorithm.setInputParameter("maxEvaluations", 100 * taskNum * 1000);

			algorithm.setInputParameter("rmp", 0.9);

			parameters = new HashMap();
			parameters.put("probability", 0.9);
			parameters.put("distributionIndex", 20.0);
			crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
			crossover1 = CrossoverFactory.getCrossoverOperator("NewCrossover", parameters);

			parameters = new HashMap();
			parameters.put("probability", 1.0 / problemSet.getMaxDimension());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

			parameters = new HashMap();
			parameters.put("comparator", new LocationComparator());
			selection = SelectionFactory.getSelectionOperator("BinaryTournament",
					parameters);

			algorithm.addOperator("crossover1", crossover1);
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);
			algorithm.addOperator("selection", selection);

			for (int t = 0; t <= times; t++) {
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

				double igd;
				for (int i = 0; i < taskNum; i++) {
					QualityIndicator indicator = new QualityIndicator(problemSet.get(i), pf[i]);
					if (resPopulation[i].size() == 0)
						continue;
					igd = indicator.getIGD(resPopulation[i]);
					ave[i] += igd;
				}
			}
			for (int i = 0; i < taskNum; i++)
				System.out.println("T" + (i + 1) + "\t" + form.format(ave[i] / times));
			System.out.println();
		}
	}
}
