package etmo.metaheuristics.matbml2;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.logging.LogIGD;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class MOEAD_T_DRA_main {
	static int MAX_POPULATION_SIZE = 100;
	static int MAX_EVALUATION_PER_INDIVIDUAL = 1000;
	static boolean LOG_IGD = true;

	public static MtoAlgorithm algorithmGenerate(Class algorithmClass, ProblemSet problemSet) throws JMException, InstantiationException, IllegalAccessException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException {
		MtoAlgorithm algorithm;
		Operator crossover;
		Operator crossover2;
		Operator mutation;
		HashMap<String, Double> parameters;

//		algorithm = (MtoAlgorithm) algorithmClass.getDeclaredConstructor().newInstance();
		algorithm = new MOEAD_T_DRA(problemSet);

		algorithm.setInputParameter("populationSize", MAX_POPULATION_SIZE);
		algorithm.setInputParameter("maxEvaluations", MAX_EVALUATION_PER_INDIVIDUAL * problemSet.size() * MAX_POPULATION_SIZE);

		// TODO: 改成相对路径
		algorithm.setInputParameter("dataDirectory", "resources/weightVectorFiles/moead");

		algorithm.setInputParameter("T", 20);
		algorithm.setInputParameter("delta", 0.9);
		algorithm.setInputParameter("nr", 2);

		algorithm.setInputParameter("aStep", 1);
		algorithm.setInputParameter("transferP", 1.0);

		parameters = new HashMap();
		parameters.put("CR", 0.6);
		parameters.put("F", 0.5);
		crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover",parameters);

//		parameters = new HashMap();
//		parameters.put("CR_LB", 0.1);
//		parameters.put("CR_UB", 0.9);
//		parameters.put("F_LB", 0.1);
//		parameters.put("F_UB", 2.0);
//		crossover = CrossoverFactory.getCrossoverOperator("RandomDECrossover",parameters);

		parameters = new HashMap();
		parameters.put("probability", 0.9);
		parameters.put("distributionIndex", 20.0);
		crossover2 = CrossoverFactory.getCrossoverOperator("SBXCrossover",parameters);

//		parameters = new HashMap();
//		parameters.put("CR", 0.5);
//		crossover2 = CrossoverFactory.getCrossoverOperator("UniformCrossover", parameters);

//		parameters = new HashMap();
//		crossover2 = CrossoverFactory.getCrossoverOperator("TransferDECrossover", parameters);

		parameters = new HashMap();
		parameters.put("probability", 1.0 / problemSet.get(0).getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("crossover2", crossover2);
		algorithm.addOperator("mutation", mutation);

		return algorithm;
	}

	public static void main(String[] args) throws JMException, SecurityException, IOException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, InstantiationException {
		ProblemSet problemSet; // The problem to solve

		String benchmarkName = "CEC2021";
		Class algorithmClass = MOEAD_T.class;
		String algorithmName = algorithmClass.getName();
		int taskStart = 25;
		int taskEnd = 32;
		int times = 10;

		// System.out.println("Algo:" + algorithmName + ".");

		System.out.println();
		String fileName = "MOEAD_T(rnd(1.0)_DE(CR0.6)_SBX_A(1)_RA(3)_T(type2)" + "_x" + times + "_" + benchmarkName;
		System.out.println("Experiment started -> " + fileName);

		long startTime = System.currentTimeMillis();

		for (int pCase = taskStart; pCase <= taskEnd; pCase++){
			problemSet = getProblemSet(benchmarkName, pCase);
			int taskNum = problemSet.size();
			String[] pf = new String[taskNum];
			List<QualityIndicator> indicators = new ArrayList<>(taskNum);
			for (int k = 0; k < taskNum; k++){
				pf[k] = "resources/PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
				indicators.add(new QualityIndicator(problemSet.get(k), pf[k]));
 			}

			String pSName = problemSet.get(0).getName();
			System.out.println(pSName + "\ttaskNum = " + taskNum + "\tfor " + times + " times.");

			List<MtoAlgorithm> algorithms = new ArrayList<>(times);
			List<SolutionSet[]> populations = new ArrayList<>(times);

			// 初始化算法
			for (int t = 0; t < times; t++){
				algorithms.add(algorithmGenerate(algorithmClass, problemSet));
			}
			// 并行执行times个算法
			algorithms.parallelStream().forEach(a -> {
				try {
					populations.add(a.execute());
				} catch (JMException | ClassNotFoundException e) {
					e.printStackTrace();
				}
			});

//			// 串行执行
//			for (int t = 0; t < times; t++){
//				populations.add(algorithms.get(t).execute());
//			}

			// 计算IGD
			double[][] igds = new double[taskNum][times];
			int t = 0;
			for (SolutionSet[] pop: populations){
				SolutionSet[] resPopulation = getEvalPopulations(pop, problemSet);
				double igd;
				for (int k = 0; k < taskNum; k++) {
					igd = indicators.get(k).getIGD(resPopulation[k]);
					igds[k][t] = igd;
				}
				t ++;
			}

//			double[][] igds = new double[taskNum][times];
//			for (int t = 0; t < times; t++) {
//				algorithm = algorithmGenerate(problemSet);
//
//				SolutionSet[] population = algorithm.execute();
//				SolutionSet[] resPopulation = getEvalPopulations(population, problemSet);
//
//				double igd;
//				for (int k = 0; k < taskNum; k++) {
//					igd = indicators[k].getIGD(resPopulation[k]);
//					igds[k][t] = igd;
//				}
//			}

			if (LOG_IGD) {
				LogIGD.LogIGD(fileName, pCase, igds);
//				LogIGD.LogIGD("MOEAD" + "_x" + times + "_" + benchmarkName, pCase, igds);
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("Total time cost: " + (endTime - startTime) / 1000 + " s.");
	}

	public static ProblemSet getProblemSet(String problemName, int problemId) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
		ProblemSet ps;
		if (problemName.equalsIgnoreCase("CEC2021")){
			ps = (ProblemSet) Class
					.forName("etmo.problems.CEC2021.ETMOF" + problemId)
					.getMethod("getProblem")
					.invoke(null, null);
		} else if (problemName.equalsIgnoreCase("WCCI2020")){
			ps = (ProblemSet) Class
					.forName("etmo.problems.WCCI2020.MATP" + problemId)
					.getMethod("getProblem")
					.invoke(null, null);
		} else if (problemName.equalsIgnoreCase("CEC2017")){
			String[] problemNames = {
					"CIHS", "CIMS", "CILS",
					"PIHS", "PIMS", "PILS",
					"NIHS", "NIMS", "NILS"
			};
			ps = (ProblemSet) Class
					.forName("etmo.problems.CEC2017." + problemNames[problemId-1])
					.getMethod("getProblem")
					.invoke(null, null);
		} else {
			System.out.println("Error: unknown benchmark type: " + problemName);
			ps = (ProblemSet) Class
					.forName("etmo.problems.CEC2021.ETMOF" + problemId)
					.getMethod("getProblem")
					.invoke(null, null);
		}

		return ps;
	}

	public static SolutionSet[] getEvalPopulations(SolutionSet[] population, ProblemSet problemSet){
		int taskNum = population.length;
		SolutionSet[] resPopulation = new SolutionSet[taskNum];

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
		}
		return resPopulation;
	}
}
