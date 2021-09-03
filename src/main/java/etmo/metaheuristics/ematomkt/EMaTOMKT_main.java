package etmo.metaheuristics.ematomkt;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import etmo.core.MtoAlgorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.logging.LogIGD;

public class EMaTOMKT_main {
    static int MAX_POPULATION_SIZE = 100;
    static int MAX_EVALUATION_PER_INDIVIDUAL = 1000;
    static boolean LOG_IGD = false;

    enum Benchmark { CEC2021, CEC2017, WCCI2020; }

    public static void main(String[] args) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, InstantiationException, JMException {
        ProblemSet problemSet;
        Benchmark benchmarkName = Benchmark.WCCI2020;
        int problemStart = 1;
        int problemEnd = 1;
        int times = 1;
        String fileName = "EMaTOMKT_x" + times + "_" + benchmarkName;

        System.out.println("\nExperiment started -> " + fileName);

        long startTime = System.currentTimeMillis();

        for (int problemID = problemStart; problemID <= problemEnd; problemID ++) {
            problemSet = getProblemSet(benchmarkName, problemID);
            int taskNum = problemSet.size();
            String[] pf = new String[taskNum];
            List<QualityIndicator> indicators = new ArrayList<>(taskNum);
            for (int k = 0; k < taskNum; k++) {
				pf[k] = "resources/PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
				indicators.add(new QualityIndicator(problemSet.get(k), pf[k]));
            }

            String pSName = problemSet.get(0).getName();
			System.out.println(pSName + "\ttaskNum = " + taskNum + "\tfor " + times + " times.");

			List<MtoAlgorithm> algorithms = new ArrayList<>(times);
			List<SolutionSet[]> populations = new ArrayList<>(times);

            			// 初始化算法
			for (int t = 0; t < times; t++){
				algorithms.add(algorithmGenerate(problemSet));
			}

			// // 并行执行times个算法
			// algorithms.parallelStream().forEach(a -> {
			// 	try {
			// 		populations.add(a.execute());
			// 	} catch (JMException | ClassNotFoundException e) {
			// 		e.printStackTrace();
			// 	}
			// });

			// 串行执行
			for (int t = 0; t < times; t++){
				populations.add(algorithms.get(t).execute());
			}

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

            if (LOG_IGD) {
				LogIGD.LogIGD(fileName, problemID, igds);
			}
        }

        long endTime = System.currentTimeMillis();
		System.out.println("Total time cost: " + (endTime - startTime) / 1000 + " s.");

    }

    public static MtoAlgorithm algorithmGenerate(ProblemSet problemSet) throws JMException, InstantiationException, IllegalAccessException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException {
		MtoAlgorithm algorithm;
		Operator crossover;
		Operator mutation;
		HashMap<String, Double> parameters;

		algorithm = new EMaTOMKT(problemSet);

		algorithm.setInputParameter("populationSize", MAX_POPULATION_SIZE);
		algorithm.setInputParameter("maxEvaluations", MAX_EVALUATION_PER_INDIVIDUAL * problemSet.size() * MAX_POPULATION_SIZE);

		parameters = new HashMap();
		parameters.put("probability", 1.0);
		parameters.put("distributionIndex", 20.0);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",parameters);

		parameters = new HashMap();
		parameters.put("probability", 1.0 / problemSet.get(0).getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);

		return algorithm;
	}

    public static ProblemSet getProblemSet(Benchmark problemName, int problemId) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
		ProblemSet problemSet;
        switch (problemName) {
            case CEC2021:
                problemSet = (ProblemSet) Class
                .forName("etmo.problems.CEC2021.ETMOF" + problemId)
                .getMethod("getProblem")
                .invoke(null, null);
                break;
            case CEC2017:
                String[] problemNames = {
                    "CIHS", "CIMS", "CILS",
                    "PIHS", "PIMS", "PILS",
                    "NIHS", "NIMS", "NILS"
                };
                problemSet = (ProblemSet) Class
                    .forName("etmo.problems.CEC2017." + problemNames[problemId-1])
                    .getMethod("getProblem")
                    .invoke(null, null);
                break;
            case WCCI2020:
                problemSet = (ProblemSet) Class
                        .forName("etmo.problems.WCCI2020.MATP" + problemId)
                        .getMethod("getProblem")
                        .invoke(null, null);
                break;
            default:
                System.out.println("Error: unknown benchmark type: " + problemName);
                problemSet = (ProblemSet) Class
                        .forName("etmo.problems.CEC2021.ETMOF" + problemId)
                        .getMethod("getProblem")
                        .invoke(null, null);
        }
        return problemSet;
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