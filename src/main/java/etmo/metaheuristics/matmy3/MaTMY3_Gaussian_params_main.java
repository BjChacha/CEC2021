package etmo.metaheuristics.matmy3;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import etmo.core.MtoAlgorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.SolutionSet;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.logging.LogIGD;

public class MaTMY3_Gaussian_params_main {
    // CONFIG
    static final Class<?> ALGORITHM_CLAZZ = MaTMY3_Gaussian.class;
    static final int MAX_POPULATION_SIZE = 100;
    static final int MAX_EVALUATION_PER_INDIVIDUAL = 1000;
    static final String CROSSOVER_TYPE = "SBX";
    static final String TRANSFER_CROSSOVER_TYPE = "SBX";
    static final double DE_CR = 0.6;
    static final double DE_F = 0.5;
    static final boolean IS_MUTATE = false;
    static final double TRANSFER_PROBABILITY = 0.5;
    static final double MUTATION_PROBABILITY = 0.5;

    static final int PLOT_TASK_ID = 31;

    static final Benchmark BENCHMARK_TYPE = Benchmark.WCCI2020;
    static final int PROBLEM_START = 1;
    static final int PROBLEM_END = 10;
    static final int PROBLEM_REPEAT_TIME = 10;

    static final boolean IGD_LOG = true;
    static final boolean IGD_PRINT = false;
    static final boolean PLOTTING = false;

    // static final String ALGO_NAME = "MaTMY3_SBX_PM0.5_Gau0.5_(0.5CMD(noSIM)_0.5EliteClosest)";
    static final String ALGO_NAME = "par(%s)_MaTMY3_SBX_PM%f_Gau%f_(0.5CMD(noSIM)_0.5EliteClosest)";

    enum Benchmark { CEC2021, CEC2017, WCCI2020 }

    public static void main(String[] args) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, InstantiationException, JMException {
        ProblemSet problemSet;

        Benchmark benchmarkName = BENCHMARK_TYPE;
        int problemStart = PROBLEM_START;
        int problemEnd = PROBLEM_END;
        int times = PROBLEM_REPEAT_TIME;

        long startTime = System.currentTimeMillis();

        HashMap<String, Double> params = new HashMap<>();
        for (int par = 0; par <= 9; par++) {
            params.clear();
            params.put("elitePartition", (par + 1) * 0.1);
            params.put("mutationProbability", 0.5);

            String fileName = String.format(ALGO_NAME, "MP", params.get("mutationProbability"), params.get("elitePartition")) + "_x" + times + "_" + benchmarkName;
            System.out.println("\nExperiment started -> " + fileName);

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
                    algorithms.add(algorithmGenerate(ALGORITHM_CLAZZ, problemSet, params));
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
                if (IGD_LOG || IGD_PRINT) {
                    int t = 0;
                    for (SolutionSet[] pop: populations){
                        double igd;
                        for (int k = 0; k < taskNum; k++) {
                            igd = indicators.get(k).getIGD(pop[k], k);
                            igds[k][t] = igd;
                        }
                        t ++;
                    }
                }
                if (IGD_LOG) {
                    LogIGD.LogIGD(fileName, problemID, igds);
                }
                if (IGD_PRINT) {
                    double[] igdMean = new double[taskNum];
                    // System.out.println("Subproblem " + problemID + ": ");
                    for (int k = 0; k < taskNum; k++) {
                        igdMean[k] = Arrays.stream(igds[k]).sum() / times;
                        System.out.println(igdMean[k]);
                    }
                }
            }
        }
        long endTime = System.currentTimeMillis();
		System.out.println("Total time cost: " + (endTime - startTime) / 1000 + " s.");

    }

    public static MtoAlgorithm algorithmGenerate(Class<?> algorithmClass, ProblemSet problemSet, HashMap<String, Double> params) throws JMException, InstantiationException, IllegalAccessException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException {
		MtoAlgorithm algorithm;
        Operator SBXCrossover;
        Operator DECrossover;
        Operator BLXAlphaCrossover;
		Operator mutation;
		HashMap<String, Double> parameters;

        algorithm = (MtoAlgorithm) algorithmClass
                    .getDeclaredConstructor(ProblemSet.class)
                    .newInstance(problemSet);

		algorithm.setInputParameter("populationSize", MAX_POPULATION_SIZE);
		algorithm.setInputParameter("maxEvaluations", MAX_EVALUATION_PER_INDIVIDUAL * problemSet.size() * MAX_POPULATION_SIZE);
        algorithm.setInputParameter("XType", CROSSOVER_TYPE);
        algorithm.setInputParameter("TXType", TRANSFER_CROSSOVER_TYPE);
        algorithm.setInputParameter("isPlot", PLOTTING);
        algorithm.setInputParameter("isMutate", IS_MUTATE);
        algorithm.setInputParameter("transferProbability", TRANSFER_PROBABILITY);
        algorithm.setInputParameter("mutationProbability", params.getOrDefault("mutationProbability", 0.5));
        algorithm.setInputParameter("plotTaskID", PLOT_TASK_ID);
        algorithm.setInputParameter("elitePartition", params.getOrDefault("elitePartition", 0.5));

        parameters = new HashMap<>();
        parameters.put("probability", 1.0);
        parameters.put("distributionIndex", 20.0);
        SBXCrossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",parameters);
        
        parameters = new HashMap<>();
        parameters.put("CR", DE_CR);
        parameters.put("F", DE_F);
        DECrossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover",parameters);

        parameters = new HashMap<>();
        parameters.put("probability", 1.0);
        parameters.put("distributionIndex", 20.0);
        parameters.put("alpha", 0.1);
        BLXAlphaCrossover = CrossoverFactory.getCrossoverOperator("BLXAlphaCrossover",parameters);

		parameters = new HashMap<>();
		parameters.put("probability", 1.0 / problemSet.get(0).getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

        algorithm.addOperator("DECrossover", DECrossover);
        algorithm.addOperator("SBXCrossover", SBXCrossover);
        algorithm.addOperator("BLXAlphaCrossover", BLXAlphaCrossover);
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
}
