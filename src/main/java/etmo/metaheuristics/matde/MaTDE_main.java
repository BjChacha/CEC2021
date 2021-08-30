package etmo.metaheuristics.matde;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import etmo.core.MtoAlgorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.operators.crossover.CrossoverFactory;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.problems.CEC2017.*;
import etmo.util.logging.LogIGD;

public class MaTDE_main {
    public static void main(String[] args) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, JMException, IOException, InstantiationException {
        int problemStart = 25;
        int problemEnd = 32;

        int times = 21;

        String benchmark_name;

        System.out.println("Algo: MaTDE.");

        long startTime = System.currentTimeMillis();

        for (int pCase = problemStart; pCase <= problemEnd; pCase++){
            ProblemSet problemSet;
          // CEC 2021
          benchmark_name = "CEC2021";
          problemSet = (ProblemSet) Class
                  .forName("etmo.problems.CEC2021.ETMOF" + pCase)
                  .getMethod("getProblem")
                  .invoke(null, null);

            // // WCCI 2020
            // benchmark_name = "WCCI2020";
            // problemSet = (ProblemSet) Class
            //         .forName("etmo.problems.WCCI2020.MATP" + pCase)
            //         .getMethod("getProblem")
            //         .invoke(null, null);

            //  // CEC2017
            //  benchmark_name = "CEC2017";
            //  ProblemSet[] cec2017 = {
            //      CIHS.getProblem(),
            //      CIMS.getProblem(),
            //      CILS.getProblem(),
            //      PIHS.getProblem(),
            //      PIMS.getProblem(),
            //      PILS.getProblem(),
            //      NIHS.getProblem(),
            //      NIMS.getProblem(),
            //      NILS.getProblem()
            //  };
            //  problemSet = cec2017[pCase-1];


            int taskNum = problemSet.size();

            List<QualityIndicator> indicators = new ArrayList<>(taskNum);
            String[] pf = new String[taskNum];
            for (int k = 0; k < pf.length; k++){
                pf[k] = "resources/PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
                indicators.add(new QualityIndicator(problemSet.get(k), pf[k]));
            }

            String pSName = problemSet.get(0).getName();
            System.out.println(pSName + "\ttaskNum = "+taskNum+"\tfor "+times+" times.");

            List<MtoAlgorithm> algorithms = new ArrayList<>(times);
			List<SolutionSet[]> populations = new ArrayList<>(times);
            
			// 初始化算法
			for (int t = 0; t < times; t++){
				algorithms.add(algorithmGenerate(problemSet));
			}
			// 并行执行times个算法
			algorithms.parallelStream().forEach(a -> {
				try {
					populations.add(a.execute());
				} catch (JMException | ClassNotFoundException e) {
					e.printStackTrace();
				}
			});

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
            LogIGD.LogIGD("MaTDE_" + "x" + times + "_" + benchmark_name, pCase, igds);
        }
        long endTime = System.currentTimeMillis();
		System.out.println("Total time cost: " + (endTime - startTime) / 1000 + " s.");

    }

    public static MtoAlgorithm algorithmGenerate(ProblemSet problemSet) throws JMException, InstantiationException, IllegalAccessException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException {
        MtoAlgorithm algorithm;
        Operator crossover1;
        Operator crossover2;

        HashMap parameters;

        int taskNum = problemSet.size();

        algorithm = new MaTDE(problemSet);

        algorithm.setInputParameter("populationSize", 100);
        algorithm.setInputParameter("archiveSize", 300);
        algorithm.setInputParameter("maxEvaluations", 1000 * taskNum * 100);

        // 迁移交叉概率
        algorithm.setInputParameter("alpha", 0.1);
        // 衰减系数
        algorithm.setInputParameter("ro", 0.8);
        // 衰减速率
        algorithm.setInputParameter("shrinkRate", 0.8);
        // Archive更新概率
        algorithm.setInputParameter("replaceRate", 0.2);

        parameters = new HashMap();
        // 原论文：CR = (0.1, 0.9)
        parameters.put("CR_LB", 0.1);
        parameters.put("CR_UB", 0.9);
        // 原论文：F = (0.1, 2)
        parameters.put("F_LB", 0.1);
        parameters.put("F_UB", 2.0);
        crossover1 = CrossoverFactory.getCrossoverOperator("RandomDECrossover",parameters);
        algorithm.addOperator("crossover1", crossover1);

        parameters = new HashMap();
        parameters.put("CR_LB", 0.1);
        parameters.put("CR_UB", 0.9);
        crossover2 = CrossoverFactory.getCrossoverOperator("RandomUniformCrossover", parameters);
        algorithm.addOperator("crossover2", crossover2);

        return algorithm;
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
