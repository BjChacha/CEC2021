package etmo.metaheuristics.momfea;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.text.DecimalFormat;
import java.util.HashMap;

import etmo.problems.CEC2017.*;
import etmo.util.comparators.LocationComparator;
import etmo.core.Algorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;

import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;

import etmo.qualityIndicator.QualityIndicator;
import etmo.util.JMException;
import etmo.util.logging.LogIGD;

public class MOMFEA_main {
	static final boolean PROCESS_LOG = true;
	public static void main(String args[]) throws IOException, JMException, ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
		ProblemSet problemSet; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection;

		HashMap parameters; // Operator parameters

		int taskStart = 1;
		int taskEnd = 6;

		int times = 20;

		DecimalFormat form = new DecimalFormat("#.####E0");

        String ALGO_NAME = "MOMFEA";
        String benchmarkName = "CEC2019";
        String fileName = ALGO_NAME + "_x" + times + "_" + benchmarkName;
        String folderPath;
        File folder;
        if (PROCESS_LOG) {
            folderPath = "./data/process";
            folder = new File(folderPath);
            if (!folder.exists()){
                folder.mkdirs();
            }
            folder = new File(folderPath + "/" + fileName);
            if (!folder.exists()){
                folder.mkdirs();
            }
        }
		System.out.println("Algo: " + ALGO_NAME + ".");

		for (int pCase = taskStart; pCase <= taskEnd; pCase++ ){
            if (PROCESS_LOG) {
                folder = new File(folderPath + "/" + fileName + "/" + Integer.toString(pCase));
                if (!folder.exists()){
                    folder.mkdirs();
                }
            }

			// // CEC2019
			// problemSet = (ProblemSet) Class
			// 		.forName("etmo.problems.CEC2021.ETMOF" + pCase)
			// 		.getMethod("getProblem")
			// 		.invoke(null, null);

			// // WCCI 2020
			// problemSet = (ProblemSet) Class
			// 		.forName("etmo.problems.WCCI2020.MATP" + pCase)
			// 		.getMethod("getProblem")
			// 		.invoke(null, null);

            // CEC2019
            problemSet = (ProblemSet) Class
                    .forName("etmo.problems.CEC2019.MATP" + pCase)
                    .getMethod("getProblem")
                    .invoke(null, null);

			// // CEC2017
			// ProblemSet[] cec2017 = {
			// 		CIHS.getProblem(),
			// 		CIMS.getProblem(),
			// 		CILS.getProblem(),
			// 		PIHS.getProblem(),
			// 		PIMS.getProblem(),
			// 		PILS.getProblem(),
			// 		NIHS.getProblem(),
			// 		NIMS.getProblem(),
			// 		NILS.getProblem()
			// };
			// problemSet = cec2017[pCase-1];

			int taskNum = problemSet.size();
			double ave[] = new double[taskNum];

			String[] pf = new String[taskNum];
			for (int i = 0; i < pf.length; i++){
				pf[i] = "resources/PF/StaticPF/" + problemSet.get(i).getHType() + "_" + problemSet.get(i).getNumberOfObjectives() + "D.pf";
			}

			String pSName = problemSet.get(0).getName();
			pSName = pSName.substring(0, pSName.length()-2);
			System.out.println(pSName + "\ttaskNum = "+taskNum+"\tfor "+times+" times.");

			algorithm = new MOMFEA(problemSet);

			algorithm.setInputParameter("populationSize",100*taskNum);
			algorithm.setInputParameter("maxEvaluations",1000 * taskNum * 100);
			algorithm.setInputParameter("rmp", 0.9);
			algorithm.setInputParameter("isProcessLog", PROCESS_LOG);

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

			double[][] IGDs = new double[taskNum][times];
			for (int t = 0; t < times; t++) {
				long startTime = System.currentTimeMillis();

				SolutionSet population = algorithm.execute();

				long endTime = System.currentTimeMillis();
				System.out.println("epoch: " + t + "\trunning: " + (endTime-startTime)/1000 + " s.");

				if (PROCESS_LOG) {
                    File src = new File("./data/tmp_mfea.txt");
                    File trg = new File(folderPath + "/" + fileName + "/" + Integer.toString(pCase) + "/" + Integer.toString(t+1) + ".txt");
                    src.renameTo(trg);
                }
				
				SolutionSet[] resPopulation = new SolutionSet[taskNum];
				for (int i = 0; i < taskNum; i++)
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
				for(int i = 0; i < taskNum; i++){
					QualityIndicator indicator = new QualityIndicator(problemSet.get(i), pf[i]);
					if(resPopulation[i].size()==0)
						continue;
					igd =  indicator.getIGD(resPopulation[i]);
					IGDs[i][t] = igd;
					ave[i] += igd;
				}
			}
			LogIGD.LogIGD("MOMFEA_x" + times + "_CEC2017", pCase, IGDs);
			for(int i=0;i<taskNum;i++)
				System.out.println("T" + (i+1) + "\t" + form.format(ave[i] / times));
			System.out.println();
		}
	}
}
