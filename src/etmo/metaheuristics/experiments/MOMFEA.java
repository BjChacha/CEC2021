package etmo.metaheuristics.experiments;

import etmo.core.*;
import etmo.metaheuristics.matbml.Utils;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PORanking;
import etmo.util.PseudoRandom;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.LocationComparator;
import etmo.util.sorting.SortingIdx;

import java.util.Arrays;
import java.util.stream.IntStream;

public class MOMFEA extends MtoAlgorithm {

	private int populationSize;
	
	private SolutionSet[] population;
	private SolutionSet[] offspringPopulation;
	private SolutionSet union;
	
	int evaluations;
	int maxEvaluations;
	
	Operator crossover;
	Operator mutation;
	Operator selection;
	
	double rmp;

	int taskNum;
	int[][] objPos;
	double[] ideals;
	double[] nadirs;

	int leaderNum;
	int[] groups;
	int[] leaders;
	double[][] distances;
	double[][] scores;

	int k1 = 1;
	int k2 = 3;

	double P = 0.1;

	// DEBUG
	int[] runTimes;
	int savedTimes;
	int BTimes;
	int NTimes;
	int WTimes;

	Distance distance = new Distance();
	
	public MOMFEA(ProblemSet problemSet) {
		super(problemSet);
	}

	@Override
	public SolutionSet[] execute() throws JMException, ClassNotFoundException {
		initState();
		initPopulation();

		savedTimes = 0;
		evaluations = 0;
		while (evaluations < maxEvaluations) {
			solelyConvergence(k1);
			transferConvergence(k2);
		}
		System.out.println("Saved times: " + savedTimes);
		return population;
	}

	void initPopulation() throws JMException, ClassNotFoundException {
		for (int k = 0; k < taskNum; k++) {
			population[k] = new SolutionSet(populationSize);
			offspringPopulation[k] = new SolutionSet(populationSize);
			for (int i = 0; i < populationSize; i++) {
				Solution newSolution = new Solution(problemSet_);
				int id = k;
				problemSet_.get(id).evaluate(newSolution);
				problemSet_.get(id).evaluateConstraints(newSolution);
				evaluations++;

				newSolution.setSkillFactor(id);
				population[k].add(newSolution);

			} // for
			assignFitness(population[k], k);
		}
		updateINPoint();
	} // initPopulation

	void initState(){
		populationSize = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();
		rmp =  ((Double) getInputParameter("rmp")).doubleValue();

		crossover = operators_.get("crossover");
		mutation = operators_.get("mutation");
		selection = operators_.get("selection");

		taskNum = problemSet_.size();
		population = new SolutionSet[taskNum];
		offspringPopulation = new SolutionSet[taskNum];

		objPos = new int[taskNum][2];

		ideals = new double[problemSet_.getTotalNumberOfObjs()];
		nadirs = new double[problemSet_.getTotalNumberOfObjs()];

		leaderNum = (int) Math.sqrt(taskNum);
		groups = new int[taskNum];
		leaders = new int[leaderNum];
		distances = new double[taskNum][taskNum];
		scores = new double[taskNum][taskNum];

		Arrays.fill(ideals, Double.POSITIVE_INFINITY);
		for (int k = 0; k < taskNum; k++) {
			objPos[k][0] = problemSet_.get(k).getStartObjPos();
			objPos[k][1] = problemSet_.get(k).getEndObjPos();
			Arrays.fill(scores[k], 3);
		}

		// DEBUG
		runTimes = new int[taskNum];
		savedTimes = 0;
		BTimes = 0;
		NTimes = 0;
		WTimes = 0;
	}

	void solelyConvergence(int times) throws JMException {
		double[] oldIdeal = ideals.clone();
		for (int t = 0; t < times; t++){
			for (int k = 0; k < taskNum; k++){
				if (runTimes[k] > populationSize * 1000)
					continue;

				createOffspring(k);
				getNextPopulation(k);
				runTimes[k] += populationSize;
			}
		}
		updateINPoint();

		Arrays.fill(leaders, -1);
		Arrays.fill(groups, -1);
		double[] improvements = new double[taskNum];
		for (int k = 0; k < taskNum; k++){
			for (int j = objPos[k][0]; j <= objPos[k][1]; j++){
				improvements[k] += (oldIdeal[j] - ideals[j]);
			}
		}

		if (Arrays.stream(improvements).sum() > 0){
			int[] idx = SortingIdx.SortingIdx(improvements, true);
			for (int i = 0; i < leaderNum; i++)
				leaders[i] = idx[i];

			UpdateDistances();
			for (int k = 0; k < taskNum; k++) {
				int finalK = k;
				if (IntStream.of(leaders).anyMatch(x -> x == finalK))
					continue;
				double[] finalScore = new double[leaders.length];
				for (int i = 0; i < finalScore.length; i++) {
					// TODO: 调整计算公式

//                double factor = scores[k][leaders_[i]] > 0 ? scores[k][leaders_[i]] : 0.5 * Math.pow(2, scores[k][leaders_[i]]);
//                finalScore[i] = factor * Math.exp(1 / (1 + distances[k][leaders_[i]]) - 1);
					finalScore[i] = scores[k][leaders[i]] > 0 ? 0.5 * scores[k][leaders[i]] * Math.exp(1 / (1 + distances[k][leaders[i]]) - 1) : 0;
//                finalScore[i] = scores[k][leaders_[i]] * Math.exp(1 / (1 + distances[k][leaders_[i]]) - 1);
				}

				if (Arrays.stream(finalScore).sum() == 0)
					continue;
				else {
					groups[k] = leaders[Utils.rouletteExceptZero(finalScore)];
//				groups[k] = leaders[Utils.roulette(finalScore)];
				}
			}
		}else{
			for (int k = 0; k < taskNum; k++){
				if (PseudoRandom.randDouble() < P){
					int idx = k;
					while (idx == k)
						idx = PseudoRandom.randInt(0, taskNum - 1);
					groups[k] = idx;
				}
			}
		}
	}

	void transferConvergence(int times) throws JMException {
		for (int t = 0; t < times; t++) {
			for (int leader : leaders) {
				if (leader < 0 || runTimes[leader] > populationSize * 1000)
					continue;

				createOffspring(leader);
				getNextPopulation(leader);
				runTimes[leader] += populationSize;
			}
		}

		for (int k = 0; k < taskNum; k++){
			if (groups[k] < 0)
				continue;

			// 用leader种群来评价本任务
			double[] tmpIdeal = ideals.clone();
			double[] tmpNadir = nadirs.clone();
			SolutionSet tmpSet = new SolutionSet(populationSize);
			int leader = groups[k];
			for (int i = 0; i < population[leader].size(); i++){
				Solution tmp = new Solution(population[leader].get(i));
				tmp.setSkillFactor(k);
				problemSet_.get(k).evaluate(tmp);
				problemSet_.get(k).evaluateConstraints(tmp);
				evaluations ++;
				tmpSet.add(tmp);
				for (int j = objPos[k][0]; j <= objPos[k][1]; j++){
					tmpIdeal[j] = Math.min(tmpIdeal[j], tmp.getObjective(j));
					tmpNadir[j] = Math.max(tmpNadir[j], tmp.getObjective(j));
				}
			}

			// 判断leader种群对本任务是否有帮助
			boolean better = false;
			boolean worse = false;
			for (int j = objPos[k][0]; j <= objPos[k][1]; j++){
				if (tmpIdeal[j] < ideals[j])
					better = true;
				if (tmpNadir[j] > nadirs[j])
					worse = true;
			}

			// 根据情况应用leader种群
			if (better) {
				// Union and selection
//				unionAndSelection(population[k], tmpSet);

				createTransferOffspring(k, tmpSet);
				getNextPopulation(k);

				scores[k][leader] += 1;

				BTimes ++;
				savedTimes += k2 * populationSize;
			}
			else if (worse){
				scores[k][leader] = 0;

				WTimes ++;
				savedTimes -= populationSize;
			}
			else{
				scores[k][leader] -= 1;

				NTimes ++;
				savedTimes -= populationSize;
			}
		}

		updateINPoint();
	}

	void createOffspring(int task) throws JMException {
		offspringPopulation[task].clear();
		Solution[] parents = new Solution[2];
		for (int i = 0; i < (populationSize / 2); i++){
			parents[0] = (Solution) selection.execute(population[task]);
			parents[1] = (Solution) selection.execute(population[task]);
			Solution[] offSpring = (Solution[]) crossover.execute(parents);
			mutation.execute(offSpring[0]);
			mutation.execute(offSpring[1]);
			resetObjectives(offSpring[0]);
			resetObjectives(offSpring[1]);
			problemSet_.get(task).evaluate(offSpring[0]);
			problemSet_.get(task).evaluate(offSpring[1]);
			problemSet_.get(task).evaluateConstraints(offSpring[0]);
			problemSet_.get(task).evaluateConstraints(offSpring[1]);
			offspringPopulation[task].add(offSpring[0]);
			offspringPopulation[task].add(offSpring[1]);
			evaluations += 2;
		}
	}

	void createTransferOffspring(int task1, SolutionSet assist) throws JMException {
		// task1 -> task2
		offspringPopulation[task1].clear();
		Solution[] parents = new Solution[2];
		for (int i = 0; i < (populationSize / 2); i++){
			parents[0] = (Solution) selection.execute(population[task1]);
			parents[1] = (Solution) selection.execute(assist);
			Solution[] offSpring = (Solution[]) crossover.execute(parents);
			mutation.execute(offSpring[0]);
			mutation.execute(offSpring[1]);
			offSpring[0].setSkillFactor(task1);
			offSpring[1].setSkillFactor(task1);
			resetObjectives(offSpring[0]);
			resetObjectives(offSpring[1]);
			problemSet_.get(task1).evaluate(offSpring[0]);
			problemSet_.get(task1).evaluate(offSpring[1]);
			problemSet_.get(task1).evaluateConstraints(offSpring[0]);
			problemSet_.get(task1).evaluateConstraints(offSpring[1]);
			offspringPopulation[task1].add(offSpring[0]);
			offspringPopulation[task1].add(offSpring[1]);
			evaluations += 2;
		}
	}

	void getNextPopulation(int task) {
//		union = population[task].union(offspringPopulation[task]);
//
//		assignFitness(union);
//		union.sort(new LocationComparator());
//
//		population[task].clear();
//
//		for (int i = 0; i < populationSize; i++)
//			population[task].add(union.get(i));
		unionAndSelection(population[task], offspringPopulation[task]);
	}

	void unionAndSelection(SolutionSet P, SolutionSet Q){
		int taskId = P.get(0).getSkillFactor();
		union = P.union(Q);
		assignFitness(union, taskId);
		union.sort(new LocationComparator());

		P.clear();
		for (int i = 0; i < populationSize; i++)
			P.add(union.get(i));
	}

	void assignFitness(SolutionSet pop, int taskId) {
		for (int i = 0; i < pop.size(); i++)
			pop.get(i).setLocation(Integer.MAX_VALUE);

		rankSolutionOnTask(pop, taskId);
	}
	
	void rankSolutionOnTask(SolutionSet pop, int taskId) {
		boolean selec[] = new boolean[problemSet_.getTotalNumberOfObjs()];

		for (int i = 0; i < selec.length; i++) {
			if (i < objPos[taskId][0] || i > objPos[taskId][1])
				selec[i] = false;
			else
				selec[i] = true;
		}
		
		PORanking pr = new PORanking(pop, selec);	
		int loc = 0;
		for (int i = 0; i < pr.getNumberOfSubfronts(); i++) {
			SolutionSet front = pr.getSubfront(i);
			distance.crowdingDistanceAssignment(front, problemSet_.getTotalNumberOfObjs(), selec);
			front.sort(new CrowdingComparator());
			for (int j = 0; j < front.size(); j++) {
				if (loc < front.get(j).getLocation())
					front.get(j).setLocation(loc);
				loc++;
			}
		}
	}

	void resetObjectives(Solution sol) {
		for (int i = 0; i < sol.getNumberOfObjectives(); i++)
			sol.setObjective(i, Double.POSITIVE_INFINITY);
	}

	void updateINPoint(){
		for (int k = 0; k < taskNum; k++){
			for (int j = objPos[k][0]; j <= objPos[k][1]; j++){
				for (int i = 0; i < population[k].size(); i++){
					ideals[j] = Math.min(ideals[j], population[k].get(i).getObjective(j));
					nadirs[j] = Math.max(nadirs[j], population[k].get(i).getObjective(j));
				}
			}
		}
	}

	void UpdateDistances() throws JMException {
		// Wasserstein Distance (psedo)
		for (int i = 0; i < taskNum - 1; i++){
			for (int j = i + 1; j < taskNum; j++){
				double d1 = WassersteinDistance.getWD(population[i].getMat(), population[j].getMat());
				distances[i][j] = distances[j][i] = d1;
			}
		}

//        // KL Diversity
//        KLD kld = new KLD(problemSet_, populations_);
//        for (int i = 0; i < taskNum_; i++)
//            distances[i] = kld.getKDL(i);

//        // random
//        for (int i = 0; i < taskNum_; i++)
//            Arrays.fill(distances[i], 0);
	}
}
