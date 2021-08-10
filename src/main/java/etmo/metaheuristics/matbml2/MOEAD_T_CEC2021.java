//  MOEAD.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package etmo.metaheuristics.matbml2;

import etmo.core.*;
import etmo.util.JMException;
import etmo.util.KLD;
import etmo.util.PseudoRandom;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.StringTokenizer;
import java.util.Vector;

public class MOEAD_T_CEC2021 extends MtoAlgorithm {
	private int populationSize_;
	private SolutionSet[] population_;

	double[][] z_;
	double[][][] lambda_;
	int T_;
	int[][][] neighborhood_;
	double delta_;
	int nr_;

	int aStep;
	double[][] transferP;
	double[][] implicitP;
	int[] preBetterCount;
	double[][] scores;
	int[] preTransferTask;
	boolean[][] bannedTask;

	int[] convergeTimes;
	int[][] transferTimes;
	int[][] transferBetterTimes;

	Solution[][] indArray_;
	String functionType_;

	int evaluations_;
	int maxEvaluations_;

	Operator crossover_;
	Operator crossover2_;
	Operator mutation_;

	String dataDirectory_;

	int taskNum_;

	public MOEAD_T_CEC2021(ProblemSet problemSet) {
		super(problemSet);

		functionType_ = "_TCHE1";

	} // MOEAD

	private void initState() {
		taskNum_ = problemSet_.size();
		evaluations_ = 0;
		maxEvaluations_ = (Integer) this.getInputParameter("maxEvaluations");
		populationSize_ = (Integer) this.getInputParameter("populationSize");
		dataDirectory_ = this.getInputParameter("dataDirectory").toString();
		T_ = (Integer) this.getInputParameter("T");
		nr_ = (Integer) this.getInputParameter("nr");
		delta_ = (Double) this.getInputParameter("delta");

		aStep = (Integer) this.getInputParameter("aStep");
		double tP = (Double) this.getInputParameter("transferP");

		crossover_ = operators_.get("crossover"); // default: DE crossover
		crossover2_ = operators_.get("crossover2");
		mutation_ = operators_.get("mutation"); // default: polynomial mutation

		population_ = new SolutionSet[taskNum_];
		indArray_ = new Solution[taskNum_][];
		neighborhood_ = new int[taskNum_][populationSize_][T_];
		z_ = new double[taskNum_][];
		lambda_ = new double[taskNum_][][];

		transferP = new double[taskNum_][taskNum_];
		implicitP = new double[taskNum_][taskNum_];
		scores = new double[taskNum_][taskNum_];
		preBetterCount = new int[taskNum_];
		preTransferTask = new int[taskNum_];
		bannedTask = new boolean[taskNum_][taskNum_];

		convergeTimes = new int[taskNum_];
		transferTimes = new int[taskNum_][taskNum_];
		transferBetterTimes = new int[taskNum_][taskNum_];

		Arrays.fill(preTransferTask, -1);
		Arrays.fill(convergeTimes, 0);

		for (int k = 0; k < taskNum_; k++){
			indArray_[k] = new Solution[problemSet_.get(k).getNumberOfObjectives()];
			z_[k] = new double[problemSet_.get(k).getNumberOfObjectives()];
			lambda_[k] = new double[populationSize_][problemSet_.get(k).getNumberOfObjectives()];

			Arrays.fill(scores[k], 3);
			Arrays.fill(transferP[k], tP);
			Arrays.fill(implicitP[k], 0.5);
			Arrays.fill(bannedTask[k], false);

			Arrays.fill(transferTimes[k], 0);
			Arrays.fill(transferBetterTimes[k], 0);
		}
	}

	private Solution normalReproduce(int taskId, int currentId, int type) throws JMException {
		Solution child;

		Vector<Integer> p = new Vector<Integer>();
		matingSelection(p, currentId, 2, type, taskId);

		Solution[] parents = new Solution[3];

		parents[0] = population_[taskId].get(p.get(0));
		parents[1] = population_[taskId].get(p.get(1));
		parents[2] = population_[taskId].get(currentId);

		child = (Solution) crossover_.execute(new Object[]{
				population_[taskId].get(currentId), parents});
		mutation_.execute(child);
		return child;
	}

	private Solution transferReproduce(int targetTaskId, int sourceTaskId, int currentId) throws JMException {
		Solution child;

//		// DE
//		Solution[] parents = new Solution[3];
//		parents[0] = population_[sourceTaskId].get(currentId);
//		parents[1] = population_[targetTaskId].get(PseudoRandom.randInt(0, populationSize_ - 1));
//		parents[2] = population_[targetTaskId].get(currentId);
//		child = (Solution) crossover_.execute(new Object[]{population_[targetTaskId].get(currentId), parents});

		// SBX
		Solution[] parents = new Solution[2];
		parents[0] = population_[sourceTaskId].get(currentId);
		parents[1] = population_[targetTaskId].get(PseudoRandom.randInt(0, populationSize_ - 1));
		child = ((Solution[]) crossover2_.execute(parents))[PseudoRandom.randInt(0, 1)];

//		// TDE
//		int r1 = PseudoRandom.randInt(0, populationSize_ - 1);
//		int r2 = PseudoRandom.randInt(0, populationSize_ - 1);
//		int r3 = PseudoRandom.randInt(0, populationSize_ - 1);
//		int r4 = PseudoRandom.randInt(0, populationSize_ - 1);
//		Solution[] parents = new Solution[6];
//		parents[0] = new Solution(population_[targetTaskId].get(currentId));
//		parents[1] = new Solution(population_[sourceTaskId].get(currentId));
//		parents[2] = new Solution(population_[sourceTaskId].get(r1));
//		parents[3] = new Solution(population_[sourceTaskId].get(r2));
//		parents[4] = new Solution(population_[targetTaskId].get(r3));
//		parents[5] = new Solution(population_[targetTaskId].get(r4));
//		((TransferDECrossover) crossover2_).adaptive(population_[targetTaskId], population_[sourceTaskId]);
//		child = (Solution) crossover2_.execute(parents);

		return child;
	}

	public SolutionSet[] execute() throws JMException, ClassNotFoundException {
		initState();
		initUniformWeight();
		initNeighborhood();
		initPopulation();
		initIdealPoint();

		while (evaluations_ < maxEvaluations_) {
			iterate();
//			if (evaluations_ % (populationSize_ * 20) == 0){
//				LogPopulation.LogPopulation("MOEAD", population_, problemSet_, evaluations_, false);
//			}
		}

		return population_;
	}

	public void iterate() throws JMException {
		for (int taskId = 0; taskId < taskNum_; taskId++) {
			int assistTask = getSourceTaskId(taskId, "random");
			if (PseudoRandom.randDouble() < 1) {
				transferConverge(taskId, assistTask);
			} else {
				solelyConverge(taskId, 1);
			}
		}
	}


	public void solelyConverge(int taskId, int times) throws JMException {
		int betterCount = 0;
		convergeTimes[taskId] ++;
		for (int t = 0; t < times; t ++){
			for (int i = 0; i < populationSize_; i++) {
				int type = PseudoRandom.randDouble() < delta_ ? 1 : 2;
				Solution child = normalReproduce(taskId, i, type);
				problemSet_.get(taskId).evaluate(child);
				evaluations_++;

				updateReference(child, taskId);
				betterCount += updateProblem(child, i, type, taskId);
			}
		}
		preBetterCount[taskId] = betterCount;
	}

	public void transferConverge(int targetTaskId, int sourceTaskId) throws JMException {
		int[] betterCount = new int[2];
		// Async
		solelyConverge(sourceTaskId, aStep);

		transferTimes[targetTaskId][sourceTaskId] ++;

		for (int i = 0; i < populationSize_; i++) {
			int type = PseudoRandom.randDouble() < delta_ ? 1 : 2;
			Solution child;

			// 0: I; 1: E
			int transferMode;
			if (PseudoRandom.randDouble() < 1) {
				// 隐式
				transferMode = 0;
				child = transferReproduce(targetTaskId, sourceTaskId, i);
			} else {
				// 显式
				transferMode = 1;
				child = transferReproduce(sourceTaskId, sourceTaskId, i);
			}

			problemSet_.get(targetTaskId).evaluate(child);
			evaluations_++;

			updateReference(child, targetTaskId);
			betterCount[transferMode] += updateProblem(child, i, type, targetTaskId);
		}

		if ((betterCount[0] + betterCount[1]) <= populationSize_ * 0.1 * transferP[targetTaskId][sourceTaskId]) {
			transferP[targetTaskId][sourceTaskId] = Math.max(transferP[targetTaskId][sourceTaskId] * 0.9, 0.05);
			scores[targetTaskId][sourceTaskId] = Math.max(scores[targetTaskId][sourceTaskId] - 1, 0);
			bannedTask[targetTaskId][sourceTaskId] = true;
		} else if ((betterCount[0] + betterCount[1]) > populationSize_ * 0.5 * transferP[targetTaskId][sourceTaskId]){
			transferP[targetTaskId][sourceTaskId] = Math.min(transferP[targetTaskId][sourceTaskId] / 0.9, 1);
			preTransferTask[targetTaskId] = sourceTaskId;
			transferBetterTimes[targetTaskId][sourceTaskId] ++;
//			implicitP[targetTaskId][sourceTaskId] = 0.1 * implicitP[targetTaskId][sourceTaskId] + 0.9 * ((double) betterCount[0] / (betterCount[0] + betterCount[1]));
		} else {
			transferBetterTimes[targetTaskId][sourceTaskId] ++;
		}
//		System.out.println(evaluations_ + ": " + assistTask + " -> " + taskId + ": " + implicitP + "\t" + transferP);
	}

	public int getSourceTaskId(int targetTaskId, String type) throws JMException {
		int sourceTaskId = targetTaskId;

		if (type.equalsIgnoreCase("random")) {
			// random
			while (sourceTaskId == targetTaskId)
				sourceTaskId = PseudoRandom.randInt(0, taskNum_ - 1);
		}
		else if (type.equalsIgnoreCase("wd")) {
			// Wasserstein Distance
			double[] distance = new double[taskNum_];
			double[] finalScore = new double[taskNum_];
			for (int k = 0; k < taskNum_; k++) {
				if (k == targetTaskId) continue;
				distance[k] = WassersteinDistance.getWD2(
						population_[targetTaskId].getMat(),
						population_[k].getMat());
				finalScore[k] = 3 * scores[targetTaskId][k] / distance[k];
			}
			sourceTaskId = Utils.rouletteExceptZero(finalScore);
		} else if (type.equalsIgnoreCase("kl")) {
			double[] distance;
			double[] finalScore = new double[taskNum_];
			KLD kld = new KLD(problemSet_, population_);
			distance = kld.getKDL(targetTaskId);
			for (int k = 0; k < taskNum_; k++){
				finalScore[k] = scores[targetTaskId][sourceTaskId] / distance[k];
			}
			sourceTaskId = Utils.rouletteExceptZero(finalScore);
		} else if (type.equalsIgnoreCase("banned")) {
			if (preTransferTask[targetTaskId] == -1) {
				while (sourceTaskId == targetTaskId && !bannedTask[targetTaskId][sourceTaskId])
					sourceTaskId = PseudoRandom.randInt(0, taskNum_ - 1);
			} else {
				sourceTaskId = preTransferTask[targetTaskId];
			}
		}

		return sourceTaskId;
	}

	public void initUniformWeight() { // init lambda vectors
		for (int k = 0; k < taskNum_; k++) {
			int nw = 0;
			if (problemSet_.get(k).getNumberOfObjectives() == 2) {
				for (int n = 0; n < populationSize_; n++) {
					double a = 1.0 * n / (populationSize_ - 1);
					lambda_[k][n][0] = a;
					lambda_[k][n][1] = 1 - a;
					nw++;
				} // for
			} // if
			else {
				String dataFileName;
				dataFileName = "W" + problemSet_.get(k).getNumberOfObjectives() + "D_"
						+ populationSize_ + ".dat";

				try {
					// Open the file
					FileInputStream fis = new FileInputStream(dataDirectory_ + "/"
							+ dataFileName);
					InputStreamReader isr = new InputStreamReader(fis);
					BufferedReader br = new BufferedReader(isr);

					int i = 0;
					int j;
					String aux = br.readLine();
					while (aux != null) {
						StringTokenizer st = new StringTokenizer(aux);
						j = 0;
						while (st.hasMoreTokens()) {
							double value = Double.parseDouble(st.nextToken());
							lambda_[k][i][j] = value;
							j++;
						}
						aux = br.readLine();
						i++;
					}
					br.close();
				} catch (Exception e) {
					System.out
							.println("initUniformWeight: failed when reading for file: "
									+ dataDirectory_ + "/" + dataFileName);
					e.printStackTrace();
				}
			} // else
		}
	} // initUniformWeight

	/**
   * 
   */
	public void initNeighborhood() {
		for (int k = 0; k < taskNum_; k++) {
			double[] x = new double[populationSize_];
			int[] idx = new int[populationSize_];

			for (int i = 0; i < populationSize_; i++) {
				// calculate the distances based on weight vectors
				for (int j = 0; j < populationSize_; j++) {
					x[j] = Utils.distVector(lambda_[k][i], lambda_[k][j]);
					// x[j] = dist_vector(population[i].namda,population[j].namda);
					idx[j] = j;
					// System.out.println("x["+j+"]: "+x[j]+
					// ". idx["+j+"]: "+idx[j]) ;
				} // for

				// find 'niche' nearest neighboring subproblems
				Utils.minFastSort(x, idx, populationSize_, T_);
				// minfastsort(x,idx,population.size(),niche);

				System.arraycopy(idx, 0, neighborhood_[k][i], 0, T_);
			} // for
		}
	} // initNeighborhood

	/**
   * 
   */
	public void initPopulation() throws JMException, ClassNotFoundException {
		for (int k = 0; k < taskNum_; k++) {
			population_[k] = new SolutionSet(populationSize_);
			for (int i = 0; i < populationSize_; i++) {
				Solution newSolution = new Solution(problemSet_);

				problemSet_.get(k).evaluate(newSolution);
				evaluations_++;
				population_[k].add(newSolution);
			} // for
		}
	} // initPopulation

	/**
   * 
   */
	void initIdealPoint() throws JMException, ClassNotFoundException {
		for (int k = 0; k < taskNum_; k++) {
			for (int i = 0; i < problemSet_.get(k).getNumberOfObjectives(); i++) {
				z_[k][i] = 1.0e+30;
				indArray_[k][i] = new Solution(problemSet_);
				problemSet_.get(k).evaluate(indArray_[k][i]);
				evaluations_++;
			} // for

			for (int i = 0; i < populationSize_; i++) {
				updateReference(population_[k].get(i), k);
			} // for
		}
	} // initIdealPoint

	/**
   * 
   */
	public void matingSelection(Vector<Integer> list, int cid, int size, int type, int taskId) {
		int ss;
		int r;
		int p;

		ss = neighborhood_[taskId][cid].length;
		while (list.size() < size) {
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood_[taskId][cid][r];
				// p = population[cid].table[r];
			} else {
				p = PseudoRandom.randInt(0, populationSize_ - 1);
			}
			boolean flag = true;
			for (int i = 0; i < list.size(); i++) {
				if (list.get(i) == p) // p is in the list
				{
					flag = false;
					break;
				}
			}

			// if (flag) list.push_back(p);
			if (flag) {
				list.addElement(p);
			}
		}
	} // matingSelection

	/**
	 * 
	 * @param individual
	 */
	boolean updateReference(Solution individual, int taskId) {
		boolean isUpdate = false;
		int objOffset = problemSet_.get(taskId).getStartObjPos();
		for (int n = 0; n < problemSet_.get(taskId).getNumberOfObjectives(); n++) {
			if (individual.getObjective(n + objOffset) < z_[taskId][n]) {
				isUpdate = true;
				z_[taskId][n] = individual.getObjective(n + objOffset);

				indArray_[taskId][n] = individual;
			}
		}
		return isUpdate;
	} // updateReference

	/**
	 * @param indiv
	 * @param id
	 * @param type
	 */
	int updateProblem(Solution indiv, int id, int type, int taskId) {
		int size;
		int time;

		time = 0;

		if (type == 1) {
			size = neighborhood_[taskId][id].length;
		} else {
			size = population_[taskId].size();
		}
		int[] perm = new int[size];

		Utils.randomPermutation(perm, size);

		for (int i = 0; i < size; i++) {
			int k;
			if (type == 1) {
				k = neighborhood_[taskId][id][perm[i]];
			} else {
				k = perm[i]; // calculate the values of objective function
								// regarding the current subproblem
			}
			double f1, f2;

			f1 = fitnessFunction(population_[taskId].get(k), lambda_[taskId][k], taskId);
			f2 = fitnessFunction(indiv, lambda_[taskId][k], taskId);

			if (f2 < f1) {
				population_[taskId].replace(k, new Solution(indiv));
				time++;
			}
			// the maximal number of solutions updated is not allowed to exceed
			// 'limit'
			if (time >= nr_) {
				break;
			}
		}

		return time;
	} // updateProblem

	double fitnessFunction(Solution individual, double[] lambda, int taskId) {
		double fitness;
		fitness = 0.0;

		int objOffset = problemSet_.get(taskId).getStartObjPos();
		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;

			for (int n = 0; n < problemSet_.get(taskId).getNumberOfObjectives(); n++) {
				double diff = Math.abs(individual.getObjective(n + objOffset) - z_[taskId][n]);

				double feval;
				if (lambda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lambda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for

			fitness = maxFun;
		} // if
		else {
			System.out.println("MOEAD.fitnessFunction: unknown type "
					+ functionType_);
			System.exit(-1);
		}
		return fitness;
	} // fitnessEvaluation
} // MOEAD

