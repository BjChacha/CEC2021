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

package etmo.metaheuristics.matmy3;

import etmo.core.*;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.math.Distance;
import etmo.util.math.Matrix;
import etmo.util.math.Probability;
import etmo.util.math.Random;
import etmo.util.math.Vector;
import etmo.util.sorting.NDSortiong;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.StringTokenizer;

public class MaTMY3_MOEAD extends MtoAlgorithm {
	private int populationSize_;
	private int taskNum;
    private int varNum;
	private SolutionSet[] population_;
	
	double[][] z_;
	double[][][] lambda_;
	int T_;
	int[][][] neighborhood_;
	double delta_;
	int nr_;
	Solution[][] indArray_;
	String functionType_;

	private int[] objStart;
	private int[] objEnd;

	double[][] transferP;
	double[][] distances1;
	double[][] distances2;

	private double mutationProbability;
    private double transferProbability;
	private double elitePart;

	int[][] transferVol;

	int generations_;
	int evaluations_;
	int maxEvaluations_;

	double[][] means;
    double[][] stds;
    double[][][] sigmas;

    private Operator DECrossover;
    private Operator SBXCrossover;
    private Operator mutation;

	String dataDirectory_;

	int taskNum_;

	public MaTMY3_MOEAD(ProblemSet problemSet) {
		super(problemSet);

		functionType_ = "_TCHE1";

	} // MOEAD

	private void initState() {
		taskNum_ = problemSet_.size();
		evaluations_ = 0;
		generations_ = 0;
		taskNum = problemSet_.size();
        varNum = problemSet_.getMaxDimension();

		maxEvaluations_ = (Integer) this.getInputParameter("maxEvaluations");
		populationSize_ = (Integer) this.getInputParameter("populationSize");
		dataDirectory_ = this.getInputParameter("dataDirectory").toString();
		T_ = (Integer) this.getInputParameter("T");
		nr_ = (Integer) this.getInputParameter("nr");
		delta_ = (Double) this.getInputParameter("delta");

        transferProbability = (Double) this.getInputParameter("transferProbability");
        mutationProbability = (Double) this.getInputParameter("mutationProbability");
		elitePart = (Double) this.getInputParameter("elitePartition");

        DECrossover = operators_.get("DECrossover");
        SBXCrossover = operators_.get("SBXCrossover");
        mutation = operators_.get("mutation");

		objStart = new int[taskNum];
        objEnd = new int[taskNum];

		population_ = new SolutionSet[taskNum_];
		indArray_ = new Solution[taskNum_][];
		neighborhood_ = new int[taskNum_][populationSize_][T_];
		z_ = new double[taskNum_][];
		lambda_ = new double[taskNum_][][];

		transferP = new double[taskNum_][taskNum_];
		distances1 = new double[taskNum_][taskNum_];
		distances2 = new double[taskNum_][taskNum_];

		means = new double[taskNum][varNum];
        stds = new double[taskNum][varNum];
        sigmas = new double[taskNum][varNum][varNum];

		for (int k = 0; k < taskNum_; k++){
			indArray_[k] = new Solution[problemSet_.get(k).getNumberOfObjectives()];
			z_[k] = new double[problemSet_.get(k).getNumberOfObjectives()];
			lambda_[k] = new double[populationSize_][problemSet_.get(k).getNumberOfObjectives()];

			objStart[k] = problemSet_.get(k).getStartObjPos();
            objEnd[k] = problemSet_.get(k).getEndObjPos();

            Arrays.fill(distances1[k], 0);
            Arrays.fill(distances2[k], 0);
			Arrays.fill(transferP[k], transferProbability);

			Arrays.fill(means[k], 0);
            Arrays.fill(stds[k], 0);
		}
	}

	private Solution normalReproduce(int taskId, int currentId, int type) throws JMException {
		Solution child;

		java.util.Vector<Integer> p = new java.util.Vector<Integer>();
		matingSelection(p, currentId, 2, type, taskId);

		Solution[] parents = new Solution[3];

		parents[0] = population_[taskId].get(p.get(0));
		parents[1] = population_[taskId].get(p.get(1));
		parents[2] = population_[taskId].get(currentId);

		child = (Solution) DECrossover.execute(new Object[]{
				population_[taskId].get(currentId), parents});
		mutateIndividual(child);
		return child;
	}

	private Solution transferReproduce(int targetTaskId, int sourceTaskId, int currentId) throws JMException {
		Solution child;
		child = new Solution(population_[sourceTaskId].get(currentId));

		double[] tmpMean = population_[sourceTaskId].getMean();
		Vector.vecSub_(tmpMean, means[sourceTaskId]);
		Vector.vecSub_(tmpMean, means[targetTaskId]);
		double[] tmpStd = population_[sourceTaskId].getStd();
		double[] newFeatures = Probability.sampleByNorm(tmpMean, tmpStd);
		Vector.vecClip_(newFeatures, 0.0, 1.0);
		
		child.setDecisionVariables(newFeatures);
		mutateIndividual(child);
		child.resetObjective();

		return child;
	}

	private void mutateIndividual(Solution child) throws JMException {
		if (PseudoRandom.randDouble() < mutationProbability)
			mutation.execute(child);
		return;
	}

	public SolutionSet[] execute() throws JMException, ClassNotFoundException {
		initState();
		initUniformWeight();
		initNeighborhood();
		initPopulation();
		initIdealPoint();

		generations_ = 1;
		while (evaluations_ < maxEvaluations_) {
			iterate();
			generations_ ++;
		}
		return population_;
	}

	public void iterate() throws JMException {
		updateDistributions(elitePart);
		updateDistances();

		for (int targetTaskId = 0; targetTaskId < taskNum_; targetTaskId++) {
			int pIdx = 0;
			int[] permutation = new int[populationSize_];
			Utils.randomPermutation(permutation, populationSize_);

			// transfer
			for (int taskId = 0; taskId < taskNum_; taskId++){
				for (int i = 0; i < population_[taskId].size(); i++){
					Solution child;
					int type;
					if (PseudoRandom.randDouble() < 0.5) {
						type = 2;
						int sourceTaskId = getSourceTaskId(taskId);
						child = transferReproduce(targetTaskId, sourceTaskId, permutation[pIdx]);
					} else {
						type = PseudoRandom.randDouble() < delta_ ? 1 : 2;
						child = normalReproduce(targetTaskId, permutation[pIdx], type);
					}
					child.setFlag(taskId);
					problemSet_.get(targetTaskId).evaluate(child);
					evaluations_ ++;
					updateReference(child, targetTaskId);
					updateProblem(child, permutation[pIdx], type, targetTaskId);
				}
			}
		}
	}

	void updateDistances() throws JMException {
        for (int k = 0; k < taskNum; k++) {
           final int srcTaskID = k;
           Arrays.parallelSetAll(distances1[k], trgTaskID -> {
                double dist = 0;
                if (trgTaskID > srcTaskID) {
                        dist = Distance.getCorrelationMatrixDistance(sigmas[srcTaskID], sigmas[trgTaskID]);    
                    } else {
                        dist = distances1[trgTaskID][srcTaskID];
                    }
                return dist;
            });
            Arrays.parallelSetAll(distances2[k], trgTaskID -> {
            double dist = 0;
            if (trgTaskID > srcTaskID) {
                    dist = Distance.getDistance(means[srcTaskID], means[trgTaskID]);    
                } else {
                    dist = distances2[trgTaskID][srcTaskID];
                }
            return dist;
            });
        }
    }

	void updateDistributions(double partition) throws JMException {
        SolutionSet[] tmpSet = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
			SolutionSet copy = population_[k].copy();
			NDSortiong.sort(copy, problemSet_, k);
			int size = (int)(population_[k].size() * partition);
			tmpSet[k] = new SolutionSet(size);
			for (int i = 0; i < size; i++) {
				tmpSet[k].add(copy.get(i));
			}
		}
        Arrays.parallelSetAll(means, k -> tmpSet[k].getMean());
        // Arrays.parallelSetAll(stds, k -> tmpSet[k].getStd());
        Arrays.parallelSetAll(sigmas, k -> {
            double[][] output = null;
            try {
				output = Matrix.getMatSigma(tmpSet[k].getMat());
			} catch (JMException e) {
				e.printStackTrace();
			}
            return output;
		});
    }

	public int getSourceTaskId(int targetTaskId) throws JMException {
		int sourceTaskId = targetTaskId;

        // CMD + EMD
        double[] scores = new double[taskNum];
        // CMD
        Arrays.setAll(scores, i -> distances1[targetTaskId][i]);
        int res1 = Random.rouletteWheel(scores, targetTaskId);
        // EMD
        Arrays.setAll(scores, i -> 1 / distances2[targetTaskId][i]);
        int res2 = Random.rouletteWheel(scores, targetTaskId);
        sourceTaskId = PseudoRandom.randDouble() < 0.5 ? res1 : res2;

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

	public void matingSelection(java.util.Vector<Integer> list, int cid, int size, int type, int taskId) {
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

