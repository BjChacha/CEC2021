//  SolutionSet.Java
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

package etmo.core;

import java.io.*;
import java.util.*;

import etmo.util.Configuration;
import etmo.util.JMException;

/**
 * Class representing a SolutionSet (a set of solutions)
 */
public class SolutionSet implements Serializable {

	/**
	 * Stores a list of <code>solution</code> objects.
	 */
	protected final List<Solution> solutionsList_;

	/**
	 * Maximum size of the solution set
	 */
	private int capacity_ = Integer.MAX_VALUE;

	/**
	 * Constructor. Creates an unbounded solution set.
	 */
	public SolutionSet() {
		solutionsList_ = new ArrayList<Solution>();
	} // SolutionSet

	/**
	 * Creates a empty solutionSet with a maximum capacity.
	 * 
	 * @param maximumSize
	 *            Maximum size.
	 */
	public SolutionSet(int maximumSize) {
		solutionsList_ = new ArrayList<Solution>();
		capacity_ = maximumSize;
	} // SolutionSet

	public SolutionSet(List<Solution> solutionList){
		solutionsList_ = solutionList;
		capacity_ = solutionList.size();
	}

	/**
	 * Inserts a new solution into the SolutionSet.
	 * 
	 * @param solution
	 *            The <code>Solution</code> to store
	 * @return True If the <code>Solution</code> has been inserted, false
	 *         otherwise.
	 */
	public boolean add(Solution solution) {
		if (solutionsList_.size() == capacity_) {
			Configuration.logger_.severe("The population is full");
			Configuration.logger_.severe("Capacity is : " + capacity_);
			Configuration.logger_.severe("\t Size is: " + this.size());
			return false;
		} // if

		solutionsList_.add(solution);
		return true;
	} // add

	public boolean add(int index, Solution solution) {
		solutionsList_.add(index, solution);
		return true;
	}

	/*
	 * public void add(Solution solution) { if (solutionsList_.size() ==
	 * capacity_) try { throw new JMException(
	 * "SolutionSet.Add(): the population is full") ; } catch (JMException e) {
	 * e.printStackTrace(); } else solutionsList_.add(solution); }
	 */
	/**
	 * Returns the ith solution in the set.
	 * 
	 * @param i
	 *            Position of the solution to obtain.
	 * @return The <code>Solution</code> at the position i.
	 * @throws IndexOutOfBoundsException
	 *             Exception
	 */
	public Solution get(int i) {
		if (i >= solutionsList_.size()) {
			throw new IndexOutOfBoundsException("Index out of Bound " + i);
		}
		return solutionsList_.get(i);
	} // get

	/**
	 * Returns the maximum capacity of the solution set
	 * 
	 * @return The maximum capacity of the solution set
	 */
	public int getMaxSize() {
		return capacity_;
	} // getMaxSize

	/**
	 * Sorts a SolutionSet using a <code>Comparator</code>.
	 * 
	 * @param comparator
	 *            <code>Comparator</code> used to sort.
	 */
	public void sort(Comparator comparator) {
		if (comparator == null) {
			Configuration.logger_.severe("No criterium for comparing exist");
			return;
		} // if
		Collections.sort(solutionsList_, comparator);
	} // sort

	/**
	 * Returns the index of the best Solution using a <code>Comparator</code>.
	 * If there are more than one occurrences, only the index of the first one
	 * is returned
	 * 
	 * @param comparator
	 *            <code>Comparator</code> used to compare solutions.
	 * @return The index of the best Solution attending to the comparator or
	 *         <code>-1<code> if the SolutionSet is empty
	 */
	int indexBest(Comparator comparator) {
		if ((solutionsList_ == null) || (this.solutionsList_.isEmpty())) {
			return -1;
		}

		int index = 0;
		Solution bestKnown = solutionsList_.get(0), candidateSolution;
		int flag;
		for (int i = 1; i < solutionsList_.size(); i++) {
			candidateSolution = solutionsList_.get(i);
			flag = comparator.compare(bestKnown, candidateSolution);
			if (flag == +1) {
				index = i;
				bestKnown = candidateSolution;
			}
		}

		return index;
	} // indexBest

	/**
	 * Returns the best Solution using a <code>Comparator</code>. If there are
	 * more than one occurrences, only the first one is returned
	 * 
	 * @param comparator
	 *            <code>Comparator</code> used to compare solutions.
	 * @return The best Solution attending to the comparator or <code>null
	 *         <code> if the SolutionSet is empty
	 */
	public Solution best(Comparator comparator) {
		int indexBest = indexBest(comparator);
		if (indexBest < 0) {
			return null;
		} else {
			return solutionsList_.get(indexBest);
		}

	} // best

	/**
	 * Returns the index of the worst Solution using a <code>Comparator</code>.
	 * If there are more than one occurrences, only the index of the first one
	 * is returned
	 * 
	 * @param comparator
	 *            <code>Comparator</code> used to compare solutions.
	 * @return The index of the worst Solution attending to the comparator or
	 *         <code>-1<code> if the SolutionSet is empty
	 */
	public int indexWorst(Comparator comparator) {
		if ((solutionsList_ == null) || (this.solutionsList_.isEmpty())) {
			return -1;
		}

		int index = 0;
		Solution worstKnown = solutionsList_.get(0), candidateSolution;
		int flag;
		for (int i = 1; i < solutionsList_.size(); i++) {
			candidateSolution = solutionsList_.get(i);
			flag = comparator.compare(worstKnown, candidateSolution);
			if (flag == -1) {
				index = i;
				worstKnown = candidateSolution;
			}
		}

		return index;

	} // indexWorst

	/**
	 * Returns the worst Solution using a <code>Comparator</code>. If there are
	 * more than one occurrences, only the first one is returned
	 * 
	 * @param comparator
	 *            <code>Comparator</code> used to compare solutions.
	 * @return The worst Solution attending to the comparator or <code>null
	 *         <code> if the SolutionSet is empty
	 */
	public Solution worst(Comparator comparator) {

		int index = indexWorst(comparator);
		if (index < 0) {
			return null;
		} else {
			return solutionsList_.get(index);
		}

	} // worst

	/**
	 * Returns the number of solutions in the SolutionSet.
	 * 
	 * @return The size of the SolutionSet.
	 */
	public int size() {
		return solutionsList_.size();
	} // size

	/**
	 * Writes the objective function values of the <code>Solution</code> objects
	 * into the set in a file.
	 * 
	 * @param path
	 *            The output file name
	 */
	public void printObjectivesToFile(String path) {
		try {
			/* Open the file */
			path = "data//" + path;
			FileOutputStream fos = new FileOutputStream(path);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			BufferedWriter bw = new BufferedWriter(osw);

			for (Solution aSolutionsList_ : solutionsList_) {
				// if (this.vector[i].getFitness()<1.0) {
				bw.write(aSolutionsList_.toString());
				bw.newLine();
				// }
			}

			/* Close the file */
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	} // printObjectivesToFile

	/**
	 * Writes the decision encodings.variable values of the
	 * <code>Solution</code> solutions objects into the set in a file.
	 * 
	 * @param path
	 *            The output file name
	 */
	public void printVariablesToFile(String path) {
		try {
			FileOutputStream fos = new FileOutputStream(path);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			BufferedWriter bw = new BufferedWriter(osw);

			if (size() > 0) {
				int numberOfVariables = solutionsList_.get(0).getDecisionVariables().length;
				for (Solution aSolutionsList_ : solutionsList_) {
					for (int j = 0; j < numberOfVariables; j++)
						bw.write(aSolutionsList_.getDecisionVariables()[j].toString() + " ");
					bw.newLine();
				}
			}
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	} // printVariablesToFile

	/**
	 * Write the function values of feasible solutions into a file
	 * 
	 * @param path
	 *            File name
	 */
	public void printFeasibleFUN(String path) {
		try {
			FileOutputStream fos = new FileOutputStream(path);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			BufferedWriter bw = new BufferedWriter(osw);

			for (Solution aSolutionsList_ : solutionsList_) {
				if (aSolutionsList_.getOverallConstraintViolation() == 0.0) {
					bw.write(aSolutionsList_.toString());
					bw.newLine();
				}
			}
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	}

	/**
	 * Write the encodings.variable values of feasible solutions into a file
	 * 
	 * @param path
	 *            File name
	 */
	public void printFeasibleVAR(String path) {
		try {
			FileOutputStream fos = new FileOutputStream(path);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			BufferedWriter bw = new BufferedWriter(osw);

			if (size() > 0) {
				int numberOfVariables = solutionsList_.get(0).getDecisionVariables().length;
				for (Solution aSolutionsList_ : solutionsList_) {
					if (aSolutionsList_.getOverallConstraintViolation() == 0.0) {
						for (int j = 0; j < numberOfVariables; j++)
							bw.write(aSolutionsList_.getDecisionVariables()[j].toString() + " ");
						bw.newLine();
					}
				}
			}
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	}

	/**
	 * Empties the SolutionSet
	 */
	public void clear() {
		solutionsList_.clear();
	} // clear

	/**
	 * Deletes the <code>Solution</code> at position i in the set.
	 * 
	 * @param i
	 *            The position of the solution to remove.
	 */
	public void remove(int i) {
		if (i > solutionsList_.size() - 1) {
			Configuration.logger_.severe("Size is: " + this.size());
		} // if
		solutionsList_.remove(i);
	} // remove

	/**
	 * Returns an <code>Iterator</code> to access to the solution set list.
	 * 
	 * @return the <code>Iterator</code>.
	 */
	public Iterator<Solution> iterator() {
		return solutionsList_.iterator();
	} // iterator

	/**
	 * Returns a new <code>SolutionSet</code> which is the result of the union
	 * between the current solution set and the one passed as a parameter.
	 * 
	 * @param solutionSet
	 *            SolutionSet to join with the current solutionSet.
	 * @return The result of the union operation.
	 */
	public SolutionSet union(SolutionSet solutionSet) {
		// Check the correct size. In development
		int newSize = this.size() + solutionSet.size();
		if (newSize < capacity_)
			newSize = capacity_;

		// Create a new population
		SolutionSet union = new SolutionSet(newSize);
		for (int i = 0; i < this.size(); i++) {
			union.add(this.get(i));
		} // for

		for (int i = this.size(); i < (this.size() + solutionSet.size()); i++) {
			union.add(solutionSet.get(i - this.size()));
		} // for

		return union;
	} // union

	/**
	 * Replaces a solution by a new one
	 * 
	 * @param position
	 *            The position of the solution to replace
	 * @param solution
	 *            The new solution
	 */
	public void replace(int position, Solution solution) {
		if (position > this.solutionsList_.size()) {
			solutionsList_.add(solution);
		} // if
		solutionsList_.remove(position);
		solutionsList_.add(position, solution);
	} // replace

	/**
	 * Copies the objectives of the solution set to a matrix
	 * 
	 * @return A matrix containing the objectives
	 */
	public double[][] writeObjectivesToMatrix() {
		if (this.size() == 0) {
			return null;
		}
		double[][] objectives;
		objectives = new double[size()][get(0).getNumberOfObjectives()];
		for (int i = 0; i < size(); i++) {
			for (int j = 0; j < get(0).getNumberOfObjectives(); j++) {
				objectives[i][j] = get(i).getObjective(j);
			}
		}
		return objectives;
	} // writeObjectivesMatrix

	public void printObjectives() {
		for (int i = 0; i < solutionsList_.size(); i++)
			System.out.println("" + solutionsList_.get(i));
	}

	public void setCapacity(int capacity) {
		capacity_ = capacity;
	}

	public int getCapacity() {
		return capacity_;
	}

//	add from moeac
	public Solution getCentroidVector() {
		int objectNumb = 0;
		if (solutionsList_.size() != 0) {
			objectNumb = solutionsList_.get(0).getNumberOfObjectives();
		} else {
			System.out.println("solutionsList size == 0");
			System.exit(0);
		}
		int objStart = 0;
		if (solutionsList_.get(0).getSkillFactor() >= 0) {
			objStart = solutionsList_.get(0).getProblemSet().get(solutionsList_.get(0).getSkillFactor()).getStartObjPos();
			objectNumb = solutionsList_.get(0).getProblemSet().get(solutionsList_.get(0).getSkillFactor()).getNumberOfObjectives();
		}

		Solution sol = new Solution(objectNumb);
		double sumValue = 0.0;
		for (int m = 0; m < objectNumb; m++) {
			double value = 0;
			for (int k = 0; k < solutionsList_.size(); k++) {
				value += solutionsList_.get(k).getNormalizedObjective(m + objStart);
			}
			value = value / solutionsList_.size();
			sol.setNormalizedObjective(m, value);
			sumValue += value;
		}

		double normDistance = 0.0;
		for (int i = 0; i < objectNumb; i++) {
			double unitValue = 0.0;
			normDistance += sol.getNormalizedObjective(i) * sol.getNormalizedObjective(i);
			unitValue = sol.getNormalizedObjective(i) / sumValue;
			sol.setUnitHyperplaneObjective(i, unitValue);
		}
		normDistance = Math.sqrt(normDistance);

		sol.setDistanceToIdealPoint(normDistance);
		return sol;
	}

    public double getMeanOfIdx(int idx) throws JMException {
		double sum = 0;
		for (int i = 0; i < solutionsList_.size(); i++){
			sum += solutionsList_.get(i).getDecisionVariables(idx);
		}
		return sum / solutionsList_.size();
	}

	public double getStdOfIdx(int idx) throws JMException {
		double sum = 0;
		double mean = getMeanOfIdx(idx);
		for (int i = 0; i < solutionsList_.size(); i++){
			sum += Math.pow(solutionsList_.get(i).getDecisionVariables(idx) - mean, 2);
		}
		return Math.sqrt(sum / solutionsList_.size());
	}

	public double[] getMean(){
		double[] mean = new double[0];
		try {
			mean = new double[solutionsList_.get(0).numberOfVariables()];
			for (int i = 0; i < mean.length; i++){
				mean[i] = getMeanOfIdx(i);
			}
			return mean;
		} catch (JMException e) {
			e.printStackTrace();
		}
		return mean;
	}

	public double[] getStd(){
		double[] std = new double[0];
		try {
			std = new double[solutionsList_.get(0).numberOfVariables()];
			for (int i = 0; i < std.length; i++){
				std[i] = getStdOfIdx(i);
			}
			return std;
		} catch (JMException e) {
			e.printStackTrace();
		}
		return std;
	}

	public double[][] getMat() throws JMException {
		int numVar = solutionsList_.get(0).numberOfVariables();
		double[][] mat = new double[solutionsList_.size()][numVar];
		for (int i = 0; i < solutionsList_.size(); i++){
			for (int j = 0; j < numVar; j++){
				mat[i][j] = solutionsList_.get(i).getDecisionVariables(j);
			}
		}
		return mat;
	}

	public double[] getBestObjectiveVector(){
		return getBestObjectiveVector(0);
	}

	public double[] getBestObjectiveVector(int task){
		Problem problem = solutionsList_.get(0).getProblemSet().get(task);
		double[] vector = new double[problem.getNumberOfObjectives()];
		for (int j = 0; j <= problem.getEndObjPos() - problem.getStartObjPos(); j++) {
			vector[j] = Double.MAX_VALUE;
			for (int i = 0; i < solutionsList_.size(); i++)
				vector[j] = Math.min(vector[j], solutionsList_.get(i).getObjective(problem.getStartObjPos() + j));
		}
		return vector;
	}

	public double[] getAverageObjectiveVector(){
		return getAverageObjectiveVector(0);
	}

	public double[] getAverageObjectiveVector(int task){
		Problem problem = solutionsList_.get(0).getProblemSet().get(task);
		double[] vector = new double[problem.getNumberOfObjectives()];
		for (int j = 0; j <= problem.getEndObjPos() - problem.getStartObjPos(); j++) {
			double sum = 0;
			for (int i = 0; i < solutionsList_.size(); i++)
				sum += solutionsList_.get(i).getObjective(problem.getStartObjPos() + j);
			vector[j] = sum / solutionsList_.size();
		}
		return vector;
	}

	public double[] getAverageDecisionVector() throws JMException {
		double[] res = new double[solutionsList_.get(0).numberOfVariables()];
		for (int i = 0; i < res.length; i++){
			double sum = 0;
			for (int j = 0; j < solutionsList_.size(); j++){
				sum += solutionsList_.get(j).getDecisionVariables(i);
			}
			res[i] = sum / solutionsList_.size();
		}
		return res;
	}

	public double[] getNDWeightedDecisionVector() throws JMException {
		double[] res = new double[solutionsList_.get(0).numberOfVariables()];
		double weightSum = 0;
		for (int i = 0; i < res.length; i++){
			double sum = 0;
			for (int j = 0; j < solutionsList_.size(); j++){
				sum += solutionsList_.get(j).getDecisionVariables(i)
						* 1 / (solutionsList_.get(j).getRank() + 1);
				weightSum += 1 / (solutionsList_.get(j).getRank() + 1);
			}
			res[i] = sum / solutionsList_.size() / weightSum;
		}
		return res;
	}

	public SolutionSet copy(){
		SolutionSet res = new SolutionSet(solutionsList_.size());
		for (int i = 0; i < solutionsList_.size(); i++){
			Solution individual = new Solution(solutionsList_.get(i));
			res.add(individual);
		}
		return res;
	}

	public boolean isFull(){
		return solutionsList_.size() == capacity_;
	}

	private boolean remove;
	public boolean isRemove() {
		return remove;
	}
	public void setRemove(boolean remove) {
		this.remove = remove;
	}
} // SolutionSet
