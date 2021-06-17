//  NSGAII.java
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

package etmo.metaheuristics.matbml.libs;

import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.Ranking;
import etmo.util.comparators.CrowdingComparator;

/**
 * Implementation of NSGA-II. This implementation of NSGA-II makes use of a
 * QualityIndicator object to obtained the convergence speed of the algorithm.
 * This version is used in the paper: A.J. Nebro, J.J. Durillo, C.A. Coello
 * Coello, F. Luna, E. Alba
 * "A Study of Convergence Speed in Multi-Objective Metaheuristics." To be
 * presented in: PPSN'08. Dortmund. September 2008.
 */

public class NSGAII extends MaTAlgorithm {
	int populationSize;
	int maxEvaluations;
	int evaluations;

	SolutionSet population;
	SolutionSet offspringPopulation;
	SolutionSet union;

	Operator mutationOperator;
	Operator crossoverOperator;
	Operator selectionOperator;

	Distance distance;

	boolean[] selected;

	/**
	 * Constructor
	 * 
	 * @param problemSet
	 *            Problem to solve
	 */
	public NSGAII(ProblemSet problemSet) {
		super(problemSet);
	//	System.out.println("sup: " + problemSet_.get(0).getHType());
	} // NSGAII

	public NSGAII(ProblemSet problemSet, SolutionSet solutionSet){
		super(problemSet);
		population = solutionSet;
		populationSize = solutionSet.size();
	}

	public NSGAII(ProblemSet problemSet, SolutionSet solutionSet, int taskIdx){
		super(problemSet);
		population = solutionSet;
		populationSize = solutionSet.size();
		taskIdx_ = taskIdx;
		selected = new boolean[problemSet.getTotalNumberOfObjs()];
		for (int i = 0; i < selected.length; i++){
			if (i < problemSet.get(taskIdx_).getStartObjPos() || i > problemSet.get(taskIdx_).getEndObjPos())
				selected[i] = false;
			else
				selected[i] = true;
		}
	}

	/**
	 * Runs the NSGA-II algorithm.
	 * 
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 *         solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initState();
		initPopulation();
		while (evaluations < maxEvaluations) {
			// Create the offSpring solutionSet
			offspringPopulation = new SolutionSet(populationSize);
			Solution[] parents = new Solution[2];
			for (int i = 0; i < (populationSize / 2); i++) {
				if (evaluations < maxEvaluations) {
					// obtain parents
					parents[0] = (Solution) selectionOperator.execute(population);
					parents[1] = (Solution) selectionOperator.execute(population);
					Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
					mutationOperator.execute(offSpring[0]);
					mutationOperator.execute(offSpring[1]);
					problemSet_.get(taskIdx_).evaluate(offSpring[0]);
					problemSet_.get(taskIdx_).evaluateConstraints(offSpring[0]);
					problemSet_.get(taskIdx_).evaluate(offSpring[1]);
					problemSet_.get(taskIdx_).evaluateConstraints(offSpring[1]);
					offspringPopulation.add(offSpring[0]);
					offspringPopulation.add(offSpring[1]);
					evaluations += 2;

//					if (evaluations % (populationSize * 20) == 0){
//						LogPopulation.LogPopulation("NSGAII", population, problemSet_, evaluations, false);
//					}
				} // if
			} // for

			// Create the solutionSet union of solutionSet and offSpring
			union = ((SolutionSet) population).union(offspringPopulation);

			// Ranking the union
			Ranking ranking = new Ranking(union);

			int remain = populationSize;
			int index = 0;
			SolutionSet front = null;
			population.clear();

			// Obtain the next front
			front = ranking.getSubfront(index);

			while ((remain > 0) && (remain >= front.size())) {
				// Assign crowding distance to individuals
				distance.crowdingDistanceAssignment(front, problemSet_.get(taskIdx_).getNumberOfObjectives());
				// Add the individuals of this front
				for (int k = 0; k < front.size(); k++) {
					population.add(front.get(k));
				} // for

				// Decrement remain
				remain = remain - front.size();

				// Obtain the next front
				index++;
				if (remain > 0) {
					front = ranking.getSubfront(index);
				} // if
			} // while

			// Remain is less than front(index).size, insert only the best one
			if (remain > 0) { // front contains individuals to insert
				distance.crowdingDistanceAssignment(front, problemSet_.get(taskIdx_).getNumberOfObjectives());
				front.sort(new CrowdingComparator());
				for (int k = 0; k < remain; k++) {
					population.add(front.get(k));
				} // for

				remain = 0;
			} // if

		} // while

		Ranking ranking = new Ranking(population);

		return ranking.getSubfront(0);
	} // execute

	void initPopulation() throws JMException, ClassNotFoundException {
		population = new SolutionSet(populationSize);
		Solution newSolution;
		for (int i = 0; i < populationSize; i++) {
			newSolution = new Solution(problemSet_);
			problemSet_.get(taskIdx_).evaluate(newSolution);
			problemSet_.get(taskIdx_).evaluateConstraints(newSolution);
			evaluations++;
			population.add(newSolution);
		}
	}

	public void initState(){
		distance = new Distance();

		populationSize = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();

		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		selectionOperator = operators_.get("selection");
		// 种群初始化在外部实现，故起始评价次数为种群初始化以后（即种群大小）
		evaluations = populationSize;
	}

	public boolean step() throws JMException {
		offspringPopulation = new SolutionSet(populationSize);
		Solution[] parents = new Solution[2];
		for (int i = 0; i < (populationSize / 2); i++) {
			if (evaluations < maxEvaluations) {
				// obtain parents
				parents[0] = (Solution) selectionOperator.execute(population);
				parents[1] = (Solution) selectionOperator.execute(population);
				Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
				mutationOperator.execute(offSpring[0]);
				mutationOperator.execute(offSpring[1]);
				problemSet_.get(taskIdx_).evaluate(offSpring[0]);
				problemSet_.get(taskIdx_).evaluateConstraints(offSpring[0]);
				problemSet_.get(taskIdx_).evaluate(offSpring[1]);
				problemSet_.get(taskIdx_).evaluateConstraints(offSpring[1]);
				offspringPopulation.add(offSpring[0]);
				offspringPopulation.add(offSpring[1]);
				evaluations += 2;
			} // if
		} // for

		union = ((SolutionSet) population).union(offspringPopulation);

		Ranking ranking = new Ranking(union);

		int remain = populationSize;
		int index = 0;
		SolutionSet front = null;
		population.clear();

		// Obtain the next front
		front = ranking.getSubfront(index);

		while ((remain > 0) && (remain >= front.size())) {
			// Assign crowding distance to individuals
			distance.crowdingDistanceAssignment(front, problemSet_.get(taskIdx_).getNumberOfObjectives(), selected);
			// Add the individuals of this front
			for (int k = 0; k < front.size(); k++) {
				population.add(front.get(k));
			} // for

			// Decrement remain
			remain = remain - front.size();

			// Obtain the next front
			index++;
			if (remain > 0) {
				front = ranking.getSubfront(index);
			} // if
		} // while

		// Remain is less than front(index).size, insert only the best one
		if (remain > 0) { // front contains individuals to insert
			distance.crowdingDistanceAssignment(front, problemSet_.get(taskIdx_).getNumberOfObjectives(), selected);
			front.sort(new CrowdingComparator());
			for (int k = 0; k < remain; k++) {
				population.add(front.get(k));
			} // for

			remain = 0;
		} // if

//		System.out.println(evaluations +"/" + maxEvaluations);
		return evaluations < maxEvaluations;
	}
} // NSGA-II
