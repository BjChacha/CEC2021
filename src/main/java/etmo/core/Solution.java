//  Solution.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Description: 
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
import etmo.util.JMException;
import java.io.Serializable;
import java.util.Arrays;

/**
 * Class representing a solution for a problem.
 */
public class Solution implements Serializable {
	/**
	 * Stores the problem
	 */
	private ProblemSet problemSet_;


	/**
	 * Stores the type of the encodings.variable
	 */
	private SolutionType type_;

	/**
	 * Stores the decision variables of the solution.
	 */
	private Variable[] variable_;

	/**
	 * Stores the objectives values of the solution.
	 */
	private double[] objective_;

	/**
	 * Stores the number of objective values of the solution
	 */
	private int numberOfObjectives_;

	/**
	 * Stores the so called fitness value. Used in some metaheuristics
	 */
	private double fitness_;

	/**
	 * Stores the skill factor
	 */
	private int skillFactor_;

	/**
	 * Stores the so called rank of the solution. Used in NSGA-II
	 */
	private int rank_;

	/**
	 * Stores the overall constraint violation of the solution.
	 */
	private double overallConstraintViolation_;

	/**
	 * Stores the number of constraints violated by the solution.
	 */
	private int numberOfViolatedConstraints_;

	/**
	 * This field is intended to be used to know the location of a solution into
	 * a <code>SolutionSet</code>. Used in MOCell
	 */
	private int location_;

	private int flag;

	private int flag2;

	/**
	 * Stores the distance to his k-nearest neighbor into a
	 * <code>SolutionSet</code>. Used in SPEA2.
	 */
	private double kDistance_;

	/**
	 * Stores the crowding distance of the the solution in a
	 * <code>SolutionSet</code>. Used in NSGA-II.
	 */
	private double crowdingDistance_;

	/**
	 * Stores the distance between this solution and a <code>SolutionSet</code>.
	 * Used in AbySS.
	 */
	private double distanceToSolutionSet_;
	
	private double gFunValue_;


	// Label, for classification training.
	private int label_;

	int LABEL_DEFAULT_VALUE = -1;

//	add from moeac
	private double[] normalizedObjective_;
	private double distanceToIdealPoint;
	private double[] unitHyperplaneObjective_;

	private double pdf_;

	/**
	 * Constructor.
	 */
	public Solution() {
		problemSet_ = null;
		skillFactor_ = -1;
		overallConstraintViolation_ = 0.0;
		numberOfViolatedConstraints_ = 0;
		type_ = null;
		variable_ = null;
		objective_ = null;

		label_ = LABEL_DEFAULT_VALUE;
	} // Solution

	/**
	 * Constructor
	 * 
	 * @param numberOfObjectives
	 *            Number of objectives of the solution
	 * 
	 *            This constructor is used mainly to read objective values from
	 *            a file to variables of a SolutionSet to apply quality
	 *            indicators
	 */
	public Solution(int numberOfObjectives) {
		numberOfObjectives_ = numberOfObjectives;
		objective_ = new double[numberOfObjectives];


		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];
		this.translatedObjectives = new double[this.numberOfObjectives_];
//		lamda_ = this;
		rank1 = -1;
		isMarked = false;

		label_ = LABEL_DEFAULT_VALUE;
	}

	public Solution(ProblemSet problemSet) throws ClassNotFoundException {
		problemSet_ = problemSet;
		type_ = problemSet.getSolutionType();
		numberOfObjectives_ = problemSet.getTotalNumberOfObjs();
		objective_ = new double[numberOfObjectives_];
		
		for (int i = 0; i < numberOfObjectives_; i++)
			objective_[i] = Double.POSITIVE_INFINITY;

		// Setting initial values
		fitness_ = 0.0;
		kDistance_ = 0.0;
		crowdingDistance_ = 0.0;
		skillFactor_ = -1;

		variable_ = type_.createVariables();

//		add from moeac
		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];
		remove = false;
		localDensity_ = -1;
		gaussianLocalDensity_ = 0.0;
		centerDistance_ = 0.0;
		convergenceDistance_ = 0.0;
		distanceToSolutionSet_ = Double.POSITIVE_INFINITY;
		distanceToIdealPoint = Double.POSITIVE_INFINITY;
		distanceToNadirPoint = 0.0;
		sumValue = Double.POSITIVE_INFINITY;

		this.translatedObjectives = new double[this.numberOfObjectives_];
//		lamda_ = this;
		rank1 = -1;
		isMarked = false;

		label_ = LABEL_DEFAULT_VALUE;
	} // Solution

	static public Solution getNewSolution(ProblemSet problemSet) throws ClassNotFoundException {
		return new Solution(problemSet);
	}

	/**
	 * Constructor
	 * 
	 * @param problemSet
	 *            The problem to solve
	 */
	public Solution(ProblemSet problemSet, Variable[] variables) {
		problemSet_ = problemSet;
		type_ = problemSet.getSolutionType();
		numberOfObjectives_ = problemSet.getTotalNumberOfObjs();
		objective_ = new double[numberOfObjectives_];
		
		for (int i = 0; i < numberOfObjectives_; i++)
			objective_[i] = Double.POSITIVE_INFINITY;

		
		// Setting initial values
		fitness_ = 0.0;
		kDistance_ = 0.0;
		crowdingDistance_ = 0.0;

		skillFactor_ = -1;

		distanceToSolutionSet_ = Double.POSITIVE_INFINITY;
		// <-

		variable_ = variables;

		label_ = LABEL_DEFAULT_VALUE;
	} // Constructor

	/**
	 * Copy constructor.
	 * 
	 * @param solution
	 *            Solution to copy.
	 */
	public Solution(Solution solution) {
		problemSet_ = solution.problemSet_;
		type_ = solution.type_;

		numberOfObjectives_ = solution.getNumberOfObjectives();
		objective_ = new double[numberOfObjectives_];
		for (int i = 0; i < objective_.length; i++) {
			objective_[i] = solution.getObjective(i);
		} // for
			// <-
		
		variable_ = type_.copyVariables(solution.variable_);
		overallConstraintViolation_ = solution.getOverallConstraintViolation();
		numberOfViolatedConstraints_ = solution.getNumberOfViolatedConstraint();
		distanceToSolutionSet_ = solution.getDistanceToSolutionSet();
		crowdingDistance_ = solution.getCrowdingDistance();
		kDistance_ = solution.getKDistance();
		fitness_ = solution.getFitness();
		rank_ = solution.getRank();
		location_ = solution.getLocation();

		skillFactor_ = solution.getSkillFactor();


		type_ = solution.type_;

		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];

		for (int i = 0; i < objective_.length; i++) {
			objective_[i] = solution.getObjective(i);
			normalizedObjective_[i] = solution.getNormalizedObjective(i);
			unitHyperplaneObjective_[i] = solution.getUnitHyperplaneObjective(i);
		} // for
		// <-
		this.translatedObjectives = new double[this.numberOfObjectives_];
		remove = false;
		variable_ = type_.copyVariables(solution.variable_);
		overallConstraintViolation_ = solution.getOverallConstraintViolation();
		numberOfViolatedConstraints_ = solution.getNumberOfViolatedConstraint();
		distanceToSolutionSet_ = solution.getDistanceToSolutionSet();
		distanceToIdealPoint = Double.POSITIVE_INFINITY;
		distanceToNadirPoint = 0.0;
		sumValue = Double.POSITIVE_INFINITY;
		crowdingDistance_ = solution.getCrowdingDistance();
		localDensity_ = solution.getLocalDensity_();
		gaussianLocalDensity_ = solution.getGaussianLocalDensity_();
		centerDistance_ = solution.getCenterDistance_();
		convergenceDistance_ = solution.getConvergenceDistance_();
		kDistance_ = solution.getKDistance();
		fitness_ = solution.getFitness();
		marked_ = solution.isMarked_();
		rank_ = solution.getRank();
		location_ = solution.getLocation();
//		max_distance_=  solution.getmaxDistance();
//		min_distance_=  solution.getminDistance();
//		cross_type=  solution.getcross_type();
//		clone_num=  solution.getclone_num();


		rank1 = solution.getRank1();
		isMarked = solution.getRank2();

		label_ = solution.getLabel();
	} // Solution

	/**
	 * Sets the distance between this solution and a <code>SolutionSet</code>.
	 * The value is stored in <code>distanceToSolutionSet_</code>.
	 * 
	 * @param distance
	 *            The distance to a solutionSet.
	 */
	public void setDistanceToSolutionSet(double distance) {
		distanceToSolutionSet_ = distance;
	} // SetDistanceToSolutionSet

	/**
	 * Gets the distance from the solution to a <code>SolutionSet</code>.
	 * <b> REQUIRE </b>: this method has to be invoked after calling
	 * <code>setDistanceToPopulation</code>.
	 * 
	 * @return the distance to a specific solutionSet.
	 */
	public double getDistanceToSolutionSet() {
		return distanceToSolutionSet_;
	} // getDistanceToSolutionSet

	/**
	 * Sets the distance between the solution and its k-nearest neighbor in a
	 * <code>SolutionSet</code>. The value is stored in <code>kDistance_</code>.
	 * 
	 * @param distance
	 *            The distance to the k-nearest neighbor.
	 */
	public void setKDistance(double distance) {
		kDistance_ = distance;
	} // setKDistance

	/**
	 * Gets the distance from the solution to his k-nearest nighbor in a
	 * <code>SolutionSet</code>. Returns the value stored in
	 * <code>kDistance_</code>. <b> REQUIRE </b>: this method has to be invoked
	 * after calling <code>setKDistance</code>.
	 * 
	 * @return the distance to k-nearest neighbor.
	 */
	double getKDistance() {
		return kDistance_;
	} // getKDistance

	/**
	 * Sets the crowding distance of a solution in a <code>SolutionSet</code>.
	 * The value is stored in <code>crowdingDistance_</code>.
	 * 
	 * @param distance
	 *            The crowding distance of the solution.
	 */
	public void setCrowdingDistance(double distance) {
		crowdingDistance_ = distance;
	} // setCrowdingDistance

	/**
	 * Gets the crowding distance of the solution into a
	 * <code>SolutionSet</code>. Returns the value stored in
	 * <code>crowdingDistance_</code>. <b> REQUIRE </b>: this method has to be
	 * invoked after calling <code>setCrowdingDistance</code>.
	 * 
	 * @return the distance crowding distance of the solution.
	 */
	public double getCrowdingDistance() {
		return crowdingDistance_;
	} // getCrowdingDistance

	/**
	 * Sets the fitness of a solution. The value is stored in
	 * <code>fitness_</code>.
	 * 
	 * @param fitness
	 *            The fitness of the solution.
	 */
	public void setFitness(double fitness) {
		fitness_ = fitness;
	} // setFitness

	/**
	 * Gets the fitness of the solution. Returns the value of stored in the
	 * encodings.variable <code>fitness_</code>. <b> REQUIRE </b>: This method
	 * has to be invoked after calling <code>setFitness()</code>.
	 * 
	 * @return the fitness.
	 */
	public double getFitness() {
		return fitness_;
	} // getFitness

	/**
	 * Sets the value of the i-th objective.
	 * 
	 * @param i
	 *            The number identifying the objective.
	 * @param value
	 *            The value to be stored.
	 */
	public void setObjective(int i, double value) {
		objective_[i] = value;
	} // setObjective

	public void setObjectives(double[] values) {
		assert values.length == objective_.length;
		for (int i = 0; i < objective_.length; i++) {
			setObjective(i, values[i]);
		}
	}

	/**
	 * Returns the value of the i-th objective.
	 * 
	 * @param i
	 *            The value of the objective.
	 */
	public double getObjective(int i) {
		return objective_[i];
	} // getObjective

	public double[] getObjectives(){
		return objective_;
	}

	public double getObjectiveWeightedSum(){
		double sum = 0;
		for (int i = 0; i < getNumberOfObjectives(); i++){
			if (Double.isFinite(getObjective(i)))
				sum += getObjective(i);
		}
		return sum;
	}

	public double getObjectiveWeightedSum(double[] w){
		double sum = 0;
		if (w.length != getNumberOfObjectives()){
			System.out.println("Error: weighted vector length not match.");
			System.exit(0);
		}
		else{
			for (int i = 0; i < getNumberOfObjectives(); i++){
				if (Double.isFinite(getObjective(i)))
					sum += w[i] * getObjective(i);
			}
		}
		return sum;
	}

	public void resetObjective(){
		Arrays.fill(objective_, Double.POSITIVE_INFINITY);
	}

	/**
	 * Returns the number of objectives.
	 * 
	 * @return The number of objectives.
	 */
	public int getNumberOfObjectives() {
		if (objective_ == null)
			return 0;
		else
			return numberOfObjectives_;
	} // getNumberOfObjectives

	/**
	 * Returns the number of decision variables of the solution.
	 * 
	 * @return The number of decision variables.
	 */
	public int numberOfVariables() {
		return problemSet_.getMaxDimension();
	} // numberOfVariables

	/**
	 * Returns a string representing the solution.
	 * 
	 * @return The string.
	 */
	public String toString() {
		String aux = "";
		for (int i = 0; i < this.numberOfObjectives_; i++)
			aux = aux + this.getObjective(i) + " ";

		return aux;
	} // toString

	/**
	 * Returns the decision variables of the solution.
	 * 
	 * @return the <code>DecisionVariables</code> object representing the
	 *         decision variables of the solution.
	 */
	public Variable[] getDecisionVariables() {
		return variable_;
	} // getDecisionVariables

	public double getDecisionVariables(int idx) throws JMException {
		return variable_[idx].getValue();
	}

	public double[] getDecisionVariablesInDouble() throws JMException {
		double[] res = new double[variable_.length];
		Arrays.setAll(res, i -> {
			double arr = 0;
			try {
				arr = variable_[i].getValue();
			} catch (JMException e) {
				e.printStackTrace();
			}
			return arr;
		});
		return res;
	}

	/**
	 * Sets the decision variables for the solution.
	 * 
	 * @param variables
	 *            The <code>DecisionVariables</code> object representing the
	 *            decision variables of the solution.
	 */
	public void setDecisionVariables(Variable[] variables) {
		variable_ = variables;
	} // setDecisionVariables

	public void setDecisionVariables(double[] variables) throws JMException {
		if (variable_.length != variables.length) {
			System.out.println("Error, length not equal");
			return;
		}
		for (int i = 0; i < variable_.length; i++){
			variable_[i].setValue(variables[i]);
		}
	}

	public void setDecisionVariables(int idx, double value) throws JMException {
		variable_[idx].setValue(value);
	} // setDecisionVariables

	public ProblemSet getProblemSet() {
		return problemSet_;
	}

	/**
	 * Sets the rank of a solution.
	 * 
	 * @param value
	 *            The rank of the solution.
	 */
	public void setRank(int value) {
		this.rank_ = value;
	} // setRank

	/**
	 * Gets the rank of the solution. <b> REQUIRE </b>: This method has to be
	 * invoked after calling <code>setRank()</code>.
	 * 
	 * @return the rank of the solution.
	 */
	public int getRank() {
		return this.rank_;
	} // getRank

	/**
	 * Sets the overall constraints violated by the solution.
	 * 
	 * @param value
	 *            The overall constraints violated by the solution.
	 */
	public void setOverallConstraintViolation(double value) {
		this.overallConstraintViolation_ = value;
	} // setOverallConstraintViolation

	/**
	 * Gets the overall constraint violated by the solution. <b> REQUIRE </b>:
	 * This method has to be invoked after calling
	 * <code>overallConstraintViolation</code>.
	 * 
	 * @return the overall constraint violation by the solution.
	 */
	public double getOverallConstraintViolation() {
		return this.overallConstraintViolation_;
	} // getOverallConstraintViolation

	/**
	 * Sets the number of constraints violated by the solution.
	 * 
	 * @param value
	 *            The number of constraints violated by the solution.
	 */
	public void setNumberOfViolatedConstraint(int value) {
		this.numberOfViolatedConstraints_ = value;
	} // setNumberOfViolatedConstraint

	/**
	 * Gets the number of constraint violated by the solution. <b> REQUIRE </b>:
	 * This method has to be invoked after calling
	 * <code>setNumberOfViolatedConstraint</code>.
	 * 
	 * @return the number of constraints violated by the solution.
	 */
	public int getNumberOfViolatedConstraint() {
		return this.numberOfViolatedConstraints_;
	} // getNumberOfViolatedConstraint

	/**
	 * Sets the location of the solution into a solutionSet.
	 * 
	 * @param location
	 *            The location of the solution.
	 */
	public void setLocation(int location) {
		this.location_ = location;
	} // setLocation

	/**
	 * Gets the location of this solution in a <code>SolutionSet</code>.
	 * <b> REQUIRE </b>: This method has to be invoked after calling
	 * <code>setLocation</code>.
	 * 
	 * @return the location of the solution into a solutionSet
	 */
	public int getLocation() {
		return this.location_;
	} // getLocation

	/**
	 * Sets the type of the encodings.variable.
	 * 
	 * @param type
	 *            The type of the encodings.variable.
	 */
	// public void setType(String type) {
	// type_ = Class.forName("") ;
	// } // setType

	/**
	 * Sets the type of the encodings.variable.
	 * 
	 * @param type
	 *            The type of the encodings.variable.
	 */
	public void setType(SolutionType type) {
		type_ = type;
	} // setType

	/**
	 * Gets the type of the encodings.variable
	 * 
	 * @return the type of the encodings.variable
	 */
	public SolutionType getType() {
		return type_;
	} // getType

	/**
	 * Returns the aggregative value of the solution
	 * 
	 * @return The aggregative value.
	 */
	public double getAggregativeValue() {
		double value = 0.0;
		for (int i = 0; i < getNumberOfObjectives(); i++) {
			value += getObjective(i);
		}
		return value;
	} // getAggregativeValue

	public void setSkillFactor(int skillFactor) {
		skillFactor_ = skillFactor;
	}

	public int getSkillFactor() {
		return skillFactor_;
	}

	public double getUnifiedUpperLimit() {
		return problemSet_.getUnifiedUpperLimit();
	}

	public double getUnifiedLowerLimit() {
		return problemSet_.getUnifiedLowerLimit();
	}
	
	
	
	public double getGFunValue() {
		return gFunValue_;
	}
	
	public void setGFunValue(double gfv) {
		gFunValue_ = gfv;
	}


//add from moeac
    public void setNormalizedObjective(int j, double val) {
		normalizedObjective_[j] = val;
    }

	public double getNormalizedObjective(int j) {
		return normalizedObjective_[j];
	}

	public void setDistanceToIdealPoint(double distanceToIdealPoint) {
		this.distanceToIdealPoint = distanceToIdealPoint;
	}
	public double getDistanceToIdealPoint() {
		return distanceToIdealPoint;
	}

	private double sumValue;
	public void setSumValue(double sumValue) {
		this.sumValue = sumValue;
	}
	public double getSumValue() {
		return sumValue;
	}

	public double getUnitHyperplaneObjective(int i) {
		return unitHyperplaneObjective_[i];
	}
	public void setUnitHyperplaneObjective(int i,double unitValue) {
		unitHyperplaneObjective_[i] = unitValue;
	}

	private boolean remove;
	public boolean isRemove() {
		return remove;
	}
	public void setRemove(boolean remove) {
		this.remove = remove;
	}

	private int localDensity_;
	public int getLocalDensity_() {
		return localDensity_;
	}
	public void setLocalDensity_(int localDensity_) {
		this.localDensity_ = localDensity_;
	}

	private double gaussianLocalDensity_;

	public double getGaussianLocalDensity_() {
		return gaussianLocalDensity_;
	}

	public void setGaussianLocalDensity_(double gaussianLocalDensity_) {
		this.gaussianLocalDensity_ = gaussianLocalDensity_;
	}

	private double centerDistance_;

	public double getCenterDistance_() {
		return centerDistance_;
	}

	public void setCenterDistance_(double centerDistance_) {
		this.centerDistance_ = centerDistance_;
	}

	private double convergenceDistance_;

	public double getConvergenceDistance_() {
		return convergenceDistance_;
	}

	public void setConvergenceDistance_(double convergenceDistance_) {
		this.convergenceDistance_ = convergenceDistance_;
	}

	private double distanceToNadirPoint;

	public double getDistanceToNadirPoint() {
		return distanceToNadirPoint;
	}

	public void setDistanceToNadirPoint(double distanceToNadirPoint) {
		this.distanceToNadirPoint = distanceToNadirPoint;
	}

	private double[] translatedObjectives;
	private int rank1;
	private boolean isMarked;
	public int getRank1() {

		return this.rank1;
	}
	public void setRank1(int id1) {
		this.rank1 = id1;
	}

	public boolean getRank2() {
		return this.isMarked;
	}
	public void setRank2(boolean id2) {
		this.isMarked = id2;
	}

	public void setLabel(int target){
		this.label_ = target;
	}

	public int getLabel(){
		return this.label_;
	}

	/**
	 * Constructor
	 *
	 * @param numberOfObjectives
	 *            Number of objectives of the solution
	 *
	 *            This constructor is used mainly to read objective values from
	 *            a file to variables of a SolutionSet to apply quality
	 *            indicators
	 */


	/**
	 * Constructor.
	 *
	 * @param problem
	 *            The problem to solve
	 * @throws ClassNotFoundException
	 */
	public Solution(Problem problem) throws ClassNotFoundException {
//		problem_ = problem;
		type_ = problem.getSolutionType();
		numberOfObjectives_ = problem.getNumberOfObjectives();
		objective_ = new double[numberOfObjectives_];
		//clone_num=new int[numberOfObjectives_];

		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];
		remove = false;
		// Setting initial values
		fitness_ = 0.0;
		kDistance_ = 0.0;
		crowdingDistance_ = 0.0;
		localDensity_ = -1;
		gaussianLocalDensity_ = 0.0;
		centerDistance_ = 0.0;
		convergenceDistance_ = 0.0;
		distanceToSolutionSet_ = Double.POSITIVE_INFINITY;
		distanceToIdealPoint = Double.POSITIVE_INFINITY;
		distanceToNadirPoint = 0.0;
		sumValue = Double.POSITIVE_INFINITY;
		// <-
		// variable_ = problem.solutionType_.createVariables() ;
		variable_ = type_.createVariables();
		this.translatedObjectives = new double[this.numberOfObjectives_];
//		lamda_ = this;
		rank1 = -1;
		isMarked = false;
	} // Solution

	public Solution(Problem problem, int groupSize) throws ClassNotFoundException {

		type_ = problem.getSolutionType();
		numberOfObjectives_ = problem.getNumberOfObjectives();
		objective_ = new double[numberOfObjectives_];
		//clone_num=new int[numberOfObjectives_];

		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];
		remove = false;
		// Setting initial values
		fitness_ = 0.0;
		kDistance_ = 0.0;
		crowdingDistance_ = 0.0;
		localDensity_ = -1;
		gaussianLocalDensity_ = 0.0;
		centerDistance_ = 0.0;
		convergenceDistance_ = 0.0;
		distanceToSolutionSet_ = Double.POSITIVE_INFINITY;
		distanceToIdealPoint = Double.POSITIVE_INFINITY;
		distanceToNadirPoint = 0.0;
		sumValue = Double.POSITIVE_INFINITY;
		// variable_ = problem.solutionType_.createVariables() ;
		variable_ = type_.createVariables();

		rank1 = -1;
		isMarked = false;

	} // Solution

	public void setProblemSet_(ProblemSet problemSet) throws JMException {
		problemSet_ = problemSet;
		type_ = problemSet.getSolutionType();
		if (numberOfObjectives_ < problemSet_.getNumberOfObjs(0)){
			double[] newObjective = new double[problemSet_.getNumberOfObjs(0)];
			double[] newNormalizedObjective = new double[problemSet_.getNumberOfObjs(0)];
			double[] newUnitHyperplaneObjective = new double[problemSet_.getNumberOfObjs(0)];

			for (int i = 0; i < numberOfObjectives_; i++){
				newObjective[i] = objective_[i];
				newNormalizedObjective[i] = normalizedObjective_[i];
				newUnitHyperplaneObjective[i] = unitHyperplaneObjective_[i];
			}
			objective_ = newObjective;
			normalizedObjective_ = newNormalizedObjective;
			unitHyperplaneObjective_ = newUnitHyperplaneObjective;

			numberOfObjectives_ = problemSet_.getNumberOfObjs(0);

		}
		problemSet.get(0).evaluate(this);
	}

	static public Solution getNewSolution(Problem problem)
			throws ClassNotFoundException {
		return new Solution(problem);
	}

	/**
	 * Constructor
	 *
	 * @param problem
	 *            The problem to solve
	 */
	public Solution(Problem problem, Variable[] variables) {

		type_ = problem.getSolutionType();
		numberOfObjectives_ = problem.getNumberOfObjectives();
		objective_ = new double[numberOfObjectives_];

		normalizedObjective_ = new double[numberOfObjectives_];
		unitHyperplaneObjective_ = new double[numberOfObjectives_];

		remove = false;
		// Setting initial values
		fitness_ = 0.0;
		kDistance_ = 0.0;
		crowdingDistance_ = 0.0;
		localDensity_ = -1;
		gaussianLocalDensity_ = 0.0;
		centerDistance_ = 0.0;
		convergenceDistance_ = 0.0;
		distanceToSolutionSet_ = Double.POSITIVE_INFINITY;
		distanceToIdealPoint = Double.POSITIVE_INFINITY;
		distanceToNadirPoint = 0.0;
		sumValue = Double.POSITIVE_INFINITY;

		variable_ = variables;
		this.translatedObjectives = new double[this.numberOfObjectives_];

		rank1 = -1;
		isMarked = false;
	} // Constructor

	/**
	 * Copy constructor.
	 *
	 * @param solution
	 *            Solution to copy.
	 */

	private boolean marked_;
	public boolean isMarked_() {
		return this.marked_;
	} // isMarked

	/*Defined in the paper entiled "Towards Adaptive Knowledge Transfer in Multifactorial Evolutionary Computation"*/
	private int isTran_;
	private int cxFactor_;
	private Solution parent_;
	private double cObj_;

	public int getIsTran() {
		return isTran_;
	}

	public void setIsTran(int isTran) {
		this.isTran_ = isTran;
	}

	public double getcObj() {
		return cObj_;
	}

	public void setcObj(double cObj) {
		this.cObj_ = cObj;
	}

	public int getCxFactor() {
		return cxFactor_;
	}

	public void setCxFactor(int cxFactor) {
		this.cxFactor_ = cxFactor;
	}

	public Solution getParent() {
		return parent_;
	}

	public void setParent(Solution parent) {
		this.parent_ = parent;
	}

	public void setFlag(int flag){
		this.flag = flag;
	}

	public int getFlag(){
		return this.flag;
	}

	public void setFlag2(int flag) {
		this.flag2 = flag;
	}

	public int getFlag2() {
		return this.flag2;
	}

    public void setPdf(double density) {this.pdf_ = density; }

	public double getPdf() {return this.pdf_; }
} // Solution
