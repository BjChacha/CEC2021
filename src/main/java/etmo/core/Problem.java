//  Problem.java
//
//  Authors:
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

import java.io.Serializable;

import etmo.util.JMException;

/**
 * Abstract class representing a multiobjective optimization problem
 */
public abstract class Problem implements Serializable {

	/**
	 * Defines the default precision of binary-coded variables
	 */
	private final static int DEFAULT_PRECISSION = 16;

	/**
	 * Stores the number of variables of the problem
	 */
	protected int numberOfVariables_;

	/**
	 * Stores the number of objectives of the problem
	 */
	protected int numberOfObjectives_;

	/**
	 * Stores the number of constraints of the problem
	 */
	protected int numberOfConstraints_;

	/**
	 * Stores the problem name
	 */
	protected String problemName_;

	protected int problemIndex_ = 0;

	
	protected int  startObjPos_;
	protected int endObjPos_;
	protected String hType_; //The type of the shape function

	/**
	 * Stores the type of the solutions of the problem
	 */
	protected SolutionType solutionType_;

	/**
	 * Stores the lower bound values for each encodings.variable (only if
	 * needed)
	 */
	protected double[] lowerLimit_;

	/**
	 * Stores the upper bound values for each encodings.variable (only if
	 * needed)
	 */
	protected double[] upperLimit_;

	/**
	 * Stores the number of bits used by binary-coded variables (e.g.,
	 * BinaryReal variables). By default, they are initialized to
	 * DEFAULT_PRECISION)
	 */
	private int[] precision_;

	/**
	 * Stores the length of each encodings.variable when applicable (e.g.,
	 * Binary and Permutation variables)
	 */
	protected int[] length_;

	/**
	 * Stores the type of each encodings.variable
	 */
	// public Class [] variableType_;

	
	
	// for g function
	protected double[] shiftValues_;
	protected double[][] rotationMatrix_;
	
	
	
	/**
	 * Constructor.
	 */
	public Problem() {
		solutionType_ = null;
	} // Problem

	/**
	 * Constructor.
	 */
	public Problem(SolutionType solutionType) {
		solutionType_ = solutionType;
	} // Problem

	/**
	 * Gets the number of decision variables of the problem.
	 * 
	 * @return the number of decision variables.
	 */
	public int getNumberOfVariables() {
		return numberOfVariables_;
	} // getNumberOfVariables

	/**
	 * Sets the number of decision variables of the problem.
	 */
	public void setNumberOfVariables(int numberOfVariables) {
		numberOfVariables_ = numberOfVariables;
	} // getNumberOfVariables

	/**
	 * Gets the the number of objectives of the problem.
	 * 
	 * @return the number of objectives.
	 */
	public int getNumberOfObjectives() {
		return numberOfObjectives_;
	} // getNumberOfObjectives

	/**
	 * Gets the lower bound of the ith encodings.variable of the problem.
	 * 
	 * @param i
	 *            The index of the encodings.variable.
	 * @return The lower bound.
	 */
	public double getLowerLimit(int i) {
		return lowerLimit_[i];
	} // getLowerLimit

	/**
	 * Gets the upper bound of the ith encodings.variable of the problem.
	 * 
	 * @param i
	 *            The index of the encodings.variable.
	 * @return The upper bound.
	 */
	public double getUpperLimit(int i) {
		return upperLimit_[i];
	} // getUpperLimit

	/**
	 * Evaluates a <code>Solution</code> object.
	 * 
	 * @param solution
	 *            The <code>Solution</code> to evaluate.
	 */
	public abstract void evaluate(Solution solution) throws JMException;
	
	
	/**
	 * Evaluates a <code>Solution</code> object.
	 * 
	 * @param solution
	 *            The <code>Solution</code> to evaluate.
	 */
	public abstract void dynamicEvaluate(Solution solution, int currentGeneration) throws JMException;

	/**
	 * Gets the number of side constraints in the problem.
	 * 
	 * @return the number of constraints.
	 */
	public int getNumberOfConstraints() {
		return numberOfConstraints_;
	} // getNumberOfConstraints

	/**
	 * Evaluates the overall constraint violation of a <code>Solution</code>
	 * object.
	 * 
	 * @param solution
	 *            The <code>Solution</code> to evaluate.
	 */
	public void evaluateConstraints(Solution solution) throws JMException {
		// The default behavior is to do nothing. Only constrained problems have
		// to
		// re-define this method
	} // evaluateConstraints

	/**
	 * Returns the number of bits that must be used to encode binary-real
	 * variables
	 * 
	 * @return the number of bits.
	 */
	public int getPrecision(int var) {
		return precision_[var];
	} // getPrecision

	/**
	 * Returns array containing the number of bits that must be used to encode
	 * binary-real variables.
	 * 
	 * @return the number of bits.
	 */
	public int[] getPrecision() {
		return precision_;
	} // getPrecision

	/**
	 * Sets the array containing the number of bits that must be used to encode
	 * binary-real variables.
	 * 
	 * @param precision
	 *            The array
	 */
	public void setPrecision(int[] precision) {
		precision_ = precision;
	} // getPrecision

	/**
	 * Returns the length of the encodings.variable.
	 * 
	 * @return the encodings.variable length.
	 */
	public int getLength(int var) {
		if (length_ == null)
			return DEFAULT_PRECISSION;
		return length_[var];
	} // getLength

	/**
	 * Sets the type of the variables of the problem.
	 * 
	 * @param type
	 *            The type of the variables
	 */
	public void setSolutionType(SolutionType type) {
		solutionType_ = type;
	} // setSolutionType

	/**
	 * Returns the type of the variables of the problem.
	 * 
	 * @return type of the variables of the problem.
	 */
	public SolutionType getSolutionType() {
		return solutionType_;
	} // getSolutionType

	/**
	 * Returns the problem name
	 * 
	 * @return The problem name
	 */
	public String getName() {
		return problemName_;
	}
	
	
	public void setName(String name) {
		this.problemName_ = name;
	}
	/**
	 * Returns the number of bits of the solutions of the problem
	 * 
	 * @return The number of bits solutions of the problem
	 */
	public int getNumberOfBits() {
		int result = 0;
		for (int var = 0; var < numberOfVariables_; var++) {
			result += getLength(var);
		}
		return result;
	} // getNumberOfBits();

	public void setProblemIndex(int index) {
		problemIndex_ = index;
	}

	public int getProblemIndex() {
		return problemIndex_;
	}
	
	

	public void setStartObjPos(int pos) {
		startObjPos_ = pos;
	}
	
	
	public void setEndObjPos(int pos) {
		endObjPos_ = pos;
	}
	
	
	public int getStartObjPos() {
		return startObjPos_;
	}
	
	public int getEndObjPos() {
		return endObjPos_;
	}
	

	public double[] scaleVariables(Solution solution) throws JMException {
		Variable[] decisionVariables = solution.getDecisionVariables();
		double[] x = new double[numberOfVariables_];

/*		System.out.println("...........");
		for (int i = 0; i < numberOfVariables_; i++) {
			x[i] = decisionVariables[i].getValue();
			System.out.print(x[i] + " ");
		}
		System.out.println();*/
		
		
	
		for (int i = 0; i < numberOfVariables_; i++) 
			x[i] = decisionVariables[i].getValue();


		for (int i = 0; i < numberOfVariables_; i++) {
			double sl = decisionVariables[i].getLowerBound();
			double su = decisionVariables[i].getUpperBound();
			double pl = lowerLimit_[i];
			double pu = upperLimit_[i];
			
/*			System.out.println("pl: " + pl);
			System.out.println("pu: " + pu);*/

			x[i] = ((x[i] - sl) * (pu - pl)) / (su - sl) + pl;

		}
		
/*		System.out.println("&&&&&&&&&");
		for (int i = 0; i < numberOfVariables_; i++) {
			System.out.print(x[i] + " ");
		}
		System.out.println();*/
		
		return x;
	}
	
	

	public void setRotationMatrix(double[][] matrix) {
		this.rotationMatrix_ = matrix;
	}

	public double[][] getRotationMatrix() {
		return rotationMatrix_;
	}

	public void setShiftValues(double[] shiftValues) {
		shiftValues_ = shiftValues;
	}

	public double[] getShiftValues() {
		return shiftValues_;
	}

	protected double[] transformVariables(double x[]) {
		shiftVariables(x);
		return rotateVariables(x);
	}

	protected void shiftVariables(double x[]) {
		for (int i = 0; i < x.length; i++)
			x[i] -= shiftValues_[i];
	}

	protected double[] rotateVariables(double x[]) {
		int len = x.length;
		double res[] = new double[len];

		for (int i = 0; i < len; i++) {
			double[] y = rotationMatrix_[i];

			double sum = 0;
			for (int j = 0; j < len; j++)
				sum += x[j] * y[j];
			res[i] = sum;
		}

		return res;
	}
	
	public void setHType(String hType) {
		hType_ = hType;
	}

	public String getHType() {
		return hType_;
	}

} // Problem
