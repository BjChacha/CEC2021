/*
 * NonUniformGenerator.java
 *
 * @author Antonio J. Nebro
 * @version 1.0
 *
 * This class returns a solution after applying SBX and Polynomial mutation
 */
package etmo.util.offspring;

import etmo.core.Operator;
import etmo.core.Solution;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.util.JMException;

import java.util.HashMap;

public class NonUniformMutationOffspring extends Offspring {

	public Operator mutation_;
	private Operator selection_;

	private double mutationProbatility_;
	private double perturbation_;
	private int maxIterations_;

	public NonUniformMutationOffspring(double mutationProbability, double perturbation, int maxIterations)
			throws JMException {
		HashMap parameters; // Operator parameters
		parameters = new HashMap();
		parameters.put("probability", mutationProbatility_ = mutationProbability);
		parameters.put("perturbation", perturbation_ = perturbation);
		parameters.put("maxIterations", maxIterations_ = maxIterations);
		mutation_ = MutationFactory.getMutationOperator("NonUniformMutation", parameters);

		selection_ = SelectionFactory.getSelectionOperator("BinaryTournament", null);
		id_ = "NonUniformMutation";
	}

	public Solution getOffspring(Solution solution) {
		Solution res = new Solution(solution);
		try {
			mutation_.execute(res);
		} catch (JMException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return res;
	}

	public String configuration() {
		String result = "-----\n";
		result += "Operator: " + id_ + "\n";
		result += "Probability: " + mutationProbatility_ + "\n";
		result += "MaxIterations: " + maxIterations_ + "\n";
		result += "Perturbation: " + perturbation_;

		return result;
	}
} // PolynomialOffspringGenerator
