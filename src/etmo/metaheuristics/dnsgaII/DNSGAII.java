//  DNSGAII.java
//
//  Author:
//       Songbai Liu
//
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package etmo.metaheuristics.dnsgaII;

import java.util.ArrayList;
import java.util.List;

import etmo.core.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.Ranking;
import etmo.util.comparators.CrowdingComparator;

/**
 * Implementation of NSGA-II to solve dynamic multiple objective optimization problems. 
 */

public class DNSGAII extends DynamicAlgorithm {
	/**
	 * Constructor
	 * 
	 * @param problemSet
	 *            Problem to solve
	 */
	public DNSGAII(ProblemSet problemSet) {
		super(problemSet);
	//	System.out.println("sup: " + problemSet_.get(0).getHType());
	} // NSGAII

	/**
	 * Runs the NSGA-II algorithm.
	 * 
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 *         solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	public List<SolutionSet> execute() throws JMException, ClassNotFoundException {
		int populationSize;
		int maxEvaluations;
		int evaluations;
		int generations;

		
		/*
		 * fc indicates the frequency of change
		 * nc number of changes in each run
		*/
		int fc = 20;
		int nc = 30;

		SolutionSet population;
		List<SolutionSet> dynamicPopulationSets;
		SolutionSet offspringPopulation;
		SolutionSet union;

		Operator mutationOperator;
		Operator crossoverOperator;
		Operator selectionOperator;

		Distance distance = new Distance();

		// Read the parameters
		populationSize = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();

		// Initialize the variables
		population = new SolutionSet(populationSize);
		dynamicPopulationSets = new ArrayList<SolutionSet>();
		evaluations = 0;
		generations = 1;


		// Read the operators
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		selectionOperator = operators_.get("selection");

		// Create the initial solutionSet
		Solution newSolution;
		for (int i = 0; i < populationSize; i++) {
			newSolution = new Solution(problemSet_);
			problemSet_.get(0).dynamicEvaluate(newSolution,generations);
			evaluations++;
			population.add(newSolution);
		} // for
		Ranking ranking;
		// Generations
		while (evaluations < maxEvaluations) {
			
			if((generations) % fc == 0){
				ranking = new Ranking(population);
				dynamicPopulationSets.add(ranking.getSubfront(0));
			}
			generations++;

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
					problemSet_.get(0).dynamicEvaluate(offSpring[0],generations);
					problemSet_.get(0).dynamicEvaluate(offSpring[1],generations);
					offspringPopulation.add(offSpring[0]);
					offspringPopulation.add(offSpring[1]);
					evaluations += 2;
				} // if
			} // for

			// Create the solutionSet union of solutionSet and offSpring
			union = ((SolutionSet) population).union(offspringPopulation);

			// Ranking the union
			ranking = new Ranking(union);

			int remain = populationSize;
			int index = 0;
			SolutionSet front = null;
			population.clear();

			// Obtain the next front
			front = ranking.getSubfront(index);

			while ((remain > 0) && (remain >= front.size())) {
				// Assign crowding distance to individuals
				distance.crowdingDistanceAssignment(front, problemSet_.get(0).getNumberOfObjectives());
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
				distance.crowdingDistanceAssignment(front, problemSet_.get(0).getNumberOfObjectives());
				front.sort(new CrowdingComparator());
				for (int k = 0; k < remain; k++) {
					population.add(front.get(k));
				} // for

				remain = 0;
			} // if
			
		} // while

		if(dynamicPopulationSets.size() != nc){
			System.out.println("The real number of changes is " + dynamicPopulationSets.size() +
			" which is not the same with the preset " + nc);
		}
		return dynamicPopulationSets;
		
	} // execute
} // DNSGA-II
