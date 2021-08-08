package etmo.metaheuristics.momfea;

import etmo.core.Algorithm;
import etmo.core.Operator;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PORanking;
import etmo.util.PseudoRandom;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.LocationComparator;
import etmo.util.logging.LogPopulation;

public class MOMFEA extends Algorithm {
	private int populationSize;
	
	private SolutionSet population;
	private SolutionSet offspringPopulation;
	private SolutionSet union;
	
	int evaluations;
	int maxEvaluations;
	
	Operator crossover;
	Operator mutation;
	Operator selection;
	
	double rmp;

	Distance distance = new Distance();
	
	public MOMFEA(ProblemSet problemSet) {
		super(problemSet);
	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		populationSize = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();
		rmp =  ((Double) getInputParameter("rmp")).doubleValue();

		crossover = operators_.get("crossover");
		mutation = operators_.get("mutation");
		selection = operators_.get("selection");

		evaluations = 0;

		initPopulation();
		while (evaluations < maxEvaluations) {
			createOffspringPopulation();
			getNextPopulation();
		}
		
		return population;
	}
	
	
	void initPopulation() throws JMException, ClassNotFoundException {
		population = new SolutionSet(populationSize);
		for (int i = 0; i < populationSize; i++) {
			Solution newSolution = new Solution(problemSet_);
			int id = PseudoRandom.randInt(0, problemSet_.size() - 1);
			problemSet_.get(id).evaluate(newSolution);
			problemSet_.get(id).evaluateConstraints(newSolution);
			evaluations++;
			
			newSolution.setSkillFactor(id);
			population.add(newSolution);
			
		} // for
		assignFitness(population);
	} // initPopulation
	
	
	
	void getNextPopulation() {
		union = population.union(offspringPopulation);
		
		//long tmpAFTime = System.currentTimeMillis();
		assignFitness(union);
		//long endTime = System.currentTimeMillis();
		union.sort(new LocationComparator());
		
		population.clear();
		
		for (int i = 0; i < populationSize; i++)
			population.add(union.get(i));
	}
	
	void createOffspringPopulation() throws JMException {
		offspringPopulation = new SolutionSet(populationSize);
		Solution[] parents = new Solution[2];
		for (int i = 0; i < (populationSize / 2); i++) {
			parents[0] = (Solution) selection.execute(population);
			parents[1] = (Solution) selection.execute(population);
			
			int[] sfs = new int[2];
			sfs[0] = parents[0].getSkillFactor();
			sfs[1] = parents[1].getSkillFactor();
			double rand = PseudoRandom.randDouble();

			Solution[] offSpring;
			if (sfs[0] == sfs[1] || rand < rmp) {
				offSpring = (Solution[]) crossover.execute(parents);
				mutation.execute(offSpring[0]);
				mutation.execute(offSpring[1]);
				
				int p0 = PseudoRandom.randInt(0, 1);
				int p1 = PseudoRandom.randInt(0, 1);
			//	int p1 = 1 - p0;
				offSpring[0].setSkillFactor(sfs[p0]);
				offSpring[1].setSkillFactor(sfs[p1]);
				
				resetObjectives(offSpring[0]);
				resetObjectives(offSpring[1]);
				problemSet_.get(sfs[p0]).evaluate(offSpring[0]);
				problemSet_.get(sfs[p1]).evaluate(offSpring[1]);
				problemSet_.get(sfs[p0]).evaluateConstraints(offSpring[0]);
				problemSet_.get(sfs[p1]).evaluateConstraints(offSpring[1]);
				evaluations += 2;

			} else {
				offSpring = new Solution[2];
				offSpring[0] = new Solution(parents[0]);
				offSpring[1] = new Solution(parents[1]);
	
				mutation.execute(offSpring[0]);
				mutation.execute(offSpring[1]);

				offSpring[0].setSkillFactor(sfs[0]);
				offSpring[1].setSkillFactor(sfs[1]);
	
				
				problemSet_.get(sfs[0]).evaluate(offSpring[0]);
				problemSet_.get(sfs[1]).evaluate(offSpring[1]);
				problemSet_.get(sfs[0]).evaluateConstraints(offSpring[0]);
				problemSet_.get(sfs[1]).evaluateConstraints(offSpring[1]);
				
				evaluations += 2;		
			}

			offspringPopulation.add(offSpring[0]);
			offspringPopulation.add(offSpring[1]);

		} // for
	}
	
	void assignFitness(SolutionSet pop) {
		for (int i = 0; i < pop.size(); i++)
			pop.get(i).setLocation(Integer.MAX_VALUE);
		for (int i = 0; i < problemSet_.size(); i++)
			rankSolutionOnTask(pop, i);
	}
	
	void rankSolutionOnTask(SolutionSet pop, int taskId) {
		int start = problemSet_.get(taskId).getStartObjPos();
		int end = problemSet_.get(taskId).getEndObjPos();
		
		boolean selec[] = new boolean[problemSet_.getTotalNumberOfObjs()];

		for (int i = 0; i < selec.length; i++) {
			if (i < start || i > end)
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

//	private void logPopulation(int eval){
//		int taskNum = problemSet_.size();
//		SolutionSet[] resPopulation = new SolutionSet[taskNum];
//		for (int i = 0; i < taskNum; i++)
//			resPopulation[i] = new SolutionSet();
//
//		for (int i = 0; i < population.size(); i++) {
//			Solution sol = population.get(i);
//
//			int pid = sol.getSkillFactor();
//
//			int start = problemSet_.get(pid).getStartObjPos();
//			int end = problemSet_.get(pid).getEndObjPos();
//
//			Solution newSolution = new Solution(end - start + 1);
//
//			for (int k = start; k <= end; k++)
//				newSolution.setObjective(k - start, sol.getObjective(k));
//
//			resPopulation[pid].add(newSolution);
//		}
//		for (int k = 0; k < taskNum; k++)
//			resPopulation[k].printObjectivesToFile("MOMFEA\\" + "MOMFEA_"+problemSet_.get(k).getNumberOfObjectives()+"Obj_"+
//					problemSet_.get(k).getName()+ "_" + problemSet_.get(k).getNumberOfVariables() + "D" + eval + ".txt");
//	}
}
