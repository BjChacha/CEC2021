package etmo.metaheuristics.maoeac;

import etmo.core.*;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.comparators.SumValueComparator;
import etmo.util.ranking.NondominatedRanking;
import etmo.util.ranking.Ranking;
import etmo.util.logging.LogPopulation;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

//mating selection with the selection restricted in the same cluster + integrated learning in recombination step
public class MaOEAC extends Algorithm{
	private SolutionSet population_;
	private SolutionSet offspringPopulation_;
	private SolutionSet union_;

	private int populationSize_;
	int generations_;
	int maxGenerations_;

	Operator selection_;
	Operator crossover_;
	Operator mutation_;
	Operator learning_;

	private final double[] zideal_; //ideal point
	private final double[] znadir_;//Nadir point

	public MaOEAC(ProblemSet problemSet) {
		super(problemSet);
		zideal_ = new double[problemSet.get(0).getNumberOfObjectives()];
		znadir_ = new double[problemSet.get(0).getNumberOfObjectives()];
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {

		/*
		 * step1: Basic Setting of this Algorithm
		 */
		baiscSetting();
		/*
		 * step2: Initialize the Population
		 */
		initPopulation();
		/*
		 * Enter the main loop��into the process of evolution
		 */
		while (generations_ < maxGenerations_) {
			/*
			 * step3 and step4:Mating Selection and Recombination to generate Offspring Populations
			 */
			generateOffspringPopulation();
			/*
			 * step5:Environmental Selection
			 */
		    environmentalSelection();

		    generations_++;

			if (generations_ % 20 == 0){
				LogPopulation.LogPopulation("MaOEAC", population_, problemSet_, generations_*populationSize_, false);
			}
		}

		Ranking ranking = new NondominatedRanking(population_);
		return ranking.getSubfront(0);
	}

	/*
	 * step1: Basic Setting of this Algorithm
	 */
	public void baiscSetting(){
		generations_ = 0;
		maxGenerations_ = ((Integer) this.getInputParameter("maxGenerations")).intValue();
		populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
		mutation_  = operators_.get("mutation");
		crossover_ = operators_.get("crossover");
		selection_ = operators_.get("selection");

		learning_ = operators_.get("learning");

	}

	/*
	 * step2: Initialize the Population
	 */
	public void initPopulation() throws JMException, ClassNotFoundException {

		population_ = new SolutionSet(populationSize_);

		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problemSet_);
			problemSet_.get(0).evaluate(newSolution);
			problemSet_.get(0).evaluateConstraints(newSolution);
			population_.add(newSolution);
		} // for
		estimateIdealPoint(population_);
		estimateNadirPoint(population_);
		normalizationObjective(population_);
		computeDistanceToIdealPoint(population_);
	} // initPopulation


	/*
	 * step3 and step4:Mating Selection and Recombination to generate Offspring Populations
	 */
	public void generateOffspringPopulation() throws JMException{
		offspringPopulation_ = new SolutionSet(populationSize_);
		SolutionSet[] offspringSolutionSets = new PartitionalSolutionSet(population_,problemSet_.get(0).getNumberOfObjectives()).partitional();
		Solution[] gbests = new Solution[problemSet_.get(0).getNumberOfObjectives()];
		//population_.sort(new SumValueComparator());
		for(int i=0;i<problemSet_.get(0).getNumberOfObjectives();i++ ){
			offspringSolutionSets[i].sort(new SumValueComparator());
			gbests[i] = offspringSolutionSets[i].get(0);
		}
		for(int i=0;i<problemSet_.get(0).getNumberOfObjectives();i++ ){
			Solution[] parents = new Solution[2];
			for(int j=0;j < offspringSolutionSets[i].size();j++){
				double rd0 = PseudoRandom.randDouble();
				if(rd0 < 0.2){
					parents = (Solution[]) selection_.execute(population_);
				}else{
					parents = (Solution[]) selection_.execute(offspringSolutionSets[i]);
				}
				Solution[] offSpring = (Solution[]) crossover_
						.execute(parents);
				/*if(rd0 < 0.2){
					mutation_.execute(offSpring[0]);
				}else if(rd0 < 0.8 && rd0 > 0.2){
					Solution[] objectSolutionSet = new Solution[4];
					objectSolutionSet[0] = offSpring[0];
					//double rd = PseudoRandom.randDouble();
					if(parents[0].getSumValue() < parents[1].getSumValue()){
						objectSolutionSet[1] = parents[0];
					}else{
						objectSolutionSet[1] = parents[1];
					}
					int rnd = PseudoRandom.randInt(0, offspringSolutionSets[i].size()/5);
					objectSolutionSet[2] = offspringSolutionSets[i].get(rnd);

					int red = PseudoRandom.randInt(0, problem_.getNumberOfObjectives()-1);
					objectSolutionSet[3] = gbests[red];
					learning_.execute(objectSolutionSet);
				}else{
					Solution[] objectSolutionSet = new Solution[4];
					objectSolutionSet[0] = offSpring[0];
					//double rd = PseudoRandom.randDouble();
					if(parents[0].getSumValue() < parents[1].getSumValue()){
						objectSolutionSet[1] = parents[0];
					}else{
						objectSolutionSet[1] = parents[1];
					}
					int rnd = PseudoRandom.randInt(0, offspringSolutionSets[i].size()/5);
					objectSolutionSet[2] = offspringSolutionSets[i].get(rnd);

					int red = PseudoRandom.randInt(0, problem_.getNumberOfObjectives()-1);
					objectSolutionSet[3] = gbests[red];
					learning_.execute(objectSolutionSet);
					mutation_.execute(offSpring[0]);
				}*/
				mutation_.execute(offSpring[0]);

				problemSet_.get(0).evaluate(offSpring[0]);
				problemSet_.get(0).evaluateConstraints(offSpring[0]);
				offspringPopulation_.add(offSpring[0]);
			}
		}
	}

	/*
	 * step5:Environmental Selection
	 */
	public void environmentalSelection(){
		/*
		 * step5.1:Combine the Population and the Offspring Population
		 */
		union_ = population_.union(offspringPopulation_);
		/*
		 * step5.2:Normalization the Combined Population
		 */
		/*Ranking nodominatedRanking = new NondominatedRanking(union_);
		SolutionSet front = nodominatedRanking.getSubfront(0);
		estimateIdealPoint(front);
		estimateNadirPoint(front);*/
		estimateIdealPoint(union_);
		estimateNadirPoint(union_);
		normalizationObjective(union_);

		/*
		 * step5.3:Compute the Convergence Distance of each Solution
		 */
		computeDistanceToIdealPoint(union_);
		population_.clear();


		/*
		 * step5.4:Partitional Clustering Based K-means
		 */
		SolutionSet[] solutionSets = new PartitionalSolutionSet(union_,problemSet_.get(0).getNumberOfObjectives()).partitional();
		for(int k=0;k<problemSet_.get(0).getNumberOfObjectives();k++){
			// temporary fix: calculate populationSize.
			populationSize_ = solutionSets[k].size() * problemSet_.get(0).getNumberOfObjectives() / 2;

			if(solutionSets[k].size() != 2*populationSize_/problemSet_.get(0).getNumberOfObjectives()){
				System.out.println("solutionSets["+k+"] = "+ solutionSets[k].size());
				System.exit(0);
			}
			SolutionSet st = getStSolutionSet(solutionSets[k],populationSize_/(problemSet_.get(0).getNumberOfObjectives()));
			List<SolutionSet> list = new <SolutionSet>ArrayList();
			for(int i=0;i<st.size();i++){
				SolutionSet sols = new SolutionSet();
				sols.add(st.get(i));
				list.add(sols);
			}
			if(list.size() < populationSize_/(problemSet_.get(0).getNumberOfObjectives())){
				System.out.println("ListSize4 = "+list.size());
				System.exit(0);
			}
			/*
			 * step5.5:Agglomerative Hierarchical Clustering Based Average-Link Method
			 * and K-Cluster Stopping Condition
			 */
			list = new HierarchicalClustering1(list).clusteringAnalysis(populationSize_/(problemSet_.get(0).getNumberOfObjectives()));
			if(list.size() != populationSize_/(problemSet_.get(0).getNumberOfObjectives())){
				System.out.println("ListSize1 = "+list.size());
				System.exit(0);
			}

			/*
			 * Step5.6:Choose the Best Solution in each Cluster and Stored it into the next Generation
			 */
			bestSolutionSelection(list,k);
		}
	}

	public double computePBIFitness(Solution s1,Solution s2){
		double fitness = 0.0;
		double d1,d2,norm;
		int objectiveSize = s1.getNumberOfObjectives();
		d1 = d2 = norm = 0.0;
		for(int i=0; i<objectiveSize; i++){
			d1 += s1.getNormalizedObjective(i)*s2.getNormalizedObjective(i);
			norm += s1.getNormalizedObjective(i)*s1.getNormalizedObjective(i);
		}
		norm = Math.sqrt(norm);
		d1 = Math.abs(d1)/norm;
		for(int j=0; j<objectiveSize; j++){
			d2 += (s2.getNormalizedObjective(j)-d1*(s1.getNormalizedObjective(j)/norm))
					*(s2.getNormalizedObjective(j)-d1*(s1.getNormalizedObjective(j)/norm));
		}
		d2 = Math.sqrt(d2);
		fitness = d1 + 2.0*d2;
		return fitness;
	}

	/*
	 * Estimate the Ideal Point
	 */
	public void estimateIdealPoint(SolutionSet solutionSet){
		for(int i=0; i<problemSet_.get(0).getNumberOfObjectives();i++){
			zideal_[i] = 1.0e+30;
			for(int j=0; j<solutionSet.size();j++){
				if(solutionSet.get(j).getObjective(i) < zideal_[i]){
					zideal_[i] = solutionSet.get(j).getObjective(i);
				}
			}

		}
	}

	/*
	 * Estimate the Nadir Point
	 */
    public void estimateNadirPoint(SolutionSet solutionSet){
    	for(int i=0; i<problemSet_.get(0).getNumberOfObjectives();i++){
			znadir_[i] = -1.0e+30;
			for(int j=0; j<solutionSet.size();j++){
				if(solutionSet.get(j).getObjective(i) > znadir_[i]){
					znadir_[i] = solutionSet.get(j).getObjective(i);
				}
			}

		}
	}

    /*
     * Normalization
     */
	public void normalizationObjective(SolutionSet solutionSet){
		for(int i=0; i<solutionSet.size(); i++){
			Solution sol = solutionSet.get(i);

			for(int j=0; j<problemSet_.get(0).getNumberOfObjectives(); j++){
				double val = 0.0;
				val = (sol.getObjective(j) - zideal_[j])/(znadir_[j]-zideal_[j] + 1e-30);
				//val = (sol.getObjective(j) - zideal_[j]);
				sol.setNormalizedObjective(j, val);
			}
		}
	}

	 /*
     * Compute the Convergence Distance of each Solutions Which use the distance of
     * each solution to the Ideal Point
     */
    public void computeDistanceToIdealPoint(SolutionSet solutionSet){
    	for(int i=0; i<solutionSet.size(); i++){
    		Solution sol = solutionSet.get(i);
    		double normDistance = 0.0;
    		double sumValue = 0.0;
    		for(int j=0; j<problemSet_.get(0).getNumberOfObjectives(); j++){
    			normDistance += sol.getNormalizedObjective(j) * sol.getNormalizedObjective(j);
    			sumValue +=  sol.getNormalizedObjective(j);
    		}
    		normDistance = Math.sqrt(normDistance);

    		sol.setDistanceToIdealPoint(normDistance);
    		sol.setSumValue(sumValue);
    	}
    }

    public double computeDistance(Solution solution){
    	double p = 1.2;
    	double normDistance = 0.0;
		for(int j=0; j<problemSet_.get(0).getNumberOfObjectives(); j++){
			normDistance += Math.pow(solution.getNormalizedObjective(j),p);
		}
		normDistance = Math.pow(normDistance,1.0/p);
    	return normDistance;
    }

	 /*
     * Compute the angle value between Solution1 and Solution2
     */
	public double computeAngle(Solution s1, Solution s2){
		double angle = 0.0;
		double distanceToidealPoint1 = s1.getDistanceToIdealPoint();
		double distanceToidealPoint2 = s2.getDistanceToIdealPoint();
		double innerProduc = 0.0;
		for(int i=0; i<problemSet_.get(0).getNumberOfObjectives(); i++){
			innerProduc += s1.getNormalizedObjective(i) * s2.getNormalizedObjective(i);
		}
		double value = innerProduc/(distanceToidealPoint1*distanceToidealPoint2);
		if(value > 1.0){
			value = 1.0;
		}
		angle = Math.acos(Math.abs(value));
		//System.out.println(Math.abs(innerProduc/(distanceToidealPoint1*distanceToidealPoint2)));
		return angle;
	}//computeAngle

	public double weightSumValue(Solution solution){
		double value = 0.0;
		for(int i=0; i<problemSet_.get(0).getNumberOfObjectives();i++){
			value += solution.getNormalizedObjective(i);
		}
		return value;
	}
	public double weightSumValue(Solution solution1,Solution solution2){
		double value = 0.0;
		double[] lamda = new double[problemSet_.get(0).getNumberOfObjectives()];
		double sum = weightSumValue(solution1);
		for(int i=0; i<problemSet_.get(0).getNumberOfObjectives();i++){
			lamda[i] = solution1.getNormalizedObjective(i)/sum;
			value += lamda[i]*solution2.getNormalizedObjective(i);
		}
		return value;
	}

	public double weightSumValue(int k,Solution solution2){
		double value = 0.0;
		for(int i=0; i<problemSet_.get(0).getNumberOfObjectives();i++){
			if(i==k){
				value += (1.0/(2*Math.sqrt(problemSet_.get(0).getNumberOfObjectives()))+0.5)*solution2.getNormalizedObjective(i);
			}else{
				value += (1.0/(2*Math.sqrt(problemSet_.get(0).getNumberOfObjectives())))*solution2.getNormalizedObjective(i);
			}
		}
		return value;
	}

	public double computeChebyshev(Solution individua2, Solution individua1){
		double fitness;
		fitness = 0.0;

		double maxFun = -1.0e+30;

		for (int n = 0; n < problemSet_.get(0).getNumberOfObjectives(); n++) {
			double diff = Math.abs(individua1.getNormalizedObjective(n));

			double feval;
			if (individua2.getNormalizedObjective(n) == 0) {
				feval = diff / 0.000001;
			} else {
				feval = diff / individua2.getNormalizedObjective(n);
			}
			if (feval > maxFun) {
				maxFun = feval;
			}
		} // for

		fitness = maxFun;
		return fitness;
	}

	public double computeASFFitness(Solution s1,Solution s2){
		double fitness = -1.0e+30;
		int objectiveSize = s1.getNumberOfObjectives();
		double[] lambda = new double[objectiveSize];
		double sumValue = 0.0;
		for(int i=0; i<objectiveSize; i++){
			sumValue += s1.getNormalizedObjective(i);
		}
		for(int j=0; j<objectiveSize; j++){
			lambda[j] = s1.getNormalizedObjective(j)/sumValue;
			if(lambda[j] == 0){
				lambda[j] = 0.000001;
			}
		}
		for(int k=0; k<objectiveSize; k++){
			double sb = s2.getNormalizedObjective(k)/lambda[k];
			if(fitness < sb){
				fitness = sb;
			}
		}
		return fitness;
	}

   public SolutionSet getStSolutionSet(SolutionSet ss,int size) {

		SolutionSet sets = new SolutionSet();
		Ranking ranking = new NondominatedRanking(ss);

		int remain = size;
		int index = 0;
		SolutionSet front = null;
		SolutionSet mgPopulation = new SolutionSet();
		front = ranking.getSubfront(index);
		while ((remain > 0) && (remain >= front.size())) {

			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			} // for

			// Decrement remain
			remain = remain - front.size();

			// Obtain the next front
			index++;
			if (remain > 0) {
				front = ranking.getSubfront(index);
			} // if
		}
		if (remain > 0) { // front contains individuals to insert
			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			}
		}

		sets = mgPopulation;

		return sets;
	 }
   public void bestSolutionSelection(List<SolutionSet> list,int k) {
	   double minClustering2Axis = 1.0e+30;
	   int minClustering2AxisID = -1;
	   for (int i = 0; i < list.size(); i++) {
		   SolutionSet sols = list.get(i);
		   if (sols.size() == 0) {
			   System.out.println("SolsSize1 = " + sols.size());
			   System.exit(0);
		   }

		   double angle1 = Math.acos(Math.abs(sols.getCentroidVector().getNormalizedObjective(k) / sols.getCentroidVector().getDistanceToIdealPoint()));
		   // TODO: for now
		   if (Double.isNaN(angle1)) {
			   angle1 = 0.0;
			   System.out.println("Warning: angle being Nan.");
		   }

		   if (angle1 < minClustering2Axis) {
			   minClustering2Axis = angle1;
			   minClustering2AxisID = i;
		   }//if
	   }//for
	   double minSolution2Axis = 1.0e+30;
	   int minSolution2AxisID = -1;
	   for (int j = 0; j < list.get(minClustering2AxisID).size(); j++) {
		   Solution sol = list.get(minClustering2AxisID).get(j);
		   double ang = Math.acos(list.get(minClustering2AxisID).get(j).getNormalizedObjective(k) / list.get(minClustering2AxisID).get(j).getDistanceToIdealPoint());
		   if (ang < minSolution2Axis) {
			   minSolution2Axis = ang;
			   minSolution2AxisID = j;
		   }
	   }//for
	   // 保险措施：老是会等于-1
	   minSolution2AxisID = Math.max(0, minSolution2AxisID);

	   population_.add(list.get(minClustering2AxisID).get(minSolution2AxisID));
	   list.remove(minClustering2AxisID);

	   double min2CenterLine = 1.0e+30;
	   int min2CenterLineId = -1;
	   for (int i = 0; i < list.size(); i++) {
		   SolutionSet sols = list.get(i);
			/*if(sols.size() != 0){
				System.out.println("SolsSize1 = "+sols.size());
				//System.exit(0);
			}*/
		   double sumValue = 0.0;
		   for (int j = 0; j < problemSet_.get(0).getNumberOfObjectives(); j++) {
			   sumValue += sols.getCentroidVector().getNormalizedObjective(j) * 1.0;
		   }
		   //System.out.println("value = "+sumValue/(sols.getCentroidVector().getDistanceToIdealPoint()));
		   //norm2 = Math.sqrt(norm2);
		   double angle2 = Math.acos(sumValue / (sols.getCentroidVector().getDistanceToIdealPoint() * Math.sqrt(problemSet_.get(0).getNumberOfObjectives())));
		   //System.out.println(angle2);
		   if (angle2 < min2CenterLine) {
			   min2CenterLine = angle2;
			   min2CenterLineId = i;
		   }
	   }

	   // 保险措施：老是会等于-1
	   min2CenterLineId = Math.max(0, min2CenterLineId);

	   //System.out.println(min2CenterLineId);
	   double minS2CenterLine = 1.0e+30;
	   int minId = -1;
	   for (int i = 0; i < list.get(min2CenterLineId).size(); i++) {
		   Solution sol = list.get(min2CenterLineId).get(i);
		   double sumValue = 0.0;
		   for (int j = 0; j < problemSet_.get(0).getNumberOfObjectives(); j++) {
			   sumValue += sol.getNormalizedObjective(j);
		   }
		   double ang = Math.acos(Math.abs(sumValue / (sol.getDistanceToIdealPoint() * Math.sqrt(problemSet_.get(0).getNumberOfObjectives()))));
		   if (ang < minS2CenterLine) {
			   minS2CenterLine = ang;
			   minId = i;
		   }
	   }

	   // 保险措施：老是会等于-1
	   minId = Math.max(0, minId);

	   if (PseudoRandom.randDouble() < 0.5) {
		   population_.add(list.get(min2CenterLineId).get(minId));
		   list.remove(min2CenterLineId);
	   }



		/*if(list.size() != populationSize_/(problemSet_.get(0).getNumberOfObjectives()) - 2){
			System.out.println("ListSize3 = "+list.size());
			System.exit(0);
		}*/
	   double delta_ = 1.0;
	   Iterator<SolutionSet> it = list.iterator();
	   while (it.hasNext()) {
		   int type = -1;
		   double rnd = PseudoRandom.randDouble();
		   if (rnd < delta_) {
			   type = 1;
		   } else {
			   type = 2;
		   }
		   SolutionSet sols = it.next();
		   if (sols.size() == 0) {
			   System.out.println("SolsSize2 = " + sols.size());
			   System.exit(0);
		   }
		   double minFitness = 1.0e+30;
		   int minFitnessID = -1;
		   if (sols.size() == 0) {
			   System.out.println("size = 0!");
			   System.exit(0);
		   }
		   Solution sol1 = sols.getCentroidVector();
		   for (int j = 0; j < sols.size(); j++) {
			   Solution sol2 = sols.get(j);
			   if (type == 1) {
				   //double fitness = sol2.getDistanceToIdealPoint();
				   double fitness = sol2.getSumValue();
				   //double fitness = computeDistance(sol2);
				   //double fitness = computeASFFitness(sol1,sol2);
				   //double fitness = weightSumValue(sol1,sol2);
				   //double fitness = computeChebyshev(sol1,sol2);
				   //double fitness = computePBIFitness(sol1,sol2);
				   if (minFitness > fitness) {
					   minFitness = fitness;
					   minFitnessID = j;
				   }
			   } else {
				   double fitness = sol2.getSumValue();
				   //double fitness = weightSumValue(k,sol2);
				   //double fitness = computePBIFitness(sol1,sol2);
				   //double fitness = computeAngle(sol1,sol2);
				   //double fitness = computeChebyshev(sol1,sol2);
				   //double fitness = computeASFFitness(sol1,sol2);
				   //System.out.println(fitness);
				   if (minFitness > fitness) {
					   minFitness = fitness;
					   minFitnessID = j;
				   }
			   }
		   }//for
		   population_.add(sols.get(minFitnessID));
		   it.remove();
	   }//while
	   if (list.size() != 0) {
		   System.out.println("ListSize2 = " + list.size());
		   System.exit(0);
	   }
   }
}
