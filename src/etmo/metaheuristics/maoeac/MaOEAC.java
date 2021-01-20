package etmo.metaheuristics.maoeac;

import etmo.core.*;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.Ranking;
//import jmetal.util.PseudoRandom;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public class MaOEAC extends Algorithm {
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

    private double[] zideal_; //ideal point
    private double[] znadir_;//Nadir point

    /**
     * Constructor
     *
     * @param problemSet The problem to be solved
     */
    public MaOEAC(ProblemSet problemSet) {
        super(problemSet);
        zideal_ = new double[problemSet.get(0).getNumberOfObjectives()];
        znadir_ = new double[problemSet.get(0).getNumberOfObjectives()];
    }

    @Override
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
         * Enter the main loop，into the process of evolution
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
        }

        Ranking ranking = new Ranking(population_);
        return ranking.getSubfront(0);
    }

    private void environmentalSelection() {

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
            if(solutionSets[k].size() != 2*populationSize_/problemSet_.get(0).getNumberOfObjectives()){
                System.out.println("solutionSets["+k+"] = "+ solutionSets[k].size());
                System.exit(0);
            }
//			按照支配关系层取前面的
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

    private void bestSolutionSelection(List<SolutionSet> list, int k) {
        double minClustering2Axis = 1.0e+30;
        int minClustering2AxisID = -1;
        for(int i=0;i<list.size();i++){
            SolutionSet sols = list.get(i);
            if(sols.size() == 0){
                System.out.println("SolsSize1 = "+sols.size());
                System.exit(0);
            }

            double angle1 = Math.acos(Math.abs(sols.getCentroidVector().getNormalizedObjective(k)/sols.getCentroidVector().getDistanceToIdealPoint()));


            if(angle1 < minClustering2Axis){
                minClustering2Axis = angle1;
                minClustering2AxisID = i;
            }//if
        }//for
        double minSolution2Axis = 1.0e+30;
        int minSolution2AxisID = -1;
        for(int j=0;j<list.get(minClustering2AxisID).size();j++){
            Solution sol = list.get(minClustering2AxisID).get(j);
            double ang = Math.acos(list.get(minClustering2AxisID).get(j).getNormalizedObjective(k)/list.get(minClustering2AxisID).get(j).getDistanceToIdealPoint());
            if(ang < minSolution2Axis){
                minSolution2Axis = ang;
                minSolution2AxisID = j;
            }
        }//for
        population_.add(list.get(minClustering2AxisID).get(minSolution2AxisID));
        list.remove(minClustering2AxisID);

        double min2CenterLine = 1.0e+30;
        int min2CenterLineId = -1;
        for(int i=0;i<list.size();i++){
            SolutionSet sols = list.get(i);
			/*if(sols.size() != 0){
				System.out.println("SolsSize1 = "+sols.size());
				//System.exit(0);
			}*/
            double sumValue = 0.0;
            for(int j=0;j<problemSet_.get(0).getNumberOfObjectives();j++){
                sumValue += sols.getCentroidVector().getNormalizedObjective(j)*1.0;
            }
            //System.out.println("value = "+sumValue/(sols.getCentroidVector().getDistanceToIdealPoint()));
            //norm2 = Math.sqrt(norm2);
            double angle2 = Math.acos(sumValue/(sols.getCentroidVector().getDistanceToIdealPoint()*Math.sqrt(problemSet_.get(0).getNumberOfObjectives())));
            //System.out.println(angle2);
            if(angle2 < min2CenterLine){
                min2CenterLine = angle2;
                min2CenterLineId = i;
            }
        }
        //System.out.println(min2CenterLineId);
        double minS2CenterLine = 1.0e+30;
        int minId = -1;
        for(int i=0;i<list.get(min2CenterLineId).size();i++){
            Solution sol = list.get(min2CenterLineId).get(i);
            double sumValue = 0.0;
            for(int j=0;j<problemSet_.get(0).getNumberOfObjectives();j++){
                sumValue += sol.getNormalizedObjective(j);
            }
            double ang = Math.acos(Math.abs(sumValue/(sol.getDistanceToIdealPoint()*Math.sqrt(problemSet_.get(0).getNumberOfObjectives()))));
            if(ang < minS2CenterLine){
                minS2CenterLine = ang;
                minId = i;
            }
        }
        if(PseudoRandom.randDouble() < 0.5){
            population_.add(list.get(min2CenterLineId).get(minId));
            list.remove(min2CenterLineId);
        }



		/*if(list.size() != populationSize_/(problem_.getNumberOfObjectives()) - 2){
			System.out.println("ListSize3 = "+list.size());
			System.exit(0);
		}*/
        double delta_ = 1.0;
        Iterator<SolutionSet> it = list.iterator();
        while(it.hasNext()){
            int type = -1;
            double rnd = PseudoRandom.randDouble();
            if (rnd < delta_)
            {
                type = 1;
            } else {
                type = 2;
            }
            SolutionSet sols = it.next();
            if(sols.size() == 0){
                System.out.println("SolsSize2 = "+sols.size());
                System.exit(0);
            }
            double minFitness = 1.0e+30;
            int minFitnessID = -1;
            if(sols.size() == 0){
                System.out.println("size = 0!");
                System.exit(0);
            }
            Solution sol1 = sols.getCentroidVector();
            for(int j=0;j<sols.size();j++){
                Solution sol2 = sols.get(j);
                if(type == 1){
                    //double fitness = sol2.getDistanceToIdealPoint();
                    double fitness = sol2.getSumValue();
                    //double fitness = computeDistance(sol2);
                    //double fitness = computeASFFitness(sol1,sol2);
                    //double fitness = weightSumValue(sol1,sol2);
                    //double fitness = computeChebyshev(sol1,sol2);
                    //double fitness = computePBIFitness(sol1,sol2);
                    if(minFitness > fitness){
                        minFitness = fitness;
                        minFitnessID = j;
                    }
                }else{
                    double fitness = sol2.getSumValue();
                    //double fitness = weightSumValue(k,sol2);
                    //double fitness = computePBIFitness(sol1,sol2);
                    //double fitness = computeAngle(sol1,sol2);
                    //double fitness = computeChebyshev(sol1,sol2);
                    //double fitness = computeASFFitness(sol1,sol2);
                    //System.out.println(fitness);
                    if(minFitness > fitness){
                        minFitness = fitness;
                        minFitnessID = j;
                    }
                }
            }//for
            population_.add(sols.get(minFitnessID));
            it.remove();
        }//while
        if(list.size() != 0){
            System.out.println("ListSize2 = "+list.size());
            System.exit(0);
        }
    }

    private SolutionSet getStSolutionSet(SolutionSet ss, int size) {
        SolutionSet sets = new SolutionSet();
        Ranking ranking = new Ranking(ss);

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

    private void computeDistanceToIdealPoint(SolutionSet solutionSet) {
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

    private void normalizationObjective(SolutionSet solutionSet) {
        for(int i=0; i<solutionSet.size(); i++){
            Solution sol = solutionSet.get(i);

            for(int j=0; j<problemSet_.get(0).getNumberOfObjectives(); j++){
                double val = 0.0;
                val = (sol.getObjective(j) - zideal_[j])/(znadir_[j]-zideal_[j]);
                //val = (sol.getObjective(j) - zideal_[j]);
                sol.setNormalizedObjective(j, val);
            }
        }


    }

    private void estimateNadirPoint(SolutionSet solutionSet) {
        for(int i=0; i<problemSet_.get(0).getNumberOfObjectives();i++){
            znadir_[i] = -1.0e+30;
            for(int j=0; j<solutionSet.size();j++){
                if(solutionSet.get(j).getObjective(i) > znadir_[i]){
                    znadir_[i] = solutionSet.get(j).getObjective(i);
                }
            }

        }
    }

    private void estimateIdealPoint(SolutionSet solutionSet) {
        for(int i=0; i< problemSet_.get(0).getNumberOfObjectives();i++){
            zideal_[i] = 1.0e+30;
            for(int j=0; j<solutionSet.size();j++){
                if(solutionSet.get(j).getObjective(i) < zideal_[i]){
                    zideal_[i] = solutionSet.get(j).getObjective(i);
                }
            }

        }
    }

    private void generateOffspringPopulation() throws JMException {
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
                    parents[0] = (Solution) selection_.execute(population_);
                    parents[1] = (Solution) selection_.execute(population_);
                }else{
                    parents[0] = (Solution) selection_.execute(population_);
                    parents[1] = (Solution) selection_.execute(population_);
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

    private void initPopulation() throws ClassNotFoundException, JMException {
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

    }

    private void baiscSetting() {
        generations_ = 0;
        maxGenerations_ = ((Integer) this.getInputParameter("maxGenerations")).intValue();
        populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
        mutation_  = operators_.get("mutation");
        crossover_ = operators_.get("crossover");
        selection_ = operators_.get("selection");

        learning_ = operators_.get("learning");
    }


}
