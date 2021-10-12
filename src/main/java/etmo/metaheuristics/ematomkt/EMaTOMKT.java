package etmo.metaheuristics.ematomkt;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import etmo.core.MtoAlgorithm;
import etmo.core.ProblemSet;
import etmo.core.Solution;
import etmo.core.Operator;
import etmo.core.SolutionSet;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PORanking;
import etmo.util.PseudoRandom;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.LocationComparator;
import etmo.util.math.Matrix;
import etmo.util.math.Probability;
import etmo.util.math.Vector;
import etmo.util.sorting.SortingIdx;

import smile.clustering.KMeans;
import smile.clustering.PartitionClustering;
import etmo.qualityIndicator.QualityIndicator;

public class EMaTOMKT extends MtoAlgorithm{
    private SolutionSet[] population;

    private int populationSize;
    private int taskNum;
    private int[] objStart;
    private int[] objEnd;

    private double initTransP;
    private double[] amp;
    private double[][] tObjectives;
    private double[][] differences;

    private int generations;
    private int evaluations;
    private int maxEvaluations;

    Operator crossover;
    Operator mutation;

    // DEBUG
    String[] pf;
    List<QualityIndicator> indicators;

    public EMaTOMKT(ProblemSet problemSet) {
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws ClassNotFoundException, JMException {
        initState();
        while (evaluations < maxEvaluations) {
            iterate();
        }
        return population;
    }

    private void initState() throws ClassNotFoundException, JMException {
        // initial parameters
        evaluations = 0;
        generations = 0;
        taskNum = problemSet_.size();
        maxEvaluations = (Integer) this.getInputParameter("maxEvaluations");
        populationSize = (Integer) this.getInputParameter("populationSize");

        crossover = operators_.get("crossover");
        mutation = operators_.get("mutation");

        initTransP = 0.9;
        objStart = new int[taskNum];
        objEnd = new int[taskNum];
        amp = new double[taskNum];
        tObjectives = new double[taskNum][];
        differences = new double[taskNum][];

        Arrays.fill(amp, initTransP);

        for (int k = 0; k < taskNum; k++) {
            objStart[k] = problemSet_.get(k).getStartObjPos();
            objEnd[k] = problemSet_.get(k).getEndObjPos();

            tObjectives[k] = new double[maxEvaluations / populationSize / taskNum];
            differences[k] = new double[taskNum];
            Arrays.fill(tObjectives[k], 0);
            Arrays.fill(differences[k], Double.MAX_VALUE);
        }

        // initial population
        population = new SolutionSet[taskNum];
        for (int k = 0; k < taskNum; k++) {
            population[k] = new SolutionSet(populationSize);
            for (int i = 0; i < populationSize; i++) {
                Solution solution = new Solution(problemSet_);
                solution.setSkillFactor(k);
                problemSet_.get(k).evaluate(solution);
                evaluations ++;
                population[k].add(solution);
            }
            // Non dominated sorting
            NDSort(population[k], k);
        }

        // DEBUG
        pf = new String[taskNum];
        indicators = new ArrayList<>(taskNum);
        for (int k = 0; k < taskNum; k++) {
            pf[k] = "resources/PF/StaticPF/" + problemSet_.get(k).getHType() + "_" + problemSet_.get(k).getNumberOfObjectives() + "D.pf";
            indicators.add(new QualityIndicator(problemSet_.get(k), pf[k]));
        }

        calIGD();
    }

    private void iterate() throws JMException {
        adaptiveMatingProbability();
        calculateMMD();
        // TODO: adaptive with CEC2021
        double[][][][] models = LEKT(10, 5);
        updatePopulation(models);
        // DEBUG: pring IGD
        calIGD();
    }

    private void NDSort(SolutionSet pop, int taskID) {
        for (int i = 0; i < pop.size(); i++)
            pop.get(i).setLocation(Integer.MAX_VALUE);

        boolean selec[] = new boolean[problemSet_.getTotalNumberOfObjs()];
        for (int i = objStart[taskID]; i <= objEnd[taskID]; i++) {
            selec[i] = true;
        }
        Distance distance = new Distance();
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
        pop.sort(new LocationComparator());
    }

    private double getMinDistance(SolutionSet pop) {
        int n = pop.size();
        int sf = pop.get(0).getSkillFactor();
        int d = objEnd[sf] - objStart[sf] + 1;
        double minDistance = Double.MAX_VALUE;
        for (int i = 0; i < d; i++) {
            double distance = 0;
            for (int j = 0; j < n; j ++){
                distance += Math.pow(pop.get(j).getObjective(i + objStart[sf]), 2);
            }
            minDistance = Math.min(Math.sqrt(distance), minDistance);
        }
        return minDistance;
    }

    private void adaptiveMatingProbability() {
        for (int k = 0; k < taskNum; k++) {
            tObjectives[k][generations] = getMinDistance(population[k]);
            if (generations >= 2) {
                if (tObjectives[k][generations] == tObjectives[k][generations-1]) {
                    amp[k] = initTransP;
                }
                else {
                    double d1 = Math.abs(tObjectives[k][generations] - tObjectives[k][generations-1]);
                    double d2 = Math.abs(tObjectives[k][generations-1] - tObjectives[k][generations-2]);
                    amp[k] = d1 / (d1 + d2);
                    if (Double.isNaN(amp[k])) {
                        amp[k] = initTransP;
                    }
                }
            }
        }
    }

    private void calculateMMD() throws JMException {
        double sigma = 1.0;
        for (int i = 0; i < taskNum; i++) {
            final int ii = i;
            Arrays.parallelSetAll(differences[ii], jj -> {
                double dist = 0;
                if (jj > ii) {
                    try {
                        dist = Utils.MMD(population[ii].getMat(), population[jj].getMat(), sigma);
                    } catch (JMException e) {
                        e.printStackTrace();
                    }
                } else {
                    dist = differences[jj][ii];
                }
                return dist;
            });

            // for (int j = i + 1; j < taskNum; j++) {
            //     differences[i][j] = differences[j][i] = Utils.MMD(population[i].getMat(), population[j].getMat(), sigma);
            // }
        }
    }

    private double[][][][] LEKT(int clusterNum, int selectedTaskNum) throws JMException {
        // models[taskNum][clusterNum][means][std]
        double[][][][] models = new double[taskNum][clusterNum][2][];
        
        for (int k = 0; k < taskNum; k++) {
            int[] idx = SortingIdx.sort(differences[k], false);
            SolutionSet selectedPopulations = new SolutionSet(population[k]);
            for (int i = 0; i < selectedTaskNum; i++) {
                selectedPopulations = selectedPopulations.union(population[idx[i]]);
            }
            
            // long startTime = System.currentTimeMillis();
            // double[][][] clusters = Utils.kmeans(selectedPopulations, clusterNum, populationSize);
            double[][] mat = selectedPopulations.getMat();
            var kmeans = PartitionClustering.run(1, () -> KMeans.fit(mat, clusterNum));

            // System.out.println("kmeans: " + (System.currentTimeMillis() - startTime) + "ms.");
            
            ArrayList<ArrayList<double[]>> clusters = new ArrayList<>();
            for (int i = 0; i < clusterNum; i++) {
                clusters.add(new ArrayList<double[]>());
            }
            for (int i  = 0; i < selectedPopulations.size(); i++) {
                double[] variables = selectedPopulations.get(i).getDecisionVariablesInDouble();
                int clusterID = kmeans.predict(variables);
                clusters.get(clusterID).add(variables);
            }
            for (int i = 0; i < clusterNum; i++) {
                if (clusters.get(i).size() == 0) {
                    // System.out.println("Warning: Empty cluster detected!");
                    clusters.get(i).add(kmeans.centroids[i]);
                }
                models[k][i][0] = Matrix.getMeanOfMat(clusters.get(i));
                models[k][i][1] = Matrix.getStdOfMat(clusters.get(i));
            }

            for (int i = 0; i < population[k].size(); i++) {
                double[] variables = population[k].get(i).getDecisionVariablesInDouble();
                int clusterID = kmeans.predict(variables);
                population[k].get(i).setFlag(clusterID);
            }
        }
        return models;
    }

    private void updatePopulation(double[][][][] models) throws JMException {
        for (int k = 0; k < taskNum; k++) {
            SolutionSet offspring = new SolutionSet(populationSize);
            for (int i = 0; i < population[k].size() / 2; i++) {
                int j = i;
                while (j == i) 
                    j = PseudoRandom.randInt(0, populationSize - 1);

                if (PseudoRandom.randDouble() < amp[k]) {
                    Solution[] parents = new Solution[2];
                    parents[0] = new Solution(population[k].get(i));
                    parents[1] = new Solution(population[k].get(j));
                    Solution[] children = (Solution[]) crossover.execute(parents);
                    mutation.execute(children[0]);
                    mutation.execute(children[1]);
                    children[0].setSkillFactor(k);
                    children[1].setSkillFactor(k);
                    problemSet_.get(k).evaluate(children[0]);
                    problemSet_.get(k).evaluate(children[1]);
                    offspring.add(children[0]);
                    offspring.add(children[1]);
                    evaluations += 2;
                }
                else {
                    int clusterID1 = population[k].get(i).getFlag();
                    int clusterID2 = population[k].get(j).getFlag();
                    double[] variables1 = Probability.sampleByNorm(models[k][clusterID1][0], models[k][clusterID1][1]);
                    double[] variables2 = Probability.sampleByNorm(models[k][clusterID2][0], models[k][clusterID2][1]);
                    Vector.vecClip_(variables1, 0.0, 1.0);
                    Vector.vecClip_(variables2, 0.0, 1.0);
                    Solution child1 = new Solution(population[k].get(i));
                    Solution child2 = new Solution(population[k].get(j));
                    child1.setSkillFactor(k);
                    child2.setSkillFactor(k);
                    child1.setDecisionVariables(variables1);
                    child2.setDecisionVariables(variables2);
                    child1.resetObjective();
                    child2.resetObjective();
                    problemSet_.get(k).evaluate(child1);
                    problemSet_.get(k).evaluate(child2);
                    offspring.add(child1);
                    offspring.add(child2);
                    evaluations += 2;
                }
            }

            SolutionSet union = population[k].union(offspring);
            NDSort(union, k);
            for (int i = 0; i < populationSize; i++) {
                population[k].replace(i, union.get(i));
            }
        }
        generations ++;
    }

    private void calIGD() {
        SolutionSet[] resPopulation = new SolutionSet[taskNum];
		for (int k = 0; k < taskNum; k++) {
			resPopulation[k] = new SolutionSet();
			for (int i = 0; i < population[k].size(); i++) {
				Solution sol = population[k].get(i);
				int start = problemSet_.get(k).getStartObjPos();
				int end = problemSet_.get(k).getEndObjPos();
				Solution newSolution = new Solution(end - start + 1);
				for (int kk = start; kk <= end; kk++)
					newSolution.setObjective(kk - start, sol.getObjective(kk));
				resPopulation[k].add(newSolution);
			}
		}

        double[] igd = new double[taskNum];
        for (int k = 0; k < taskNum; k++) {
            igd[k] = indicators.get(k).getIGD(resPopulation[k]);
        }

        // System.out.println("Generation " + generations + ": " + Arrays.toString(igd));
    }

}
