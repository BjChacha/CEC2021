package etmo.metaheuristics.mfeaddra;

import etmo.core.*;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PORanking;
import etmo.util.PseudoRandom;
import etmo.util.comparators.CrowdingComparator;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

public class MFEADDRA extends Algorithm {
    private int populationSize;
    private int taskNum;
    private SolutionSet population;
    
    double beta;
    int nr;
    int T;
    double rmp;
    int subSize;

    int rate;
    int neighborType;

    double[] idealPoint;
    double[][][] lambda;
    int[][][] neighborhood;

    String dataDirectory_;
    
    List<Solution> savedValues;

    int evaluations;
    int generation;
    int maxEvaluations;

    Operator crossover;
    Operator mutation;

    Distance distance = new Distance();

    public MFEADDRA(ProblemSet problemSet) {
        super(problemSet);
    }

    @Override
    public SolutionSet execute() throws JMException, ClassNotFoundException {
        initState();
        while (evaluations < maxEvaluations) {
            iteration();
        }

        return population;
    }

    protected void initState() throws JMException, ClassNotFoundException {
        populationSize = (Integer) getInputParameter("populationSize");
        maxEvaluations = (Integer) getInputParameter("maxEvaluations");
        rmp = (Double) getInputParameter("rmp");
        beta = (Double) getInputParameter("beta");
        nr = (Integer) getInputParameter("nr");
        T = (Integer) getInputParameter("T");

        taskNum = problemSet_.size();
        subSize = populationSize / taskNum;

        dataDirectory_ = this.getInputParameter("dataDirectory").toString();

        crossover = operators_.get("crossover");
        mutation = operators_.get("mutation");

        savedValues = new ArrayList<>(populationSize);

        rate = 30;

        evaluations = 0;
        generation = 0;

        initUniformWeight();
        initNeighborhood();
        initPopulation();
        initIdealPoint();
    }

    void iteration() throws JMException {
        List<Integer> order = tourSelection(10);
        for (int i = 0; i < order.size(); i++) {
            int subProblemId = order.get(i);

            int skillFactor = population.get(subProblemId).getSkillFactor();
            int subProblem = subProblemId % subSize;
            chooseNeighborType();

            Solution[] parents = parentSelection(skillFactor, subProblem);
            Solution child = (Solution) crossover.execute(new Object[] {parents[2], parents});
            mutation.execute(child);

            int id = chooseTask(skillFactor);
            child.setSkillFactor(id);
            child.resetObjective();
            problemSet_.get(id).evaluate(child);

            this.evaluations++;
            updateIdealPoint(child, id);
            updateNeighborhood(child, id, subProblem);
        }

        generation ++;

        if (rate != 0 && generation % rate == 0) {
            updateUtility();
        }
    }

    private List<Integer> tourSelection(int depth) {
        List<Integer> selected = new ArrayList<Integer>();
        List<Integer> candidate = new ArrayList<Integer>();

        for (int n = 0; n < populationSize; n++) {
            candidate.add(n);
        }

        while (selected.size() < (int) (populationSize / 5.0)) {
            int best_idd = (int) (PseudoRandom.randDouble() * candidate.size());
            int i2;
            int best_sub = candidate.get(best_idd);
            int s2;
            for (int i = 1; i < depth; i++) {
                i2 = (int) (PseudoRandom.randDouble() * candidate.size());
                s2 = candidate.get(i2);
                if (population.get(s2).getFitness() > population.get(best_idd).getFitness()) {
                    best_idd = i2;
                    best_sub = s2;
                }
            }
            selected.add(best_sub);
            candidate.remove(best_idd);
        }
        return selected;
    }

    protected void chooseNeighborType() {
        double rnd = PseudoRandom.randDouble();
        if (rnd < this.beta) {
            this.neighborType = 1;
        } else {
            this.neighborType = 2;
        }
    }

    protected Solution[] parentSelection(int skillFactor, int subProblem) {
        List<Integer> matingPool = matingSelection(skillFactor, subProblem, 2);

        Solution[] parents = new Solution[3];

        parents[0] = (population.get(matingPool.get(0)));
        parents[1] = (population.get(matingPool.get(1)));
        parents[2] = (population.get(skillFactor * subSize + subProblem));

        return parents;
    }

    protected List<Integer> matingSelection(int skillFactor, int subproblem, int numberOfSolutionsToSelect) {
        int neighbourSize;
        int selectedSolution;

        List<Integer> listOfSolutions = new ArrayList<>(numberOfSolutionsToSelect);

        neighbourSize = neighborhood[skillFactor][subproblem].length;
        while (listOfSolutions.size() < numberOfSolutionsToSelect) {
            int random;
            if (neighborType == 1) {
                random = PseudoRandom.randInt(0, neighbourSize - 1);
                selectedSolution = neighborhood[skillFactor][subproblem][random];
            } else {
                selectedSolution = subSize * skillFactor + PseudoRandom.randInt(0, subSize - 1);
            }
            boolean flag = true;

            if (skillFactor * subSize + subproblem == selectedSolution) {
                flag = false;
            }

            for (Integer individualId : listOfSolutions) {
                if (individualId == selectedSolution) {
                    flag = false;
                    break;
                }
            }

            if (flag) {
                listOfSolutions.add(selectedSolution);
            }
        }

        return listOfSolutions;
    }

    private int chooseTask(int skillFactor) {
        double rd = PseudoRandom.randDouble();

        if (rd < rmp) {
            // redefine neighborType
            neighborType = 2;

            int randTask = PseudoRandom.randInt(0, taskNum - 1);
            while (randTask == skillFactor) {
                randTask = PseudoRandom.randInt(0, taskNum - 1);
            }

            return randTask;
        } else {
            return skillFactor;
        }
    }

    protected void updateIdealPoint(Solution child, int id){
        int offset = problemSet_.get(id).getStartObjPos();
        for (int i = 0; i < problemSet_.get(id).getNumberOfObjectives(); i++){
            idealPoint[i+offset] = Math.min(idealPoint[i+offset], child.getObjective(i+offset));
        }
    }

    protected void updateNeighborhood(Solution child, int id, int subProblem) {
        int size;
        int time = 0;

        if (neighborType == 1) {
            size = neighborhood[id][subProblem].length;
        } else {
            size = subSize;
        }
        int[] perm = new int[size];

        Utils.randomPermutation(perm, size);

        for (int i = 0; i < size; i++) {
            int k;
            if (neighborType == 1) {
                k = neighborhood[id][subProblem][perm[i]];
            } else {
                k = id * subSize + perm[i];
            }
            double f1, f2;

            f1 = fitnessFunction(id, population.get(k), lambda[id][k % subSize]);
            f2 = fitnessFunction(id, child, lambda[id][k % subSize]);

            if (f2 < f1) {
                population.replace(k, new Solution(child));
                time++;
            }

            if (time >= nr) {
                return;
            }
        }
    }

    protected double fitnessFunction(int id, Solution individual, double[] lambda) {
        double maxFun = -1.0e+30;
        int offset = problemSet_.get(id).getStartObjPos();
        for (int n = 0; n < problemSet_.get(id).getNumberOfObjectives(); n++) {
            double diff = Math.abs(
                    individual.getObjective(n+offset) -
                    idealPoint[n+offset]);

            double feval;
            if (lambda[n] == 0) {
                feval = diff / 0.000001;
            } else {
                feval = diff / lambda[n];
            }
            if (feval > maxFun) {
                maxFun = feval;
            }
        }

        return maxFun;
    }

    private void updateUtility() {
        double f1, f2, uti, delta;
        for (int n = 0; n < populationSize; n++) {
            int skillFactor = population.get(n).getSkillFactor();
            int subproblem = n % subSize;

            f1 = fitnessFunction(skillFactor, population.get(n), lambda[skillFactor][subproblem]);
            f2 = fitnessFunction(skillFactor, savedValues.get(n), lambda[skillFactor][subproblem]);

            // delta = f2 - f1;
            delta = (f2 - f1) / f2;


            if (delta > 0.001) {
                population.get(n).setFitness(1.0);
            } else {
//                uti = (0.95 + (0.05 * delta / 0.001)) * utility[n];
                uti = (0.95 + (0.05 * delta / 0.001)) * ((double) population.get(n).getFitness());
                population.get(n).setFitness(Math.min(uti, 1.0));
            }
            savedValues.set(n, new Solution(population.get(n)));
        }
    }

    public void initUniformWeight() { // init lambda vectors
        lambda = new double[taskNum][][];
        for (int k = 0; k < taskNum; k++) {
            lambda[k] = new double[subSize][problemSet_.get(k).getNumberOfObjectives()];
            if (problemSet_.get(k).getNumberOfObjectives() == 2) {
                for (int n = 0; n < subSize; n++) {
                    double a = 1.0 * n / (subSize - 1);
                    lambda[k][n][0] = a;
                    lambda[k][n][1] = 1 - a;
                } // for
            } // if
            else {
                String dataFileName;
                dataFileName = "W" + problemSet_.get(k).getNumberOfObjectives() + "D_"
                        + subSize + ".dat";

                try {
                    // Open the file
                    FileInputStream fis = new FileInputStream(dataDirectory_ + "/"
                            + dataFileName);
                    InputStreamReader isr = new InputStreamReader(fis);
                    BufferedReader br = new BufferedReader(isr);

                    int i = 0;
                    int j;
                    String aux = br.readLine();
                    while (aux != null) {
                        StringTokenizer st = new StringTokenizer(aux);
                        j = 0;
                        while (st.hasMoreTokens()) {
                            double value = Double.parseDouble(st.nextToken());
                            lambda[k][i][j] = value;
                            j++;
                        }
                        aux = br.readLine();
                        i++;
                    }
                    br.close();
                } catch (Exception e) {
                    System.out
                            .println("initUniformWeight: failed when reading for file: "
                                    + dataDirectory_ + "/" + dataFileName);
                    e.printStackTrace();
                }
            } // else
        }
    } // initUniformWeight

    public void initNeighborhood() {
        neighborhood = new int[taskNum][subSize][T];
        for (int k = 0; k < taskNum; k++) {
            double[] x = new double[subSize];
            int[] idx = new int[subSize];

            for (int i = 0; i < subSize; i++) {
                for (int j = 0; j < subSize; j++) {
                    x[j] = Utils.distVector(lambda[k][i], lambda[k][j]);
                    idx[j] = subSize * k + j;
                }

                Utils.minFastSort(x, idx, subSize, T);
                System.arraycopy(idx, 0, neighborhood[k][i], 0, T);
            } // for
        }
    } // initNeighborhood

    void initIdealPoint() throws JMException, ClassNotFoundException {
        idealPoint = new double[problemSet_.getTotalNumberOfObjs()];
        Arrays.fill(idealPoint, Double.MAX_VALUE);
        for (int i = 0; i < populationSize; i++){
            int sf = population.get(i).getSkillFactor();
            int offset = problemSet_.get(sf).getStartObjPos();
            for (int j = 0; j < problemSet_.get(sf).getNumberOfObjectives(); j++){
                idealPoint[j+offset] = Math.min(idealPoint[j+offset], population.get(i).getObjective(j+offset));
            }
        }
    } // initIdealPoint


    void initPopulation() throws JMException, ClassNotFoundException {
        population = new SolutionSet(populationSize);
        for (int i = 0; i < populationSize; i++) {
            int id = i / subSize;
            Solution newSolution = new Solution(problemSet_);
            problemSet_.get(id).evaluate(newSolution);
            problemSet_.get(id).evaluateConstraints(newSolution);
            evaluations++;
            newSolution.setSkillFactor(id);
            population.add(newSolution);

            savedValues.add(new Solution(newSolution));
        } // for
        assignFitness(population);
    } // initPopulation

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
}
