package etmo.metaheuristics.matbml;

import etmo.core.*;
import etmo.metaheuristics.matmy2.libs.MOEAD;
import etmo.metaheuristics.matmy2.libs.MaOEAC;
import etmo.metaheuristics.matmy2.libs.MaTAlgorithm;
import etmo.metaheuristics.matmy2.libs.NSGAII;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PORanking;
import etmo.util.PseudoRandom;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.LocationComparator;
import etmo.util.logging.LogPopulation;

import javax.swing.plaf.basic.BasicSliderUI;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class MaTBML extends MtoAlgorithm {
    MaTAlgorithm[] optimizers_;
    SolutionSet[] populations_;
    Operator crossover_;
    Operator selection_;

    String algoName_;
    int populationSize_;
    int maxEvaluations_;
    int evaluations_;
    int taskNum_;
    int objNum_;
    int minGroupNum_;
    int k1_;
    int k2_;
    int implicitTransferNum_;
    // implicit transfer probability
    double P_;
    
    int[] objStart_;
    int[] objEnd_;
    boolean[] canTransfer_;
    boolean[] isFinished_;
    double[] ideals_;
    double[] nadirs_;

    double[][] scores;

    List<Integer> leaderTaskIdx_;

    List<List<Integer>> groups;

    public MaTBML(ProblemSet problemSet){
        super(problemSet);
    }

    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        // Initialization
        initState();
        initPopulations();
        initOptimizers();

        // Algorithm execution
        while (evaluations_ < maxEvaluations_){
             solelyConverge(k1_);
            transferConverge(k2_);
        }

        return populations_;
    }

    private void solelyConverge(int k1) throws JMException, ClassNotFoundException {
        double[] oldIdeal = ideals_.clone();
        Converge(k1);
        updateIdealPoint();
        updateNadirPoint();

//        clusters = Utils.AGNES(populations_, minGroupNum_);
        groups = Grouping();

//        // DEBUG
//        for (int c = 0; c < groups.size(); c++)
//            System.out.println(evaluations_ + ": " + groups.get(c));

        leaderTaskIdx_.clear();
        Arrays.fill(canTransfer_, true);
        for (int g = 0; g < groups.size(); g++) {
            double maxImprovement = 0;
            int leader = -1;
            for (int k: groups.get(g)) {
                double improvement = 0;
                for (int j = objStart_[k]; j <= objEnd_[k]; j++) {
                    improvement += (oldIdeal[j] - ideals_[j]);
                }
                if (improvement > maxImprovement) {
                    leader = k;
                    maxImprovement = improvement;
                }
            }
            leaderTaskIdx_.add(leader);
        }
    }

    private void transferConverge(int k2) throws JMException, ClassNotFoundException {
        for (int g = 0; g < groups.size(); g++) {
            if (leaderTaskIdx_.get(g) >= 0) {
                // DEBUG
                int cbnw = 0;  // complete better and not worse
                int cbw = 0;   // complete better and worse
                int pbnw = 0;  // partly better and not worse
                int pbw = 0;   // partly better and worse
                int nbnw = 0;  // not better and not worse
                int nbw = 0;   // not better and worse

                int[] runTimes = new int[taskNum_];
                Arrays.fill(runTimes, 0);
                runTimes[leaderTaskIdx_.get(g)] = k2;
                Converge(runTimes);

                for (int k : groups.get(g)) {
                    if (k == leaderTaskIdx_.get(g) || isFinished_[k])
                        continue;

                    double[] tmpIdeal = ideals_.clone();
                    double[] tmpNadir = nadirs_.clone();

                    SolutionSet tmpSet = new SolutionSet(populationSize_);

                    // evaluate the leader population on other task
                    for (int i = 0; i < populations_[leaderTaskIdx_.get(g)].size(); i++) {
                        Solution tmp = new Solution(populations_[leaderTaskIdx_.get(g)].get(i));
                        tmp.setSkillFactor(k);
                        problemSet_.get(k).evaluate(tmp);
                        evaluations_++;
                        for (int j = objStart_[k]; j <= objEnd_[k]; j++) {
                            tmpIdeal[j] = Math.min(tmpIdeal[j], tmp.getObjective(j));
                            tmpNadir[j] = Math.max(tmpNadir[j], tmp.getObjective(j));
                        }
                        tmpSet.add(tmp);
                    }

                    // apply leader population or not
                    boolean completelyBetter = true;
                    boolean partlyBetter = false;

                    boolean isWorse = false;
                    for (int j = objStart_[k]; j <= objEnd_[k]; j++) {
                        if (tmpIdeal[j] < ideals_[j]) {
                            partlyBetter = true;
                        }else{
                            completelyBetter = false;
                        }
                        if (tmpNadir[j] > nadirs_[j])
                            isWorse = true;
                    }

                    if (completelyBetter) {
                        // apply the whole population
                        for (int i = 0; i < tmpSet.size(); i++) {
                            populations_[k].replace(i, tmpSet.get(i));
                        }
                        scores[k][leaderTaskIdx_.get(g)] += 2;
                        cbnw ++;
                    }
                    else if (completelyBetter && isWorse){
                        cbw ++;
                    }
                    else if (partlyBetter) {
                        // Union and selection
                        SolutionSet union = populations_[k].union(tmpSet);
                        rankSolutionOnTask(union, k);
                        union.sort(new LocationComparator());
                        for (int i = 0; i < populations_[k].size(); i++) {
                            populations_[k].replace(i, union.get(i));
                        }
                        scores[k][leaderTaskIdx_.get(g)] += 1;
                        if (!isWorse){
                            pbnw ++;
                        }else{
                            pbw ++;
                        }
                    }
                    else if (!isWorse){
                        nbnw ++;
//                        scores[k][leaderTaskIdx_.get(g)] = 0;
                        if (groups.get(g).size() < 2) continue;
                        if (PseudoRandom.randDouble() < P_) {
                            // select transfer source task base on score
//                            double[] scoreList = new double[groups.get(g).size()];
//                            for (int i = 0; i < scoreList.length; i++) {
//                                scoreList[i] = scores[k][groups.get(g).get(i)];
//                            }
//                            int transferIdx = k;
//                            while (transferIdx == k)
//                                transferIdx = Utils.roulette(scoreList);
                            int transferIdx = leaderTaskIdx_.get(g);

                            // get individuals used for transfer
                            SolutionSet transferPopulation = new SolutionSet(implicitTransferNum_);
                            for (int i = 0; i < implicitTransferNum_; i++) {
                                Solution newOne = new Solution((Solution) selection_.execute(populations_[transferIdx]));
                                transferPopulation.add(newOne);
                            }

                            SolutionSet offspring = new SolutionSet(populationSize_);
                            for (int i = 0; i < populations_[transferIdx].size(); i++){
                                Solution[] parents = new Solution[2];
                                parents[0] = populations_[k].get(i);
                                parents[1] = populations_[transferIdx].get(i);
                                Solution[] children = (Solution[]) crossover_.execute(parents);
                                Solution child = new Solution(children[PseudoRandom.randInt(0, children.length - 1)]);
                                child.setSkillFactor(k);
                                problemSet_.get(k).evaluate(child);
                                evaluations_++;
                                offspring.add(child);
                            }
                            SolutionSet union = populations_[k].union(offspring);
                            rankSolutionOnTask(union, k);
                            union.sort(new LocationComparator());
                            for (int i = 0; i < populations_[k].size(); i++) {
                                populations_[k].replace(i, union.get(i));
                            }
                        }
                    }
                    else{
                        scores[k][leaderTaskIdx_.get(g)] += -1;
                        nbw ++;
                    }


                    if (((evaluations_ / 100) * 100) % (20 * taskNum_ * populationSize_) == 0) {
                        LogPopulation.LogPopulation("MaTBML", populations_, problemSet_, evaluations_);
                    }
                }

//                // DEBUG
//                System.out.println("Leader Task: " + leaderTaskIdx_.get(g));
//                System.out.println("\tCompletely better and not worse: " + cbnw +
//                        "\tCompletely better and worse: " + cbw +
//                        "\tPartly better and not worse: " + pbnw +
//                        "\tPartly better and worse: " + pbw +
//                        "\tNot better and not worse: " + nbnw +
//                        "\tNot better and worse: " + nbw
//                        );
            }
            else{

            }
        }
        updateIdealPoint();
        updateNadirPoint();
    }

    private void rankSolutionOnTask(SolutionSet pop, int taskId) {
        for (int i = 0; i < pop.size(); i++)
            pop.get(i).setLocation(Integer.MAX_VALUE);

        boolean selec[] = new boolean[problemSet_.getTotalNumberOfObjs()];

        for (int i = 0; i < selec.length; i++) {
            if (i < objStart_[taskId] || i > objEnd_[taskId])
                selec[i] = false;
            else
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
    }

    private List<List<Integer>> Grouping(){
        // 1. Initialize groups
        List<List<Integer>> groups = new ArrayList<>();
        for (int i = 0; i < minGroupNum_; i++)
            groups.add(new ArrayList<Integer>());

        // 2. iterate begin with a random task
        int[] perm = Utils.permutation(taskNum_, taskNum_);
        for (int k: perm){
            // the score list of groups for task k
            List<Double> kthScores = new ArrayList<>();
            for (List<Integer> group: groups){
                double tmpScore = 1;
                // calculate the total score of this group
                for (int i = 0; i < group.size(); i++){
                    tmpScore += scores[k][group.get(i)];
                }
                tmpScore /= Math.max(1, group.size() - taskNum_ / minGroupNum_ + 2);
                kthScores.add(tmpScore);
            }
            double[] scoreList = kthScores.stream().mapToDouble(Double::doubleValue).toArray();
            // if the score of all groups is negative, then create a new group
            boolean needNewGroup = true;
            for (double score: scoreList){
                if (score >= 0) {
                    needNewGroup = false;
                    break;
                }
            }
            int groupIdx;
            if (needNewGroup){
                groupIdx = groups.size();
                groups.add(new ArrayList<>());
            } else {
                groupIdx = Utils.roulette(scoreList);
            }
            groups.get(groupIdx).add(k);
        }

        return groups;
    }

    private void updateIdealPoint() {
        for (int k = 0; k < taskNum_; k++){
            for (int j = objStart_[k]; j <= objEnd_[k]; j++){
                for (int i = 0; i < populations_[k].size(); i++){
                    ideals_[j] = Math.min(ideals_[j], populations_[k].get(i).getObjective(j));
                }
            }
        }
    }

    private void updateNadirPoint() {
        for (int k = 0; k < taskNum_; k++){
            for (int j = objStart_[k]; j <= objEnd_[k]; j++){
                for (int i = 0; i < populations_[k].size(); i++){
                    nadirs_[j] = Math.max(nadirs_[j], populations_[k].get(i).getObjective(j));
                }
            }
        }
    }

    private void Converge(int times) throws JMException, ClassNotFoundException {
        int[] tmpTimes = new int[taskNum_];
        Arrays.fill(tmpTimes, times);
        Converge(tmpTimes);
    }

    private void Converge(int[] times) throws JMException, ClassNotFoundException {
        for (int k = 0; k < taskNum_; k++) {
            for (int t = 0; t < times[k]; t++) {
                if (!isFinished_[k]) {
                    isFinished_[k] = !optimizers_[k].step();
                    evaluations_ += populationSize_;

                    if (((evaluations_ / 100) * 100) % (20 * taskNum_ * populationSize_) == 0) {
                        LogPopulation.LogPopulation("MaTBML", populations_, problemSet_, evaluations_);
                    }
                }
            }
        }
    }

    private void initState() {
        evaluations_ = 0;
        populationSize_ = (Integer) getInputParameter("populationSize");
        maxEvaluations_ = (Integer) getInputParameter("maxEvaluations");
        k1_ = (Integer) getInputParameter("k1");
        k2_ = (Integer) getInputParameter("k2");
        implicitTransferNum_ = (Integer) getInputParameter("implicitTransferNum");
        P_ = (double) getInputParameter("P_");
        algoName_ = (String) getInputParameter("algoName");

        crossover_ = operators_.get("crossover");
        selection_ = operators_.get("selection");

        taskNum_ = problemSet_.size();
        objNum_ = problemSet_.getTotalNumberOfObjs();
        minGroupNum_ = (int) Math.sqrt(taskNum_);
        isFinished_ = new boolean[taskNum_];
        canTransfer_ = new boolean[minGroupNum_];
        ideals_ = new double[objNum_];
        nadirs_ = new double[objNum_];
        objStart_ = new int[taskNum_];
        objEnd_ = new int[taskNum_];

        scores = new double[taskNum_][taskNum_];

        leaderTaskIdx_ = new ArrayList<>();
        
//        Arrays.fill(isFinished_, false);
        Arrays.fill(ideals_, Double.MAX_VALUE);
//        Arrays.fill(nadirs_, 0);

        for (int k = 0; k < taskNum_; k++){
            objStart_[k] = problemSet_.get(k).getStartObjPos();
            objEnd_[k] = problemSet_.get(k).getEndObjPos();

            Arrays.fill(scores[k], 0);
        }
    }

    private void initPopulations() throws JMException, ClassNotFoundException {
        populations_ = new SolutionSet[taskNum_];
        for (int k = 0; k < taskNum_; k++){
            populations_[k] = new SolutionSet(populationSize_);
            for (int i = 0; i < populationSize_; i++){
                Solution individual = new Solution(problemSet_);
                individual.setSkillFactor(k);
                problemSet_.get(k).evaluate(individual);
                evaluations_ ++;
                populations_[k].add(individual);
            }
        }
        updateIdealPoint();
        updateNadirPoint();

        LogPopulation.LogPopulation("MaTBML", populations_, problemSet_, evaluations_);
    }

    private void initOptimizers() throws JMException, ClassNotFoundException {
        optimizers_ = new MaTAlgorithm[taskNum_];
        if (algoName_.equalsIgnoreCase("MOEAD")){
            Operator crossover;
            Operator mutation;

            HashMap parameters;

            for (int k = 0; k < taskNum_; k++) {
                optimizers_[k] = new MOEAD(problemSet_, populations_[k], k);
                optimizers_[k].setInputParameter("populationSize", 100);
                optimizers_[k].setInputParameter("maxEvaluations", 1000 * 100);

                optimizers_[k].setInputParameter("dataDirectory", "D:\\_r\\EA\\ETMO\\MTO-cec2021-\\resources\\weightVectorFiles\\moead");

                optimizers_[k].setInputParameter("T", 20);
                optimizers_[k].setInputParameter("delta", 0.9);
                optimizers_[k].setInputParameter("nr", 2);

                parameters = new HashMap();
                parameters.put("CR", 1.0);
                parameters.put("F", 0.5);
                crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet_.get(k).getNumberOfVariables());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                optimizers_[k].addOperator("crossover", crossover);
                optimizers_[k].addOperator("mutation", mutation);

                optimizers_[k].initState();
            }
        }
        else if (algoName_.equalsIgnoreCase("MaOEAC")){
            Operator crossover;
            Operator mutation;
            Operator selection;
            HashMap parameters;

            for (int k = 0; k < taskNum_; k++) {
                optimizers_[k] = new MaOEAC(problemSet_, populations_[k], k);

                optimizers_[k].setInputParameter("populationSize",100);
                optimizers_[k].setInputParameter("maxGenerations",1000);

                parameters = new HashMap();
                parameters.put("probability", 1.0);
                parameters.put("distributionIndex", 30.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet_.get(k).getNumberOfVariables());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                parameters = null;
                selection = SelectionFactory.getSelectionOperator("RandomSelection", parameters);

                optimizers_[k].addOperator("crossover", crossover);
                optimizers_[k].addOperator("mutation", mutation);
                optimizers_[k].addOperator("selection", selection);

                optimizers_[k].initState();
            }
        }
        else if (algoName_.equalsIgnoreCase("NSGAII")){
            Operator crossover;
            Operator mutation;
            Operator selection;

            HashMap parameters;

            for (int k = 0; k < taskNum_; k++){
                ProblemSet pS = problemSet_.getTask(k);
                optimizers_[k] = new NSGAII(pS, populations_[k]);

                optimizers_[k].setInputParameter("populationSize", 100);
                optimizers_[k].setInputParameter("maxEvaluations", 100 * 1000);

                parameters = new HashMap();
                parameters.put("probability", 0.9);
                parameters.put("distributionIndex", 20.0);
                crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

                parameters = new HashMap();
                parameters.put("probability", 1.0 / pS.getMaxDimension());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                parameters = null;
                selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters);

                optimizers_[k].addOperator("crossover", crossover);
                optimizers_[k].addOperator("mutation", mutation);
                optimizers_[k].addOperator("selection", selection);

                optimizers_[k].initState();
            }
        }
        else{
            throw new JMException("Exception in " + algoName_ + ": Not such algorithm.");
        }
    }
}
