package etmo.metaheuristics.experiments;

import etmo.core.*;
import etmo.util.*;
import etmo.util.comparators.DominanceComparator;

import java.util.Arrays;

public class findPF extends MtoAlgorithm {
    SolutionSet[] population;
    SolutionSet archive;
    int taskNum;

    double step;


    public findPF(ProblemSet problemSet) {
        super(problemSet);
    }


    @Override
    public SolutionSet[] execute() throws JMException, ClassNotFoundException {
        initState();

        Solution s = new Solution(problemSet_.getTask(0));
        int len = s.getDecisionVariables().length;
        int[] idx = new int[len];
        Arrays.fill(idx, 0);

        while (idx[len-1] <= 10){
            for (int i = 0; i < len; i++){
                s.setDecisionVariables(i, idx[i] * step);
            }

            problemSet_.get(0).evaluate(s);

            boolean added = false;
            for (int i = 0; i < archive.size(); i++){
                int flag = new DominanceComparator().compare(s, archive.get(i));
                if (flag > 0){
                    added = true;
                    break;
                }
                else if (flag < 0){
                    if (!added){
                        archive.replace(i, s);
                        added = true;
                    }else{
                        archive.remove(i);
                    }
                }
            }
            if (!added)
                archive.add(s);

            int id = 0;
            idx[id] += 1;
            while (id < len-1 && idx[id] > 10){
                idx[id] = 0;
                idx[++id] ++;
            }

            System.out.println(Arrays.toString(idx));
        }
        return population;
    }


    private void initState() {
        taskNum = problemSet_.size();
        population = new SolutionSet[taskNum];
        archive = new SolutionSet();
        step = 0.1;
    }
}
