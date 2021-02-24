package etmo.metaheuristics.matmy2.libs;

import etmo.core.Algorithm;
import etmo.core.ProblemSet;
import etmo.core.SolutionSet;
import etmo.util.JMException;

public abstract class MaTAlgorithm extends Algorithm {

    protected SolutionSet solutionSet_;

    /**
     * Constructor
     *
     * @param problemSet The problem to be solved
     */
    public MaTAlgorithm(ProblemSet problemSet) {
        super(problemSet);
    }

    // TODO
    public MaTAlgorithm(ProblemSet problemSet, SolutionSet solutionSet){
        super(problemSet);
    }

    public abstract void initState() throws JMException, ClassNotFoundException;
    public abstract boolean step() throws JMException, ClassNotFoundException;
}
