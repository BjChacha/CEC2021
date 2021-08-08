package etmo.metaheuristics.experiments;

import etmo.core.*;
import etmo.util.JMException;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;

public class findPF_main {
    public static void main(String[] args) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, JMException, IOException {
        ProblemSet problemSet;
        MtoAlgorithm algorithm;

        int problemStart = 25;
        int problemEnd = 32;
        for (int pCase = problemStart; pCase <= problemEnd; pCase++) {
            System.out.println("Processing Task " + pCase + " ...");
            problemSet = (ProblemSet) Class
                    .forName("etmo.problems.benchmarks_ETMO.ETMOF" + pCase)
                    .getMethod("getProblem")
                    .invoke(null, null);

            int taskNum = problemSet.size();
            String[] pf = new String[taskNum];
            for (int k = 0; k < pf.length; k++) {
                pf[k] = "PF/StaticPF/" + problemSet.get(k).getHType() + "_" + problemSet.get(k).getNumberOfObjectives() + "D.pf";
            }
            algorithm = new findPF(problemSet);
            algorithm.execute();
        }
    }
}
