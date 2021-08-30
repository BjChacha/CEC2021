package etmo.problems.WCCI2020.SO;

import java.io.IOException;

import etmo.core.ProblemSet;

public class MATP1 {
    public static ProblemSet getProblem() throws IOException {
        int taskNumber = 50;
        ProblemSet problemSet = new ProblemSet(taskNumber);

        for (int i = 0; i < taskNumber; i++) {
            problemSet.add(getT(i).get(0));
        }

        return problemSet;
    }

    public static ProblemSet getT(int taskID) throws IOException {
        ProblemSet problemSet = new ProblemSet(1);

        return problemSet;
    }
}
