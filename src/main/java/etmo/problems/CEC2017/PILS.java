package etmo.problems.CEC2017;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2017.base.*;


public class PILS {
	public static ProblemSet getProblem() throws IOException {
		ProblemSet ps1 = getT1();
		ProblemSet ps2 = getT2();
		ProblemSet problemSet = new ProblemSet(2);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		return problemSet;

	}

	public static ProblemSet getT1() {
		ProblemSet problemSet = new ProblemSet(1);

		MMDTLZ prob = new MMDTLZ(2, 50, 1, -50,50);
		prob.setGType("griewank");

		((Problem)prob).setName("PILS1");
		
		problemSet.add(prob);
		return problemSet;
	}

	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);

		MMDTLZ prob = new MMDTLZ(2, 50, 1, -100,100);
		prob.setGType("ackley");

		double[] shiftValues= IO.readShiftValuesFromFile("resources/SVData/S_PILS_2.txt");
		prob.setShiftValues(shiftValues);
		
		((Problem)prob).setName("PILS2");
		
		problemSet.add(prob);
		return problemSet;
	}
	
}
