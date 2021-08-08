package etmo.problems.CEC2017;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2017.base.*;


public class NILS {

	public static ProblemSet getProblem() throws IOException {
		ProblemSet ps1 = getT1();
		ProblemSet ps2 = getT2();
		ProblemSet problemSet = new ProblemSet(2);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		return problemSet;

	}

	public static ProblemSet getT1() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		MMDTLZ prob = new MMDTLZ(3, 25, 1, -50,50);
		prob.setGType("griewank");

		
		double shiftValues[] = IO.readShiftValuesFromFile("resources/SVData/S_NILS_1.txt");
		prob.setShiftValues(shiftValues);
				
		((Problem)prob).setName("NILS1");
		
		problemSet.add(prob);
		return problemSet;
	}

	public static ProblemSet getT2() {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMZDT prob = new MMZDT(50, 2,  -100,100);
		prob.setGType("ackley");
		prob.setHType("concave");
		
		((Problem)prob).setName("NILS2");

		problemSet.add(prob);
		return problemSet;
	}
}
