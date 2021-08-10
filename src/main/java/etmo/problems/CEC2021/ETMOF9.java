package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.IO;
import etmo.problems.CEC2021.base.staticBase.MMIDTLZ;

import java.io.IOException;

public class ETMOF9 {
	
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
		
		MMIDTLZ prob = new MMIDTLZ(5, 25, 1, -10,10);
		prob.setGType("F13");
		prob.setHType("inverted_lineoid");
		
		((Problem)prob).setName("ETMOF9_1");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		
		MMIDTLZ prob = new MMIDTLZ(5, 53, 1, -100,100);
		prob.setGType("F13");
		prob.setHType("inverted_sphere");
		
		double[] shiftValues = IO.readShiftValuesFromFile("resources/MData/CEC2021/benchmark_9/bias_2");
		prob.setShiftValues(shiftValues);	
		
		
		((Problem)prob).setName("ETMOF9_2");
		
		problemSet.add(prob);
		return problemSet;
	}

}
