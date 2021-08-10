package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.MMLZ;

import java.io.IOException;

public class ETMOF11 {
	
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
		
		MMLZ prob = new MMLZ(10, 50, -20,20);
		prob.setGenType("multiplication");
		prob.setHType("convex");
		prob.setGType("LF6");
		
		((Problem)prob).setName("ETMOF11_1");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(10, 50, -10,10);
		prob.setGenType("addition");
		prob.setHType("convex");
		prob.setGType("LF7");
		
		((Problem)prob).setName("ETMOF11_2");
		
		problemSet.add(prob);
		return problemSet;
	}
}
