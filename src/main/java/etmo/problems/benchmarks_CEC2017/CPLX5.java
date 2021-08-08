package etmo.problems.benchmarks_CEC2017;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.LZ09.*;


public class CPLX5 {
	
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
		
		LZ09_F3 prob = new LZ09_F3(2,30);
		
		((Problem)prob).setName("CPLX5-1");
		prob.setHType("LZ09_F3");		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		LZ09_F6 prob = new LZ09_F6(3,10);

		((Problem)prob).setName("CPLX5-2");
		prob.setHType("LZ09_F6");		
		problemSet.add(prob);
		
		return problemSet;
	}
}
