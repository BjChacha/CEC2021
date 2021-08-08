package etmo.problems.CEC2017;


import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2017.base.*;


public class CIHS {
	
	public static ProblemSet getProblem() {
		ProblemSet ps1 = getT1();
		ProblemSet ps2 = getT2();
		ProblemSet problemSet = new ProblemSet(2);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		return problemSet;

	}
	
	
	public static ProblemSet getT1() {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMDTLZ prob = new MMDTLZ(2, 50, 1, -100,100);
		prob.setGType("sphere");
		
		((Problem)prob).setName("CIHS1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	
	public static ProblemSet getT2() {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMZDT prob = new MMZDT(50, 1,  -100,100);
		prob.setGType("mean");
		prob.setHType("concave");

		((Problem)prob).setName("CIHS2");
		
		problemSet.add(prob);
		return problemSet;
	}
}
