package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.MMLMOP;

import java.io.IOException;

public class ETMOF18 {

	public static ProblemSet getProblem() throws IOException {
		ProblemSet ps1 = getT1();
		ProblemSet ps2 = getT2();
		ProblemSet problemSet = new ProblemSet(2);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		return problemSet;

	}
	
	public static ProblemSet getT(int tsk) throws IOException {
		if(tsk == 0){
			return getT1();
		}else{
			return getT2();
		}
	}
	
	public static ProblemSet getT1() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLMOP prob = new MMLMOP(2, 512, -10,10);
		prob.setGenType("addition");//Formulation model
		prob.setHType("circle"); //Shape Function
		prob.setGType("DF1"); //Landscape Function
		prob.setLinkageType("linear");

		((Problem)prob).setName("ETOMF18_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLMOP prob = new MMLMOP(2, 512, -10,10);
		prob.setGenType("addition");//Formulation model
		prob.setHType("convex"); //Shape Function
		prob.setGType("DF2"); //Landscape Function
		prob.setLinkageType("linear");

		((Problem)prob).setName("ETOMF18_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}
}
