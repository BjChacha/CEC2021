package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
//import etmo.problems.base.*;
import etmo.problems.base.dynamicBase.DMZDT;

public class ETMOF33 {
	
	/*
	 * fc indicates the frequency of change
	 * sc indicates the severity of change
	*/
	
	public static ProblemSet getProblem(int fc, int sc) throws IOException {
		ProblemSet ps1 = getT1(fc, sc);
		ProblemSet ps2 = getT2(fc, sc);
		ProblemSet problemSet = new ProblemSet(2);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		return problemSet;

	}
	
	public static ProblemSet getT1(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMZDT prob = new DMZDT(256, 1, -1, 1, fc, sc);
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHtType("dynamic_t1"); //The type of time-related function used to control the change of the PF
		prob.setHType("dynamic");  //The type of the shape function
		prob.setGType("sphere");   //The type of the landscape function
		prob.setF1Type("linear");  //The type of processing the position-related variables for f1 in ZDT-based formulation
			
		((Problem)prob).setName("ETMOF33_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMZDT prob = new DMZDT(256, 1, -1, 2, fc, sc);
		
		prob.setGtType("dynamic_t2"); //The type of time-related function used to control the change of the PS
		prob.setHtType("dynamic_t2"); //The type of time-related function used to control the change of the PF
		prob.setHType("dynamic");  //The type of the shape function
		prob.setGType("sphere_linkage");  //The type of the landscape function
		prob.setF1Type("linear");  //The type of processing the position-related variables for f1 in ZDT-based formulation
			
		((Problem)prob).setName("ETMOF33_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}

}
