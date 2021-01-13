package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
//import etmo.problems.base.*;
import etmo.problems.base.dynamicBase.DMDTLZ;
import etmo.problems.base.dynamicBase.DMZDT;

public class ETMOF38 {
	
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
		DMDTLZ prob = new DMDTLZ(3, 50, -1, 1, fc, sc);
		
		prob.setGtType("dynamic_t3"); //The type of time-related function used to control the change of the PS
		//prob.setHtType("dynamic_t1"); //The type of time-related function used to control the change of the PF
		//prob.setHType("inverse");  //The type of the shape function
		prob.setGType("sphere_linkage4");  //The type of the landscape function
		//prob.setPvType("dynamic");
			
		((Problem)prob).setName("ETMOF38_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMDTLZ prob = new DMDTLZ(3, 50, 0, 1, fc, sc);

		prob.setGtType("dynamic_t2"); //The type of time-related function used to control the change of the PS
		//prob.setHtType("dynamic_t1"); //The type of time-related function used to control the change of the PF
		prob.setHType("inverse");  //The type of the shape function
		prob.setGType("sphere_linkage3");  //The type of the landscape function
		prob.setPvType("dynamic");
			
		((Problem)prob).setName("ETMOF38_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}

}
