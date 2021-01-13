package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
//import etmo.problems.base.*;
import etmo.problems.base.dynamicBase.DMLZ;
import etmo.problems.base.dynamicBase.DMZDT;
import etmo.problems.base.staticBase.IO;

public class ETMOF35 {
	
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
		DMLZ prob = new DMLZ(50,-1, 1, fc, sc);
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHt1Type("dynamic_t1"); //The type of time-related function used to control the change of the PF
		prob.setHt2Type("static"); //The type of time-related function used to control the change of the PF
		prob.setGType("sphere");   //The type of the landscape function
			
		((Problem)prob).setName("ETMOF35_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMLZ prob = new DMLZ(50,-1, 1, fc, sc);
		
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHt1Type("static"); //The type of time-related function used to control the change of the PF
		prob.setHt2Type("dynamic_t2"); //The type of time-related function used to control the change of the PF
		prob.setGType("rastrigin");   //The type of the landscape function
			
		((Problem)prob).setName("ETMOF35_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}

}
