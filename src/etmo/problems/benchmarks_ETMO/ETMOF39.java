package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
//import etmo.problems.base.*;
import etmo.problems.base.dynamicBase.DMLZ;
import etmo.problems.base.dynamicBase.DMZDT;
import etmo.problems.base.staticBase.IO;

public class ETMOF39 {
	
	/*
	 * fc indicates the frequency of change
	 * sc indicates the severity of change
	*/
	
	public static ProblemSet getProblem(int fc, int sc) throws IOException {
		ProblemSet ps1 = getT1(fc, sc);
		ProblemSet ps2 = getT2(fc, sc);
		ProblemSet ps3 = getT3(fc, sc);
		ProblemSet problemSet = new ProblemSet(3);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		problemSet.add(ps3.get(0));
		return problemSet;

	}
	
	public static ProblemSet getT1(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMLZ prob = new DMLZ(50,-1, 1, fc, sc);
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHt1Type("dynamic_t1"); //The type of time-related function used to control the change of the PF
		prob.setHt2Type("static"); //The type of time-related function used to control the change of the PF
		prob.setGType("sphere");   //The type of the landscape function
			
		((Problem)prob).setName("ETMOF39_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMLZ prob = new DMLZ(50,-1, 1, fc, sc);
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHt1Type("dynamic_t1"); //The type of time-related function used to control the change of the PF
		prob.setHt2Type("static"); //The type of time-related function used to control the change of the PF
		prob.setGType("rosenbrock");   //The type of the landscape function	
		((Problem)prob).setName("ETMOF39_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT3(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMLZ prob = new DMLZ(50,-1, 1, fc, sc);
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHt1Type("dynamic_t1"); //The type of time-related function used to control the change of the PF
		prob.setHt2Type("static"); //The type of time-related function used to control the change of the PF
		prob.setGType("griewank");   //The type of the landscape function	
		((Problem)prob).setName("ETMOF39_3");
		
		problemSet.add(prob);
		
		return problemSet;
	}

}
