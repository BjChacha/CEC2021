package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.dynamicBase.DMLZ;

import java.io.IOException;

public class ETMOF40 {
	
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
		prob.setHt1Type("static"); //The type of time-related function used to control the change of the PF
		prob.setHt2Type("dynamic_t2"); //The type of time-related function used to control the change of the PF
		prob.setGType("rosenbrock");   //The type of the landscape function
			
		((Problem)prob).setName("ETMOF40_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMLZ prob = new DMLZ(50,-1, 1, fc, sc);
		
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHt1Type("static"); //The type of time-related function used to control the change of the PF
		prob.setHt2Type("dynamic_t2"); //The type of time-related function used to control the change of the PF
		prob.setGType("ackley");   //The type of the landscape function
			
		((Problem)prob).setName("ETMOF40_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT3(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMLZ prob = new DMLZ(50,-1, 1, fc, sc);
		
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHt1Type("static"); //The type of time-related function used to control the change of the PF
		prob.setHt2Type("dynamic_t2"); //The type of time-related function used to control the change of the PF
		prob.setGType("rastrigin");   //The type of the landscape function
			
		((Problem)prob).setName("ETMOF40_3");
		
		problemSet.add(prob);
		
		return problemSet;
	}

}
