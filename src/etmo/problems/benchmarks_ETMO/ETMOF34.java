package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
//import etmo.problems.base.*;
import etmo.problems.base.dynamicBase.DMZDT;
import etmo.problems.base.staticBase.IO;

public class ETMOF34 {
	
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
		DMZDT prob = new DMZDT(50, 1, 0, 1, fc, sc);
		
		prob.setPvType("dynamic"); //The type of the position-related variables
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHtType("static"); //The type of time-related function used to control the change of the PF
		prob.setHType("static_convex");  //The type of the shape function
		prob.setGType("sphere");   //The type of the landscape function
		prob.setF1Type("linear");  //The type of processing the position-related variables for f1 in ZDT-based formulation
			
		((Problem)prob).setName("ETMOF34_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2(int fc, int sc) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		DMZDT prob = new DMZDT(51, 2, 0, 1, fc, sc);
		
		//double[][] matrix = IO.readMatrixFromFile("MData/CEC2021/benchmark_34/matrix_1");
		//prob.setRotationMatrix(matrix);
		
		prob.setPvType("dynamic"); //The type of the position-related variables
		prob.setGtType("dynamic_t1"); //The type of time-related function used to control the change of the PS
		prob.setHtType("static"); //The type of time-related function used to control the change of the PF
		prob.setHType("static_convex");  //The type of the shape function
		prob.setGType("rastrigin");   //The type of the landscape function
		prob.setF1Type("linear");  //The type of processing the position-related variables for f1 in ZDT-based formulation
			
		((Problem)prob).setName("ETMOF34_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}

}
