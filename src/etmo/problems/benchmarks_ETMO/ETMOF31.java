package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
//import etmo.problems.base.*;
import etmo.problems.base.staticBase.IO;
import etmo.problems.base.staticBase.MMLZ;
import etmo.problems.base.staticBase.MMZDT;


public class ETMOF31 {
	
	public static ProblemSet getProblem() throws IOException {
		
		int taskNumber=50;
		
		ProblemSet problemSet = new ProblemSet(taskNumber);
		
		for(int i=0;i<taskNumber;i++)
			problemSet.add(getT(i).get(0));
		
		return problemSet;

	}
	
	
	public static ProblemSet getT(int taskID) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		MMLZ prob;
		switch(taskID%6){
		case 0:
			prob = new MMLZ(2, 50,  -60,60);
			prob.setGType("LF1");
			break;
		case 1:
			prob = new MMLZ(2, 50,  -50,50);
			prob.setGType("LF2");
			break;		
		case 2:
			prob = new MMLZ(2, 50,  -40,40);
			prob.setGType("LF3");
			break;
		case 3:
			prob = new MMLZ(2, 50, -30,30);
			prob.setGType("LF8");
			break;
		case 4:
			prob = new MMLZ(2, 50,  -20,20);
			prob.setGType("LF9");
			break;
		case 5:
			prob = new MMLZ(2, 50,  -10,10);
			prob.setGType("LF7");
			break;
		default:
			prob = new MMLZ(2, 50,  -5,5);
			prob.setGType("LF1");
		} 
		
		prob.setGenType("addition");
		prob.setHType("concave");
		
				
		double[][] matrix = IO.readMatrixFromFile("MData/CEC2021/benchmark_31/matrix_"+(taskID+1));
		
		prob.setRotationMatrix(matrix);		
		
		((Problem)prob).setName("ETMOF31_"+(taskID+1));
		
		problemSet.add(prob);
		
		return problemSet;
	}
		
}
