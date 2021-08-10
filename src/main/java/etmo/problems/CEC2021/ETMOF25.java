package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.IO;
import etmo.problems.CEC2021.base.staticBase.MMLZ;

import java.io.IOException;


public class ETMOF25 {
	
	public static ProblemSet getProblem() throws IOException {
		
		int taskNumber=5;
		
		ProblemSet problemSet = new ProblemSet(taskNumber);
		
		for(int i=0;i<taskNumber;i++)
			problemSet.add(getT(i).get(0));
		
		return problemSet;

	}
	
	
	public static ProblemSet getT(int taskID) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		MMLZ prob;
		
		switch(taskID){
		case 0:
			prob = new MMLZ(2, 50,  -10,10);
			prob.setGenType("addition");
			prob.setGType("LF1");
			prob.setHType("convex");
			break;
		case 1:
			prob = new MMLZ(2, 50,  -5,5);
			prob.setGenType("addition");
			prob.setGType("LF2");
			prob.setHType("convex");
			break;		
		case 2:
			prob = new MMLZ(2, 50,  -1,1);
			prob.setGenType("addition");
			prob.setGType("LF3");
			prob.setHType("convex");
			break;
		case 3:
			prob = new MMLZ(2, 50, -5,5);
			prob.setGenType("addition");
			prob.setGType("LF4_5");
			prob.setHType("convex");
			break;
		case 4:
			prob = new MMLZ(2, 50,  -10,10);
			prob.setGenType("addition");
			prob.setGType("LF7");
			prob.setHType("convex");
			break;					
		default:
			prob = new MMLZ(2, 50,  -10,10);
			prob.setGenType("addition");
			prob.setGType("LF7");
			prob.setHType("convex");
		}   
				
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_25/matrix_"+(taskID+1));
		
		prob.setRotationMatrix(matrix);		
		
		((Problem)prob).setName("ETMOF25_"+(taskID+1));
		
		problemSet.add(prob);
		
		return problemSet;
	}
		
}
