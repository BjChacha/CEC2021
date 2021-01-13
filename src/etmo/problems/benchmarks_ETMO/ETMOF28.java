package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
//import etmo.problems.base.*;
import etmo.problems.base.staticBase.IO;
import etmo.problems.base.staticBase.MMDTLZ;
import etmo.problems.base.staticBase.MMLZ;
import etmo.problems.base.staticBase.MMZDT;


public class ETMOF28 {
	
	public static ProblemSet getProblem() throws IOException {
		
		int taskNumber=20;
		
		ProblemSet problemSet = new ProblemSet(taskNumber);
		
		for(int i=0;i<taskNumber;i++)
			problemSet.add(getT(i).get(0));
		
		return problemSet;

	}
	
	
	public static ProblemSet getT(int taskID) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		MMDTLZ prob;
		switch(taskID%5){
		case 0:
			prob = new MMDTLZ(3, 51,1,  -50,50);
			prob.setGType("F1");
			break;
		case 1:
			prob = new MMDTLZ(3, 51,1,  -10,10);
			prob.setGType("F5");
			break;		
		case 2:
			prob = new MMDTLZ(3, 51,1,  -20,20);
			prob.setGType("F6");
			break;
		case 3:
			prob = new MMDTLZ(3, 51,1,  -30,30);
			prob.setGType("F8");
			break;
		case 4:
			prob = new MMDTLZ(3, 51,1,  -40,40);
			prob.setGType("F9");
			break;					
		default:
			prob = new MMDTLZ(3, 51,1,  -50,50);
			prob.setGType("F1");
		} 
		
		if(taskID % 2 == 0 ){
			prob.setGenType("addition");
		}else{
			prob.setGenType("multiplication");
		}
		
		switch(taskID%3){
		case 0:
			prob.setHType("lineoid");
			break;
		case 1:
			prob.setHType("sphere");
			break;		
		case 2:
			prob.setHType("convex");
			break;
		default:
			prob.setHType("sphere");
		}
				
		double[][] matrix = IO.readMatrixFromFile("MData/CEC2021/benchmark_28/matrix_"+(taskID+1));
		
		prob.setRotationMatrix(matrix);		
		
		((Problem)prob).setName("ETMOF28_"+(taskID+1));
		
		problemSet.add(prob);
		
		return problemSet;
	}
		
}
