package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.IO;
import etmo.problems.CEC2021.base.staticBase.MMIDTLZ;

import java.io.IOException;


public class ETMOF29 {
	
	public static ProblemSet getProblem() throws IOException {
		
		int taskNumber=30;
		
		ProblemSet problemSet = new ProblemSet(taskNumber);
		
		for(int i=0;i<taskNumber;i++)
			problemSet.add(getT(i).get(0));
		
		return problemSet;

	}
	
	
	public static ProblemSet getT(int taskID) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		MMIDTLZ prob;
		switch(taskID%3){
		case 0:
			prob = new MMIDTLZ(3, 51,1,  -50,50);
			prob.setGType("F1");
			break;
		case 1:
			prob = new MMIDTLZ(3, 51,1,  -100,100);
			prob.setGType("F9");
			break;		
		case 2:
			prob = new MMIDTLZ(3, 51,1,  -0.5,0.5);
			prob.setGType("F6");
			break;				
		default:
			prob = new MMIDTLZ(3, 51,1,  -50,50);
			prob.setGType("F1");
		} 
		
		switch(taskID%2){
		case 0:
			prob.setHType("inverted_lineoid");
			break;
		case 1:
			prob.setHType("inverted_sphere");
			break;		
		default:
			prob.setHType("inverted_lineoid");
		}
				
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_29/matrix_"+(taskID+1));
		
		double shiftValues[] = IO.readShiftValuesFromFile("resources/MData/CEC2021/benchmark_29/bias_"+(taskID+1));
		
		prob.setRotationMatrix(matrix);
		prob.setShiftValues(shiftValues);		
		
		((Problem)prob).setName("ETMOF29_"+(taskID+1));
		
		problemSet.add(prob);
		
		return problemSet;
	}
		
}
