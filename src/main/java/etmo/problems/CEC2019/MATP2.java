package etmo.problems.CEC2019;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2019.base.IO;
import etmo.problems.CEC2019.base.MMZDT;

import java.io.IOException;


public class MATP2 {
	
	public static ProblemSet getProblem() throws IOException {
		
		int taskNumber=50;
		
		ProblemSet problemSet = new ProblemSet(taskNumber);
		
		for(int i=0;i<taskNumber;i++)
			problemSet.add(getT(i).get(0));
		
		return problemSet;

	}
	
	public static ProblemSet getT(int taskID) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMZDT prob = new MMZDT(50, 1, -100,100);
		prob.setGType("mean");
		prob.setHType("concave");
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2019/benchmark_2/matrix_"+(taskID+1));
		
		double shiftValues[] = IO.readShiftValuesFromFile("resources/MData/CEC2019/benchmark_2/bias_"+(taskID+1));
		
		prob.setRotationMatrix(matrix);
		prob.setShiftValues(shiftValues);		
		
		((Problem)prob).setName("MATP2-"+(taskID+1));
		
		problemSet.add(prob);
		
		return problemSet;
	}
		
}
