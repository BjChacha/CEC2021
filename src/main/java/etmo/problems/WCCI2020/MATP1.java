package etmo.problems.WCCI2020;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.WCCI2020.base.IO;
import etmo.problems.WCCI2020.base.MMDTLZ;
import java.io.IOException;


public class MATP1 {
	
	public static ProblemSet getProblem() throws IOException {
		
		int taskNumber=50;
		
		ProblemSet problemSet = new ProblemSet(taskNumber);
		
		for(int i=0;i<taskNumber;i++)
			problemSet.add(getT(i).get(0));
		
		return problemSet;

	}
	
	
	public static ProblemSet getT(int taskID) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMDTLZ prob = new MMDTLZ(2, 50, 1, -100,100);
		prob.setGType("sphere");
				
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/WCCI2020/benchmark_1/matrix_"+(taskID+1));
		
		double shiftValues[] = IO.readShiftValuesFromFile("resources/MData/WCCI2020/benchmark_1/bias_"+(taskID+1));
		
		prob.setRotationMatrix(matrix);
		prob.setShiftValues(shiftValues);		
		
		((Problem)prob).setName("MATP1-"+(taskID+1));
		
		problemSet.add(prob);
		
		return problemSet;
	}
		
}
