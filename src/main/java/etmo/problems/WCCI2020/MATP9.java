package etmo.problems.WCCI2020;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.WCCI2020.base.IO;
import etmo.problems.WCCI2020.base.MMZDT;
import java.io.IOException;


public class MATP9 {
	
	public static ProblemSet getProblem() throws IOException {
		
		int taskNumber=50;
		
		ProblemSet problemSet = new ProblemSet(taskNumber);
		
		for(int i=0;i<taskNumber;i++)
			problemSet.add(getT(i).get(0));
		
		return problemSet;

	}
	
	
	public static ProblemSet getT(int taskID) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		MMZDT prob;
		
		switch(taskID%5){
		case 0:
			prob = new MMZDT(50, 1,  -50,50);
			prob.setGType("rosenbrock");
			break;
		case 1:
			prob = new MMZDT(50, 1,  -50,50);
			prob.setGType("ackley");
			break;		
		case 2:
			prob = new MMZDT(50, 1,  -50,50);
			prob.setGType("rastrigin");
			break;
		case 3:
			prob = new MMZDT(50, 1,  -100,100);
			prob.setGType("griewank");
			break;
		case 4:
			prob = new MMZDT(50, 1, -0.5, 0.5);
			prob.setGType("weierstrass");
			break;					
		default:
			prob = new MMZDT(50, 1, -100,100);
		}   
		prob.setHType("convex");
				
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/WCCI2020/benchmark_9/matrix_"+(taskID+1));
		
		double shiftValues[] = IO.readShiftValuesFromFile("resources/MData/WCCI2020/benchmark_9/bias_"+(taskID+1));
		
		prob.setRotationMatrix(matrix);
		prob.setShiftValues(shiftValues);		
		
		((Problem)prob).setName("MATP9-"+(taskID+1));
		
		problemSet.add(prob);
		
		return problemSet;
	}
		
}
