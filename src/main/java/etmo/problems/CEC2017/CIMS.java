package etmo.problems.CEC2017;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2017.base.*;


public class CIMS {
	
	public static ProblemSet getProblem() throws IOException {
		ProblemSet ps1 = getT1();
		ProblemSet ps2 = getT2();
		ProblemSet problemSet = new ProblemSet(2);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		return problemSet;

	}
	
	
	public static ProblemSet getT1() {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMZDT prob = new MMZDT(10,  1, -5,5);
		prob.setGType("rosenbrock");
		prob.setHType("concave");
		
		
		((Problem)prob).setName("CIMS1");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		
		MMDTLZ prob = new MMDTLZ(2, 10, 1, -5,5);
		prob.setGType("mean");
		
		
		double[] shiftValues = IO.readShiftValuesFromFile("resources/SVData/S_CIMS_2.txt");
		prob.setShiftValues(shiftValues);
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/M_CIMS_2.txt");
		prob.setRotationMatrix(matrix);
		
		
		((Problem)prob).setName("CIMS2");
		
		problemSet.add(prob);
		return problemSet;
	}
}
