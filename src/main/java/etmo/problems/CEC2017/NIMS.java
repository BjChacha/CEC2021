package etmo.problems.CEC2017;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2017.base.*;


public class NIMS {

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
		
		MMDTLZ prob = new MMDTLZ(3, 20, 1, -20,20);
		prob.setGType("rosenbrock");

		((Problem)prob).setName("NIMS1");
		
		problemSet.add(prob);
		return problemSet;
	}

	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMZDT prob = new MMZDT(20, 2,  -20,20);
		prob.setGType("sphere");
		prob.setHType("concave");
		
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/M_NIMS_2.txt");
		prob.setRotationMatrix(matrix);
		
		((Problem)prob).setName("NIMS2");
		
		
		problemSet.add(prob);
		return problemSet;
	}
}
