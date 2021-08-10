package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.IO;
import etmo.problems.CEC2021.base.staticBase.MMDTLZ;
import etmo.problems.CEC2021.base.staticBase.MMLZ;

import java.io.IOException;

public class ETMOF13 {
	
	public static ProblemSet getProblem() throws IOException {
		ProblemSet ps1 = getT1();
		ProblemSet ps2 = getT2();
		ProblemSet ps3 = getT3();
		ProblemSet problemSet = new ProblemSet(3);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		problemSet.add(ps3.get(0));
		return problemSet;

	}
	
	
	public static ProblemSet getT1() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMDTLZ prob = new MMDTLZ(5, 53,1,-10,10);
		prob.setGenType("multiplication");
		prob.setHType("degenerate");
		prob.setGType("F9");
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_13/matrix_1");
		prob.setRotationMatrix(matrix);	
		
		((Problem)prob).setName("ETMOF13_1");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(8, 56, -10,10);
		prob.setGenType("addition");
		prob.setHType("sphere");
		prob.setGType("LF2");
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_13/matrix_2");
		prob.setRotationMatrix(matrix);	
		
		
		((Problem)prob).setName("ETMOF13_2");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	public static ProblemSet getT3() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(10, 58, -10,10);
		prob.setGenType("addition");
		prob.setHType("irconcave");
		prob.setGType("LF3");
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_13/matrix_3");
		prob.setRotationMatrix(matrix);	
		
		
		((Problem)prob).setName("ETMOF13_3");
		
		problemSet.add(prob);
		return problemSet;
	}
}
