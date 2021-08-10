package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.IO;
import etmo.problems.CEC2021.base.staticBase.MMDTLZ;

import java.io.IOException;

public class ETMOF14 {
	
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
		
		MMDTLZ prob = new MMDTLZ(5, 53,1,-1,1);
		prob.setGenType("addition");
		prob.setHType("disconnect");
		prob.setGType("HF6");
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_14/matrix_1");
		prob.setRotationMatrix(matrix);	
		
		((Problem)prob).setName("ETMOF14_1");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMDTLZ prob = new MMDTLZ(8, 56,1, -1,1);
		prob.setGenType("addition");
		prob.setHType("disconnect");
		prob.setGType("HF7");
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_14/matrix_2");
		prob.setRotationMatrix(matrix);	
		
		
		((Problem)prob).setName("ETMOF14_2");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	public static ProblemSet getT3() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMDTLZ prob = new MMDTLZ(10, 58,1,-1,1);
		prob.setGenType("multiplication");
		prob.setHType("disconnect");
		prob.setGType("HF8");
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_14/matrix_3");
		prob.setRotationMatrix(matrix);	
		
		((Problem)prob).setName("ETMOF14_3");
		
		problemSet.add(prob);
		return problemSet;
	}
}
