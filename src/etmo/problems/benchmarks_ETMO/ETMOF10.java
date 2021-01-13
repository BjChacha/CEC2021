package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.base.staticBase.IO;
import etmo.problems.base.staticBase.MMDTLZ;
import etmo.problems.base.staticBase.MMIDTLZ;
import etmo.problems.base.staticBase.MMLZ;

public class ETMOF10 {
	
	public static ProblemSet getProblem() throws IOException {
		ProblemSet ps1 = getT1();
		ProblemSet ps2 = getT2();
		ProblemSet problemSet = new ProblemSet(2);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		return problemSet;
	}
	
	
	public static ProblemSet getT1() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(8, 56, -20,20);
		prob.setGType("F14");
		prob.setHType("irconcave");
		prob.setGenType("multiplication");
		
		double[][] matrix = IO.readMatrixFromFile("MData/CEC2021/benchmark_10/matrix_1");
		prob.setRotationMatrix(matrix);	
		
		((Problem)prob).setName("ETMOF10_1");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(8, 56, -10,10);
		prob.setGType("F15");
		prob.setHType("irconcave");
		prob.setGenType("addition");
		
		double[][] matrix = IO.readMatrixFromFile("MData/CEC2021/benchmark_10/matrix_2");
		prob.setRotationMatrix(matrix);	
		
		((Problem)prob).setName("ETMOF10_2");
		
		problemSet.add(prob);
		return problemSet;
	}

}
