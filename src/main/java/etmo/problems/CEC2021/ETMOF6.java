package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.IO;
import etmo.problems.CEC2021.base.staticBase.MMDTLZ;
import etmo.problems.CEC2021.base.staticBase.MMZDT;

import java.io.IOException;

public class ETMOF6 {
	
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
		
		MMDTLZ prob = new MMDTLZ(2, 50, 1, -100,100);
		prob.setGType("F13");

		
		double[] shiftValues = IO.readShiftValuesFromFile("resources/MData/CEC2021/benchmark_6/bias_1");
		prob.setShiftValues(shiftValues);
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_6/matrix_1");
		prob.setRotationMatrix(matrix);	
		
		((Problem)prob).setName("ETMOF6_1");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMZDT prob = new MMZDT(50, 1, -100,100);
		prob.setHType("concave");//Shape Function
		prob.setGType("HF5");//Landscape Function
		
		double[] shiftValues = IO.readShiftValuesFromFile("resources/MData/CEC2021/benchmark_6/bias_2");
		prob.setShiftValues(shiftValues);
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_6/matrix_2");
		prob.setRotationMatrix(matrix);	
		
		
		((Problem)prob).setName("ETMOF6_2");
		
		problemSet.add(prob);
		return problemSet;
	}

}
