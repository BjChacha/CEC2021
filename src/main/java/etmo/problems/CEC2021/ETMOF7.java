package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.IO;
import etmo.problems.CEC2021.base.staticBase.MMLZ;

import java.io.IOException;

public class ETMOF7 {
	
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
		
		MMLZ prob = new MMLZ(2, 50, -50,50);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("convex"); //Shape Function
		prob.setGType("LF2"); //Landscape Function
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_7/matrix_1");
		prob.setRotationMatrix(matrix);	
		
		((Problem)prob).setName("ETMOF7_1");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(2, 50, -50,50);
		prob.setGenType("addition");//Formulation model
		prob.setHType("concave"); //Shape Function
		prob.setGType("LF3"); //Landscape Function
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_7/matrix_2");
		prob.setRotationMatrix(matrix);	
		
		
		((Problem)prob).setName("ETMOF7_2");
		
		problemSet.add(prob);
		return problemSet;
	}
	
	public static ProblemSet getT3() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(2, 50, -50,50);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("lineoid"); //Shape Function
		prob.setGType("LF10"); //Landscape Function
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_7/matrix_2");
		prob.setRotationMatrix(matrix);	
		
		
		((Problem)prob).setName("ETMOF7_3");
		
		problemSet.add(prob);
		return problemSet;
	}

}
