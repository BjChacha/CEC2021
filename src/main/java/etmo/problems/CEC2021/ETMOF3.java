package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.IO;
import etmo.problems.CEC2021.base.staticBase.MMLZ;

import java.io.IOException;

public class ETMOF3 {

	public static ProblemSet getProblem() throws IOException {
		ProblemSet ps1 = getT1();
		ProblemSet ps2 = getT2();
		ProblemSet problemSet = new ProblemSet(2);

		problemSet.add(ps1.get(0));
		problemSet.add(ps2.get(0));
		return problemSet;

	}
	
	public static ProblemSet getT(int tsk) throws IOException {
		if(tsk == 0){
			return getT1();
		}else{
			return getT2();
		}
	}
	
	public static ProblemSet getT1() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(2, 50, -10,10);
		prob.setGenType("addition");//Formulation model
		prob.setHType("convex"); //Shape Function
		prob.setGType("LF4_5"); //Landscape Function
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_3/matrix_1");
		prob.setRotationMatrix(matrix);	
		((Problem)prob).setName("ETOMF3_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(3, 51, -10,10);
		prob.setGenType("addition");//Formulation model
		prob.setHType("sphere"); //Shape Function
		prob.setGType("LF6"); //Landscape Function
		
		double[][] matrix = IO.readMatrixFromFile("resources/MData/CEC2021/benchmark_3/matrix_2");
		prob.setRotationMatrix(matrix);	
		((Problem)prob).setName("ETOMF3_2");
		
		problemSet.add(prob);
		return problemSet;
	}
}
