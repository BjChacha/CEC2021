package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.base.staticBase.IO;
import etmo.problems.base.staticBase.MMDTLZ;
import etmo.problems.base.staticBase.MMLZ;
import etmo.problems.base.staticBase.MMZDT;

public class ETMOF1 {

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
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("convex"); //Shape Function
		prob.setGType("LF1"); //Landscape Function

		double[][] matrix = IO.readMatrixFromFile("MData/CEC2021/benchmark_1/matrix_1");
		prob.setRotationMatrix(matrix);	
		((Problem)prob).setName("ETOMF1_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLZ prob = new MMLZ(2, 50, -10,10);
		prob.setGenType("addition");//Formulation model
		prob.setHType("concave");//Shape Function
		prob.setGType("LF1");//Landscape Function
		
		double[][] matrix = IO.readMatrixFromFile("MData/CEC2021/benchmark_1/matrix_2");
		prob.setRotationMatrix(matrix);
		((Problem)prob).setName("ETOMF1_2");
		
		problemSet.add(prob);//add中修改start和end值

		return problemSet;
	}
}
