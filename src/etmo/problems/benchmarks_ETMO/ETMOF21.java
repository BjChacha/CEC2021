package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.base.staticBase.IO;
import etmo.problems.base.staticBase.MMDTLZ;
import etmo.problems.base.staticBase.MMLMOP;
import etmo.problems.base.staticBase.MMLZ;
import etmo.problems.base.staticBase.MMZDT;

public class ETMOF21 {

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
		
		MMLMOP prob = new MMLMOP(3, 512, -10,10);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("sphere"); //Shape Function
		prob.setGType("DF6"); //Landscape Function
		prob.setLinkageType("nonlinear");

		((Problem)prob).setName("ETOMF21_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLMOP prob = new MMLMOP(3, 512, -10,10);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("sphere"); //Shape Function
		prob.setGType("DF7"); //Landscape Function
		prob.setLinkageType("nonlinear");
		((Problem)prob).setName("ETOMF21_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT3() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLMOP prob = new MMLMOP(3, 512, -10,10);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("sphere"); //Shape Function
		prob.setGType("DF8"); //Landscape Function
		prob.setLinkageType("nonlinear");
		((Problem)prob).setName("ETOMF21_3");
		
		problemSet.add(prob);
		
		return problemSet;
	}
}
