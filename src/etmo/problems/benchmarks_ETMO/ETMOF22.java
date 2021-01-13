package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.base.staticBase.IO;
import etmo.problems.base.staticBase.MMDTLZ;
import etmo.problems.base.staticBase.MMLMOP;
import etmo.problems.base.staticBase.MMLZ;
import etmo.problems.base.staticBase.MMZDT;

public class ETMOF22 {

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
		
		MMLMOP prob = new MMLMOP(3, 256, -10,10);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("inverted_lineoid"); //Shape Function
		prob.setGType("DF9"); //Landscape Function
		prob.setLinkageType("linear");

		((Problem)prob).setName("ETOMF22_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLMOP prob = new MMLMOP(3, 512, -10,10);
		prob.setGenType("addition");//Formulation model
		prob.setHType("inverted_lineoid"); //Shape Function
		prob.setGType("DF9"); //Landscape Function
		prob.setLinkageType("nonlinear");
		((Problem)prob).setName("ETOMF22_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT3() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLMOP prob = new MMLMOP(3, 1024, -10,10);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("inverted_lineoid"); //Shape Function
		prob.setGType("DF9"); //Landscape Function
		prob.setLinkageType("nonlinear");
		((Problem)prob).setName("ETOMF22_3");
		
		problemSet.add(prob);
		
		return problemSet;
	}
}
