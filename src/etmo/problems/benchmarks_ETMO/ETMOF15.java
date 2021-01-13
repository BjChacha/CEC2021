package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.base.staticBase.IO;
import etmo.problems.base.staticBase.MMDTLZ;
import etmo.problems.base.staticBase.MMIDTLZ;
import etmo.problems.base.staticBase.MMLMOP;
import etmo.problems.base.staticBase.MMLZ;

public class ETMOF15 {
	
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
		
		MMLMOP prob = new MMLMOP(10, 99, 3, -10,10);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("sphere"); //Shape Function
		prob.setGType("DF1"); //Landscape Function
		prob.setLinkageType("nonLinkage");

		((Problem)prob).setName("ETOMF15_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLMOP prob = new MMLMOP(10, 99, 3, -10,10);
		prob.setGenType("addition");//Formulation model
		prob.setHType("sphere"); //Shape Function
		prob.setGType("DF2"); //Landscape Function
		prob.setLinkageType("nonLinkage");

		((Problem)prob).setName("ETOMF15_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}

}
