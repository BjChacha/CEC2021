package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.base.staticBase.IO;
import etmo.problems.base.staticBase.MMDTLZ;
import etmo.problems.base.staticBase.MMLMOP;
import etmo.problems.base.staticBase.MMLZ;
import etmo.problems.base.staticBase.MMZDT;

public class ETMOF17 {

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
		
		MMLMOP prob = new MMLMOP(3, 256, -10,10);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("lineoid"); //Shape Function
		prob.setGType("DF1"); //Landscape Function
		prob.setLinkageType("linear");

		((Problem)prob).setName("ETOMF17_1");
		
		problemSet.add(prob);
		
		return problemSet;
	}
	
	public static ProblemSet getT2() throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		
		MMLMOP prob = new MMLMOP(3, 256, -10,10);
		prob.setGenType("multiplication");//Formulation model
		prob.setHType("inverted_lineoid"); //Shape Function
		prob.setGType("DF1"); //Landscape Function
		prob.setLinkageType("nonlinear");

		((Problem)prob).setName("ETOMF17_2");
		
		problemSet.add(prob);
		
		return problemSet;
	}
}
