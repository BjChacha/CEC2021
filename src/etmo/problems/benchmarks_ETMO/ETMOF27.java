package etmo.problems.benchmarks_ETMO;

import java.io.IOException;

import etmo.core.Problem;
import etmo.core.ProblemSet;
//import etmo.problems.base.*;
import etmo.problems.base.staticBase.IO;
import etmo.problems.base.staticBase.MMLMOP;
import etmo.problems.base.staticBase.MMLZ;
import etmo.problems.base.staticBase.MMZDT;


public class ETMOF27 {
	
	public static ProblemSet getProblem() throws IOException {
		
		int taskNumber=10;
		
		ProblemSet problemSet = new ProblemSet(taskNumber);
		
		for(int i=0;i<taskNumber;i++)
			problemSet.add(getT(i).get(0));
		
		return problemSet;

	}
	
	
	public static ProblemSet getT(int taskID) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		MMLMOP prob;
		switch(taskID%5){
		case 0:
			prob = new MMLMOP(3, 50, 3, -10,10);
			prob.setGType("DF14");
			break;
		case 1:
			prob = new MMLMOP(3, 50, 3, -10,10);
			prob.setGType("DF15");
			break;		
		case 2:
			prob = new MMLMOP(3, 50, 3, -10,10);
			prob.setGType("DF16");
			break;
		case 3:
			prob = new MMLMOP(3, 50, 3, -10,10);
			prob.setGType("DF17");
			break;
		case 4:
			prob = new MMLMOP(3, 50, 3, -10,10);
			prob.setGType("DF18");
			break;					
		default:
			prob = new MMLMOP(3, 50, 3, -10,10);
			prob.setGType("DF1");
		} 
		
		if((double)taskID/5 <= 1){
			prob.setGenType("multiplication");
			prob.setLinkageType("linear");
		}else{
			prob.setGenType("addition");
			prob.setLinkageType("nonlinear");
		}
		prob.setHType("lineoid");
		
				
		((Problem)prob).setName("ETMOF27_"+(taskID+1));
		
		problemSet.add(prob);
		
		return problemSet;
	}
		
}
