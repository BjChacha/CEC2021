package etmo.problems.CEC2021;

import etmo.core.Problem;
import etmo.core.ProblemSet;
import etmo.problems.CEC2021.base.staticBase.MMLMOP;

import java.io.IOException;


public class ETMOF32 {
	
	public static ProblemSet getProblem() throws IOException {
		
		int taskNumber=28;
		
		ProblemSet problemSet = new ProblemSet(taskNumber);
		
		for(int i=0;i<taskNumber;i++)
			problemSet.add(getT(i).get(0));
		
		return problemSet;

	}
	
	public static ProblemSet getT(int taskID) throws IOException {
		ProblemSet problemSet = new ProblemSet(1);
		MMLMOP prob;
		switch(taskID%14){
		case 0:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF1");
			break;
		case 1:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF19");
			break;		
		case 2:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF20");
			break;
		case 3:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF21");
			break;
		case 4:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF14");
			break;
		case 5:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF15");
			break;
		case 6:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF16");
			break;
		case 7:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF17");
			break;
		case 8:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF18");
			break;
		case 9:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF22");
			break;
		case 10:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF23");
			break;
		case 11:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF24");
			break;
		case 12:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF25");
			break;
		case 13:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF2");
			break;
		default:
			prob = new MMLMOP(3, 80, 3, -10,10);
			prob.setGType("DF1");
		} 
		
		prob.setGenType("multiplication");
		prob.setHType("sphere");
		if((double)taskID/14 < 1){
			prob.setLinkageType("linear");
		}else{
			prob.setLinkageType("nonlinear");
		}
		
		((Problem)prob).setName("ETMOF32_"+(taskID+1));
		
		problemSet.add(prob);
		
		return problemSet;
	}
		
}
