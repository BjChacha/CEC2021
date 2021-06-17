package etmo.problems.base.staticBase;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.util.JMException;

public class MMLMOP extends Problem {
	String gType_;
	String genType_;
	String linkageType_;
	
	int[] nip;//Number of independent variables in each position-related variable group
    int[] nop;//Number of overlapping variables in each position-related variable group
    int nsp;//Number of shared variables in each position-related variable group
    int[] nid;//Number of independent variables in each distance-related variable group
    int[] nod;//Number of overlapping variables in each distance-related variable group
    int nsd;// Number of shared variables in each distance-related variable group
    int K;//Number of position-related variables
    int L;//Number of distance-related variables
    int[][] gp;
    int[][] gd;
    int[][][] dgd;

	public MMLMOP() {
		numberOfObjectives_ = 2;
	}

	public MMLMOP(int numberOfObjectives, int numberOfVariables, double lg, double ug) {
		numberOfObjectives_ = numberOfObjectives;
		numberOfVariables_ = numberOfVariables;
	
		gType_ = "sphere";
		hType_ = "sphere";
		genType_ = "multiplication";
		linkageType_ = "linear";
		
		nip = new int[numberOfObjectives_-1];
		nop = new int[numberOfObjectives_-1];
		int sumNIP = 0;
		for(int i=0; i<numberOfObjectives_-1; i++){
			nip[i] = 5;
			if(numberOfObjectives_ == 2){
				nop[i] = 0;
			}else{
				nop[i] = 1;
			}
			sumNIP += nip[i];
		}
		nsp = 1;
		K = nsp + sumNIP;
		L = numberOfVariables_ - K;
		
		double[] c = new double[numberOfObjectives_];
		c[0] = 3.8*0.1*(1.0-0.1);
		double sumC = c[0];
		for(int i = 1; i < numberOfObjectives_; i++){
			c[i] = 3.8*c[i-1]*(1.0-c[i-1]);
			sumC += c[i];
		}
		nid = new int[numberOfObjectives_];
		int sumNID = 0;
		for(int i = 0; i < numberOfObjectives_; i++){
			nid[i] = (int)Math.floor((c[i]/sumC)*L);
			sumNID += nid[i];
		}
		nod = new int[numberOfObjectives_];
		for(int i = 0; i < numberOfObjectives_; i++){
			if(i == 0){
				nod[i] = nid[numberOfObjectives_-1]/5;
			}else{
				nod[i] = nid[i-1]/5;
			}
		}
		nsd = L - sumNID;
		
		gp = GroupingStrategy.overlapGrouping(nip, nop, nsp, numberOfObjectives_-1, 0);
		gd = GroupingStrategy.overlapGrouping(nid, nod, nsd, numberOfObjectives_, K);
		dgd = new int[numberOfObjectives_][][];
		for(int i=0;i<numberOfObjectives_;i++){
			if(i%2 == 0){
				dgd[i] = GroupingStrategy.deepGrouping(gd[i], 10, 3);
			}else{
				dgd[i] = GroupingStrategy.deepGrouping(gd[i], 10, 2);
			}
		}
		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		for (int var = 0; var < K; var++) {
			lowerLimit_[var] = -1.0;
			upperLimit_[var] = 1.0;
		} // for
		for (int var = K; var < numberOfVariables; var++) {
			lowerLimit_[var] = lg;
			upperLimit_[var] = ug;
		}

		shiftValues_ = new double[L];
		for (int i = 0; i < shiftValues_.length; i++)
			shiftValues_[i] = 0;

		rotationMatrix_ = new double[L][L];
		for (int i = 0; i < rotationMatrix_.length; i++) {
			for (int j = 0; j < rotationMatrix_.length; j++) {
				if (i != j)
					rotationMatrix_[i][j] = 0;
				else
					rotationMatrix_[i][j] = 1;
			}
		}
	}
	
	public MMLMOP(int numberOfObjectives, int numberOfVariables, int number_nip, double lg, double ug) {
		numberOfObjectives_ = numberOfObjectives;
		numberOfVariables_ = numberOfVariables;
	
		gType_ = "sphere";
		hType_ = "sphere";
		genType_ = "multiplication";
		linkageType_ = "linear";
		
		nip = new int[numberOfObjectives_-1];
		nop = new int[numberOfObjectives_-1];
		int sumNIP = 0;
		for(int i=0; i<numberOfObjectives_-1; i++){
			nip[i] = number_nip;
			if(numberOfObjectives_ == 2){
				nop[i] = 0;
			}else{
				nop[i] = 1;
			}
			sumNIP += nip[i];
		}
		nsp = 1;
		K = nsp + sumNIP;
		L = numberOfVariables_ - K;
		
		double[] c = new double[numberOfObjectives_];
		c[0] = 3.8*0.1*(1.0-0.1);
		double sumC = c[0];
		for(int i = 1; i < numberOfObjectives_; i++){
			c[i] = 3.8*c[i-1]*(1.0-c[i-1]);
			sumC += c[i];
		}
		nid = new int[numberOfObjectives_];
		int sumNID = 0;
		for(int i = 0; i < numberOfObjectives_; i++){
			nid[i] = (int)Math.floor((c[i]/sumC)*L);
			sumNID += nid[i];
		}
		nod = new int[numberOfObjectives_];
		for(int i = 0; i < numberOfObjectives_; i++){
			if(i == 0){
				nod[i] = nid[numberOfObjectives_-1]/5;
			}else{
				nod[i] = nid[i-1]/5;
			}
		}
		nsd = L - sumNID;
		
		gp = GroupingStrategy.overlapGrouping(nip, nop, nsp, numberOfObjectives_-1, 0);
		gd = GroupingStrategy.overlapGrouping(nid, nod, nsd, numberOfObjectives_, K);
		dgd = new int[numberOfObjectives_][][];
		for(int i=0;i<numberOfObjectives_;i++){
			if(i%2 == 0){
				dgd[i] = GroupingStrategy.deepGrouping(gd[i], 4, 2);
			}else{
				dgd[i] = GroupingStrategy.deepGrouping(gd[i], 4, 1);
			}
		}
		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		for (int var = 0; var < K; var++) {
			lowerLimit_[var] = -1.0;
			upperLimit_[var] = 1.0;
		} // for
		for (int var = K; var < numberOfVariables; var++) {
			lowerLimit_[var] = lg;
			upperLimit_[var] = ug;
		}

		shiftValues_ = new double[L];
		for (int i = 0; i < shiftValues_.length; i++)
			shiftValues_[i] = 0;

		rotationMatrix_ = new double[L][L];
		for (int i = 0; i < rotationMatrix_.length; i++) {
			for (int j = 0; j < rotationMatrix_.length; j++) {
				if (i != j)
					rotationMatrix_[i][j] = 0;
				else
					rotationMatrix_[i][j] = 1;
			}
		}
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		double vars[] = scaleVariables(solution);
		
		double[][] xI = new double[numberOfObjectives_ - 1][];
		double[][][] xII = new double[numberOfObjectives_][][];

		for (int i = 0; i < numberOfObjectives_ - 1; i++){
			xI[i] = new double[gp[i].length];
			for(int j=0;j<gp[i].length;j++){
				xI[i][j] = vars[gp[i][j]];
			}
		}
		
		
		double[] y = new double[numberOfObjectives_-1];
		for(int i=0;i<numberOfObjectives_-1;i++){
			y[i] = sum_avg_abs(xI[i]);
		}

		for (int i = 0; i < numberOfObjectives_; i++){
			xII[i] = new double[dgd[i].length][];
			for(int j=0;j<dgd[i].length;j++){
				xII[i][j] = new double[dgd[i][j].length];
				for(int k=0;k<dgd[i][j].length;k++){
					xII[i][j][k] = vars[dgd[i][j][k]];
				}
			}
		}
		
		double[] h = evalH(y);
		double[] g = evalG(vars[0],xII);
		
		if(genType_.equalsIgnoreCase("addition")){
			for (int i = 0; i < numberOfObjectives_; i++)
				solution.setObjective(startObjPos_ + i, h[i] + g[i]);
		}else if(genType_.equalsIgnoreCase("multiplication"))
			for (int i = 0; i < numberOfObjectives_; i++)
				solution.setObjective(startObjPos_ + i, h[i] * (1 + g[i]));
		else{
			System.out.println("Error: generator type " + genType_ + " invalid");
			System.exit(0);
		}
	}

	double[] evalG(double xI, double[][][] xIII) throws JMException {
		double[] g = null;
		if (gType_.equalsIgnoreCase("DF1"))
			g = GFunctions.getDF1(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF2"))
			g = GFunctions.getDF2(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF3"))
			g = GFunctions.getDF3(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF4"))
			g = GFunctions.getDF4(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF5"))
			g = GFunctions.getDF5(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF6"))
			g = GFunctions.getDF6(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF7"))
			g = GFunctions.getDF7(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF8"))
			g = GFunctions.getDF8(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF9"))
			g = GFunctions.getDF9(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF10"))
			g = GFunctions.getDF10(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF11"))
			g = GFunctions.getDF11(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF12"))
			g = GFunctions.getDF12(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF13"))
			g = GFunctions.getDF13(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF14"))
			g = GFunctions.getDF14(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF15"))
			g = GFunctions.getDF15(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF16"))
			g = GFunctions.getDF16(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF17"))
			g = GFunctions.getDF17(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF18"))
			g = GFunctions.getDF18(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF19"))
			g = GFunctions.getDF19(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF20"))
			g = GFunctions.getDF20(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF21"))
			g = GFunctions.getDF21(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF22"))
			g = GFunctions.getDF22(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF23"))
			g = GFunctions.getDF23(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF24"))
			g = GFunctions.getDF24(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else if (gType_.equalsIgnoreCase("DF25"))
			g = GFunctions.getDF25(xI, xIII,dgd,numberOfVariables_,linkageType_);
		else {
			System.out.println("Error: g function type " + gType_ + " invalid");
			System.exit(0);
		}
		return g;
	}

	double[] evalH(double[] xI) {
		double[] h = new double[numberOfObjectives_];
		if (hType_.equalsIgnoreCase("lineoid")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= xI[j];
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= (1 - xI[aux]);
				} // if
			} // for
		}else if(hType_.equalsIgnoreCase("circle" ) || hType_.equalsIgnoreCase("sphere")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
				} // if
			} // for
		}else if(hType_.equalsIgnoreCase("convex")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
				} // if
				if(i != numberOfObjectives_-1)
					h[i] = Math.pow(h[i], 4);
				else
					h[i] = Math.pow(h[i], 2);
			} // for
		}if (hType_.equalsIgnoreCase("inverted_lineoid")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= xI[j];
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= (1 - xI[aux]);
				} // if
				h[i] = 1.0 - h[i];
			} // for
		}else if(hType_.equalsIgnoreCase("inverted_sphere")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
				} // if
				h[i] = 1.0 - h[i];
			} // for
		}
		return h;
	}

	public void setGType(String gType) {
		gType_ = gType;
	}
	public void setGenType(String genType) {
		genType_ = genType;
	}
	
	public void setLinkageType(String linkageType) {
		linkageType_ = linkageType;
	}
	
	protected double sum_avg_abs(double[] x){
		double avg = 0;
		for(int i=0; i<x.length; i++){
			avg += x[i];
		}
		return Math.abs(avg)/x.length;
	}

	@Override
	public void dynamicEvaluate(Solution solution, int currentGeneration) throws JMException {
		
	}
}
