package etmo.problems.base.dynamicBase;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.util.JMException;

public class DMDTLZ extends Problem {

	Double t_;//time instant
    Double ht_;//the time-related variable used to control the change of the PF
    Double gt_;//the time-related variable used to control the change of the PS

	String gType_;  //The type of the landscape function
	String htType_; //The type of time-related function used to control the change of the PF
	String gtType_; //The type of time-related function used to control the change of the PS
	String pvType_; 
	
	int fc_;
    int sc_;


	public DMDTLZ() {
		numberOfObjectives_ = 3;
	}

	/*
	 * gc indicates the generation counter
	 * fc indicates the frequency of change
	 * sc indicates the severity of change
	*/
	public DMDTLZ(int numberOfObjectives, int numberOfVariables, double lg, double ug, int fc, int sc) {
		numberOfObjectives_ = numberOfObjectives;
		numberOfVariables_ = numberOfVariables;
		
		int num = numberOfVariables_ - numberOfObjectives_ + 1;
		
		fc_ = fc;
		sc_ = sc;
	
		pvType_ = "static";
		gType_ = "sphere";
		hType_ = "normal";
		htType_ = "static";
		gtType_ = "static";

		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		for (int var = 0; var < numberOfObjectives_-1; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
		} // for
		for (int var = numberOfObjectives_-1; var < numberOfVariables; var++) {
			lowerLimit_[var] = lg;
			upperLimit_[var] = ug;
		}

		shiftValues_ = new double[num];
		for (int i = 0; i < shiftValues_.length; i++)
			shiftValues_[i] = 0;

		rotationMatrix_ = new double[num][num];
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
	public void dynamicEvaluate(Solution solution, int gc) throws JMException {
		double vars[] = scaleVariables(solution);
		t_ = (1.0/sc_)*Math.floor(gc/fc_);
		gt_ = evalGt();
		ht_ = evalHt();

		double[] xI = new double[numberOfObjectives_-1]; //The m-1 position-related variables
		double[] xII = new double[numberOfVariables_ - numberOfObjectives_ + 1]; //The n-m+1 distance-related variables
		
	    double[] pv = new double[numberOfObjectives_-1];
		
		if(pvType_.equalsIgnoreCase("static")){
			for (int i = 0; i < numberOfObjectives_ - 1; i++){
				xI[i] = vars[i];
			}
				
		}else if(pvType_.equalsIgnoreCase("dynamic")){
			for (int i = 0; i < numberOfObjectives_ - 1; i++)
				xI[i] = evalPV(vars[i]);
		}else{
			System.out.println("Error: position-related variable type " + pvType_ + " invalid");
			System.exit(0);
		}
		for (int i = numberOfObjectives_ - 1; i < numberOfVariables_; i++)
			xII[i - numberOfObjectives_ + 1] = vars[i];	
		
		for (int i = 0; i < numberOfObjectives_ - 1; i++){
			pv[i] = vars[i];
		}
		
		double[] f = new double[numberOfObjectives_];
		
		double[] h = evalH(xI);

		double g = evalG(xII, pv);
		
		for (int i = 0; i < numberOfObjectives_; i++)
			f[i] = (1 + g) * h[i];

		for (int i = 0; i < numberOfObjectives_; i++)
			solution.setObjective(startObjPos_ + i, f[i]);
	}

	double evalG(double[] xII, double[] pv) throws JMException {
		if (gType_.equalsIgnoreCase("sphere"))
			return DGFunctions.getSphere(xII, gt_);
		else if(gType_.equalsIgnoreCase("sphere_linkage2"))
			return DGFunctions.getSphereLinkage2(xII, pv, gt_);
		else if(gType_.equalsIgnoreCase("sphere_linkage3"))
			return DGFunctions.getSphereLinkage3(xII, pv, gt_);
		else if(gType_.equalsIgnoreCase("sphere_linkage4"))
			return DGFunctions.getSphereLinkage4(xII, pv, gt_);
		else {
			System.out.println("Error: G function type " + gType_ + " invalid");
			return Double.NaN;
		}
	}

	double[] evalH(double[] xI) {
		double[] h = new double[numberOfObjectives_];
		
		if (hType_.equalsIgnoreCase("normal")){
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
				} // if
			} // for
		}else if (hType_.equalsIgnoreCase("inverse")){	
			for (int i = numberOfObjectives_-1; i >=0; i--) {
				h[i] = 1.0;
				for (int j = numberOfObjectives_-(i+1)-1; j >= 0; j--)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != numberOfObjectives_-1) {
					h[i] *= Math.sin(xI[i] * 0.5 * Math.PI);
				} // if
			} // for
		}else {
			System.out.println("Error: H function type " + hType_ + " invalid");
			System.exit(0);
		}
		
		for (int i = 0; i < numberOfObjectives_; i++) {
			h[i] = Math.pow(h[i], ht_);
		}
		
		return h;
	}
	
	double evalPV(double x){
		double px;
		px = (Math.PI/6)*gt_ + (Math.PI/2 - (Math.PI/3)*gt_)*x;
		return px;
	}

	
	double evalHt(){
		if (htType_.equalsIgnoreCase("static"))
			return 1.0;
		else if (htType_.equalsIgnoreCase("dynamic_t1"))
			return 2.25 + 2.0*Math.cos(0.5*Math.PI*t_);
		else if (htType_.equalsIgnoreCase("dynamic_t2"))
			return 1.5+Math.sin(0.5*Math.PI*t_);
		else
			return 1.0;
	}
	
	double evalGt(){
		if (gtType_.equalsIgnoreCase("static"))
			return 1.0;
		else if (gtType_.equalsIgnoreCase("dynamic_t1"))
			return Math.sin(0.5*Math.PI*t_);
		else if (gtType_.equalsIgnoreCase("dynamic_t2"))
			return Math.abs(Math.sin(0.5*Math.PI*t_));
		else if (gtType_.equalsIgnoreCase("dynamic_t3"))
			return t_;
		else
			return 1.0;
	}

	public void setGType(String gType) {
		gType_ = gType;
	}
	
	public void setHType(String hType) {
		hType_ = hType;
	}
	
	public void setHtType(String htType) {
		htType_ = htType;
	}

	public void setGtType(String gtType) {
		gtType_ = gtType;
	}
	
	public void setPvType(String pvType) {
		pvType_ = pvType;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		// TODO Auto-generated method stub
		
	}
}
