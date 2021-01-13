package etmo.problems.base.dynamicBase;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.util.JMException;

public class DMLZ extends Problem {

	Double t_;//time instant
    Double ht1_;//the time-related variable used to control the change of the PF
    double[] ht2_;
    Double gt_;//the time-related variable used to control the change of the PS
    
    int fc_;
    int sc_;

	String gType_;  //The type of the landscape function
	String ht1Type_; //The type of time-related function used to control the change of the PF
	String ht2Type_; //The type of time-related function used to control the change of the PF
	String gtType_; //The type of time-related function used to control the change of the PS 


	public DMLZ() {
		numberOfObjectives_ = 2;
	}

	/*
	 * gc indicates the generation counter
	 * fc indicates the frequency of change
	 * sc indicates the severity of change
	*/
	public DMLZ(int numberOfVariables, double lg, double ug, int fc, int sc) {
		numberOfObjectives_ = 2;
		numberOfVariables_ = numberOfVariables;
		
		fc_ = fc;
		sc_ = sc;
		
		gType_ = "sphere";
		hType_ = "convex";
		ht1Type_ = "static";
		ht2Type_ = "static";
		gtType_ = "static";

		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		for (int var = 0; var < 1; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
		} // for
		for (int var = 1; var < numberOfVariables; var++) {
			lowerLimit_[var] = lg;
			upperLimit_[var] = ug;
		}

		shiftValues_ = new double[numberOfVariables_ - 1];
		for (int i = 0; i < shiftValues_.length; i++)
			shiftValues_[i] = 0;

		rotationMatrix_ = new double[numberOfVariables_ - 1][numberOfVariables_ - 1];
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
	
		double[] xI = new double[1]; //The k position-related variables
		double[] xII = new double[numberOfVariables_ - 1]; //The n-k distance-related variables
		
		for (int i = 0; i < 1; i++)
			xI[i] = vars[i];
		for (int i = 1; i < numberOfVariables_; i++)
			xII[i - 1] = vars[i];
		xII = transformVariables(xII);
		
		gt_ = evalGt();
		ht1_ = evalHt1(xI[0]);
		ht2_ = evalHt2();
		
		double[] f = new double[2];
		double[] h = evalH(xI[0]);
		double g = evalG(xII, xI[0]);

		for (int i = 0; i < 2; i++)
			f[i] = (1 + g) * h[i];

		for (int i = 0; i < 2; i++)
			solution.setObjective(startObjPos_ + i, f[i]);
	}
	
	double[] evalH(double xI) {
		double[] h = new double[2];
		h[0] = xI + ht1_;
		h[0] = Math.pow(h[0],ht2_[0]);
		h[1] = 1 - xI + ht1_;
		h[1] = Math.pow(h[1],ht2_[1]);
		return h;
	}

	double evalG(double[] xII, double x1) throws JMException {
		if (gType_.equalsIgnoreCase("sphere"))
			return DGFunctions.getSphere(xII, gt_);
		else if (gType_.equalsIgnoreCase("rastrigin"))
			return DGFunctions.getRastrigin(xII, gt_);
		else if (gType_.equalsIgnoreCase("griewank"))
			return DGFunctions.getGriewank(xII, gt_);
		else if (gType_.equalsIgnoreCase("ackley"))
			return DGFunctions.getAckley(xII, gt_);
		else if (gType_.equalsIgnoreCase("rosenbrock"))
			return DGFunctions.getRosenbrock(xII, gt_);
		else if (gType_.equalsIgnoreCase("sphere_linkage5"))
			return DGFunctions.getSphereLinkage5(xII, gt_, x1);
		else {
			System.out.println("Error: G function type " + gType_ + " invalid");
			return Double.NaN;
		}
	}
	
	double evalHt1(double x1){
		if (ht1Type_.equalsIgnoreCase("static"))
			return 0.1*Math.sin(3*Math.PI*x1);
		else if (ht1Type_.equalsIgnoreCase("dynamic_t1"))
			return 0.02*Math.sin(Math.floor(10*gt_)*Math.PI*x1);
		else if (ht1Type_.equalsIgnoreCase("dynamic_t2"))
			return 1.5+Math.sin(0.5*Math.PI*t_);
		else
			return 1.0;
	}
	
	double[] evalHt2(){
		double[] at = new double[2];
		if (ht2Type_.equalsIgnoreCase("static")){
			at[0] = at[1] = 1.0;
		}else if (ht2Type_.equalsIgnoreCase("dynamic_t1")){
			at[0] = 1.0;
			at[1] = 2.25 + 2*Math.cos(2*Math.PI*t_);
		}else if (ht2Type_.equalsIgnoreCase("dynamic_t2"))
			at[0] = at[1] = 0.2 + 2.8*Math.abs(gt_);
		else{
			System.out.println("Error: HT2 type " + gType_ + " invalid");
			System.exit(0);
		}
		return at;
	}
	
	double evalGt(){
		if (gtType_.equalsIgnoreCase("static"))
			return 1.0;
		else if (gtType_.equalsIgnoreCase("dynamic_t1"))
			return Math.abs(Math.sin(0.5*Math.PI*t_));
		else if (gtType_.equalsIgnoreCase("dynamic_t2"))
			return Math.sin(0.5*Math.PI*t_);
		else
			return 1.0;
	}

	public void setGType(String gType) {
		gType_ = gType;
	}

	public void setHType(String hType) {
		hType_ = hType;
	}

	public void setHt1Type(String htType) {
		ht1Type_ = htType;
	}
	
	public void setHt2Type(String htType) {
		ht2Type_ = htType;
	}

	public void setGtType(String gtType) {
		gtType_ = gtType;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		// TODO Auto-generated method stub
		
	}

}
