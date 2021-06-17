package etmo.problems.base.dynamicBase;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.util.JMException;

public class DMZDT extends Problem {

	Integer k_; //number of position-related variables
	Double t_;//time instant
    Double ht_;//the time-related variable used to control the change of the PF
    Double gt_;//the time-related variable used to control the change of the PS
    
    int fc_;
    int sc_;

	String gType_;  //The type of the landscape function
	String f1Type_; //The type of processing the position-related variables for f1 in ZDT-based formulation
	String htType_; //The type of time-related function used to control the change of the PF
	String gtType_; //The type of time-related function used to control the change of the PS
	String pvType_; 


	public DMZDT() {
		numberOfObjectives_ = 2;
	}

	/*
	 * gc indicates the generation counter
	 * fc indicates the frequency of change
	 * sc indicates the severity of change
	*/
	public DMZDT(int numberOfVariables, int k, double lg, double ug, int fc, int sc) {
		numberOfObjectives_ = 2;
		numberOfVariables_ = numberOfVariables;
		k_ = k;
		
		fc_ = fc;
		sc_ = sc;
		
		pvType_ = "static";
		gType_ = "sphere";
		f1Type_ = "linear";
		hType_ = "convex";
		htType_ = "constant";
		gtType_ = "constant";

		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		for (int var = 0; var < k_; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
		} // for
		for (int var = k_; var < numberOfVariables; var++) {
			lowerLimit_[var] = lg;
			upperLimit_[var] = ug;
		}

		shiftValues_ = new double[numberOfVariables_ - k_];
		for (int i = 0; i < shiftValues_.length; i++)
			shiftValues_[i] = 0;

		rotationMatrix_ = new double[numberOfVariables_ - k_][numberOfVariables_ - k_];
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
		
		ht_ = evalHt();
		gt_ = evalGt();
		
		double[] xI = new double[k_]; //The k position-related variables
		double[] xII = new double[numberOfVariables_ - k_]; //The n-k distance-related variables
		
		if(pvType_.equalsIgnoreCase("static")){
			for (int i = 0; i < k_; i++)
				xI[i] = vars[i];
			for (int i = k_; i < numberOfVariables_; i++)
				xII[i - k_] = vars[i];
		}else if(pvType_.equalsIgnoreCase("dynamic")){
			int pos = (int) Math.floor((numberOfVariables_-k_)*gt_);
			int iter = pos;
			for (int i = 0; i < k_; i++){
				xI[i] = vars[iter];
				iter++;
			}
			iter = 0;
			for (int i = 0; i < numberOfVariables_-k_; i++){
				if(iter != pos){
					xII[i] = vars[iter];
					iter++;
				}else{
					iter += k_;
					xII[i] = vars[iter];
				}
			}	
		}else{
			System.out.println("Error: position-related variable type " + pvType_ + " invalid");
			System.exit(0);
		}
		
		
		
		xII = transformVariables(xII);
		double f1 = evalF1(xI);
		double g = evalG(xII, f1) + 1;
		double f2 = g * evalH(f1, g);

		solution.setGFunValue(g);
		// System.out.println("g: " + g);

		solution.setObjective(startObjPos_, f1);
		solution.setObjective(startObjPos_ + 1, f2);
	}

	double evalG(double[] xII, double f1) throws JMException {//Rastrigin
		if (gType_.equalsIgnoreCase("sphere"))
			return DGFunctions.getSphere(xII, gt_);
		if (gType_.equalsIgnoreCase("rosenbrock"))
			return DGFunctions.getRosenbrock(xII, gt_);
		if (gType_.equalsIgnoreCase("rastrigin"))
			return DGFunctions.getRastrigin(xII, gt_);
		else if (gType_.equalsIgnoreCase("sphere_linkage"))
			return DGFunctions.getSphereLinkage1(xII, gt_, f1, ht_);
		else {
			System.out.println("Error: G function type " + gType_ + " invalid");
			return Double.NaN;
		}
	}

	double evalF1(double[] xI) {
		if (f1Type_.equalsIgnoreCase("linear"))
			return F1_linear(xI);
		else if (f1Type_.equalsIgnoreCase("nonlinear"))
			return F1_nonlinear(xI);
		else {
			System.out.println("Error: f1 function type " + f1Type_ + " invalid");
			return Double.NaN;
		}
	}

	double evalH(double f1, double g) {
		if (hType_.equalsIgnoreCase("static_convex"))
			return H_convex(f1, g);
		else if (hType_.equalsIgnoreCase("static_concave"))
			return H_nonconvex(f1, g);
		else if (hType_.equalsIgnoreCase("dynamic"))
			return H_dynamic(f1, g);
		else {
			System.out.println("Error: H function type " + hType_ + " invalid");
			return Double.NaN;
		}
	}

	double H_convex(double f1, double g) {
		return 1 - Math.pow(f1 / g, 0.5);
	}

	double H_nonconvex(double f1, double g) {
		return 1 - Math.pow(f1 / g, 2);
	}
	
	double H_dynamic(double f1, double g) {
		return 1 - Math.pow(f1 / g, ht_);
	}
	
	double evalHt(){
		if (htType_.equalsIgnoreCase("static"))
			return 1.0;
		else if (htType_.equalsIgnoreCase("dynamic_t1"))
			return 0.75*Math.sin(0.5*Math.PI*t_);
		else if (htType_.equalsIgnoreCase("dynamic_t2"))
			return 1.5+Math.sin(0.5*Math.PI*t_);
		else
			return 1.0;
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

	double F1_linear(double xI[]) {
		double sum = 0;
		for (int i = 0; i < xI.length; i++)
			sum += xI[i];

		return sum / xI.length;
	}

	double F1_nonlinear(double xI[]) {
		double r = 0;

		for (int i = 0; i < xI.length; i++)
			r += (xI[i] * xI[i]);

		r = Math.sqrt(r);

		return 1 - Math.exp(-4 * r) * Math.pow(Math.sin(5 * Math.PI * r), 4);
	}

	public void setGType(String gType) {
		gType_ = gType;
	}

	public void setF1Type(String f1Type) {
		f1Type_ = f1Type;
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
