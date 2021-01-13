package etmo.problems.base.staticBase;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.util.JMException;

public class MMZDT extends Problem {

	Integer k_;


	String gType_;
	String f1Type_;


	public MMZDT() {
		numberOfObjectives_ = 2;
	}

	public MMZDT(int numberOfVariables, int k, double lg, double ug) {
		numberOfObjectives_ = 2;
		numberOfVariables_ = numberOfVariables;
		k_ = k;
	
		gType_ = "sphere";
		f1Type_ = "linear";
		hType_ = "convex";

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

	public MMZDT(int numberOfVariables, int k, double lg, double ug, String gType, String f1Type, String hType,
			double[] shiftValues, double[][] rotationMatrix) {
		numberOfObjectives_ = 2;
		numberOfVariables_ = numberOfVariables;
		k_ = k;

		gType_ = gType;
		f1Type_ = f1Type;
		hType_ = hType;

		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;

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

	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		double vars[] = scaleVariables(solution);

		double[] xI = new double[k_];
		double[] xII = new double[numberOfVariables_ - k_];

		for (int i = 0; i < k_; i++)
			xI[i] = vars[i];

		for (int i = k_; i < numberOfVariables_; i++)
			xII[i - k_] = vars[i];

		//xII = transformVariables(xII);

		double f1 = evalF1(xI);
		double g = evalG(xII) + 1;
		double f2 = g * evalH(f1, g);

		solution.setGFunValue(g);
		// System.out.println("g: " + g);

		solution.setObjective(startObjPos_, f1);
		solution.setObjective(startObjPos_ + 1, f2);
	}

	double evalG(double[] xII) throws JMException {
		if(gType_.equalsIgnoreCase("F1"))
			return GFunctions.getF1(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F2"))
			return GFunctions.getF2(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F3"))
			return GFunctions.getF3(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F4"))
			return GFunctions.getF4(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F5"))
			return GFunctions.getF5(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F6"))
			return GFunctions.getF6(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F7"))
			return GFunctions.getF7(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F8"))
			return GFunctions.getF8(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F9"))
			return GFunctions.getF9(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F10"))
			return GFunctions.getF10(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F11"))
			return GFunctions.getF11(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F12"))
			return GFunctions.getF12(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F13"))
			return GFunctions.getF13(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("HF1"))
			return GFunctions.getHF1(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("HF2"))
			return GFunctions.getHF2(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("HF3"))
			return GFunctions.getHF3(xII, shiftValues_, rotationMatrix_);		
		else if (gType_.equalsIgnoreCase("HF4"))
			return GFunctions.getHF4(xII, shiftValues_, rotationMatrix_);			
		else if (gType_.equalsIgnoreCase("HF5"))
			return GFunctions.getHF5(xII, shiftValues_, rotationMatrix_);		
		else {
			System.out.println("Error: g function type " + gType_ + " invalid");
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
		if (hType_.equalsIgnoreCase("convex"))
			return H_convex(f1, g);
		else if (hType_.equalsIgnoreCase("concave"))
			return H_nonconvex(f1, g);
		else {
			System.out.println("Error: f1 function type " + f1Type_ + " invalid");
			return Double.NaN;
		}
	}

	double H_convex(double f1, double g) {
		return 1 - Math.pow(f1 / g, 0.5);
	}

	double H_nonconvex(double f1, double g) {
		return 1 - Math.pow(f1 / g, 2);
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

	@Override
	public void dynamicEvaluate(Solution solution, int currentGeneration) throws JMException {
		double vars[] = scaleVariables(solution);

		double[] xI = new double[k_];
		double[] xII = new double[numberOfVariables_ - k_];

		for (int i = 0; i < k_; i++)
			xI[i] = vars[i];

		for (int i = k_; i < numberOfVariables_; i++)
			xII[i - k_] = vars[i];

		//xII = transformVariables(xII);

		double f1 = evalF1(xI);
		double g = evalG(xII) + 1;
		double f2 = g * evalH(f1, g);

		solution.setGFunValue(g);
		// System.out.println("g: " + g);

		solution.setObjective(startObjPos_, f1);
		solution.setObjective(startObjPos_ + 1, f2);
		
	}
}
