package etmo.problems.base.staticBase;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.util.JMException;

public class MMDTLZ extends Problem {
	String gType_;
	String genType_;
	Integer alpha_;

	
	public MMDTLZ(int numberOfObjectives, int numberOfVariables, int alpha, double lg, double ug) {
		numberOfObjectives_ = numberOfObjectives;
		numberOfVariables_ = numberOfVariables;

		if (numberOfObjectives == 2)
			hType_ = "circle";
		else
			hType_ = "sphere";
		gType_ = "sphere";
		genType_ = "multiplication";

		alpha_ = alpha;

		int num = numberOfVariables_ - numberOfObjectives_ + 1;

		// System.out.println(num);

		shiftValues_ = new double[num];
		rotationMatrix_ = new double[num][num];

		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		for (int var = 0; var < numberOfObjectives_ - 1; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
		} // for

		for (int var = numberOfObjectives_ - 1; var < numberOfVariables; var++) {
			lowerLimit_[var] = lg;
			upperLimit_[var] = ug;
		}

		for (int i = 0; i < num; i++)
			shiftValues_[i] = 0;

		for (int i = 0; i < num; i++) {
			for (int j = 0; j < num; j++) {
				if (i != j)
					rotationMatrix_[i][j] = 0;
				else
					rotationMatrix_[i][j] = 1;
			}
		}
	}

	public void evaluate(Solution solution) throws JMException {
		double vars[] = scaleVariables(solution);

		double[] xI = new double[numberOfObjectives_ - 1];
		double[] xII = new double[numberOfVariables_ - numberOfObjectives_ + 1];

		for (int i = 0; i < numberOfObjectives_ - 1; i++)
			xI[i] = vars[i];

		for (int i = numberOfObjectives_ - 1; i < numberOfVariables_; i++)
			xII[i - numberOfObjectives_ + 1] = vars[i];	

		double[] f = new double[numberOfObjectives_];

		double g = evalG(xII);
		
		double[] h = evalH(xI, g);
		
		if(genType_.equalsIgnoreCase("addition")){
			for (int i = 0; i < numberOfObjectives_; i++)
				solution.setObjective(startObjPos_ + i, h[i] + g);
		}else if(genType_.equalsIgnoreCase("multiplication"))
			for (int i = 0; i < numberOfObjectives_; i++)
				solution.setObjective(startObjPos_ + i, h[i] * (1 + g));
		else{
			System.out.println("Error: generator type " + genType_ + " invalid");
			System.exit(0);
		}
	}
	
	double[] evalH(double[] xI, double g) {
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
		}else if(hType_.equalsIgnoreCase("degenerate")){
			degenerateEvalPV(xI, g);
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
				} // if
			} // for
		}else if(hType_.equalsIgnoreCase("irconcave")){
			shiftEvalPV(xI);
			for (int i = 0; i < numberOfObjectives_; i++) {
				h[i] = 1.0;
				for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
					h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
				if (i != 0) {
					int aux = numberOfObjectives_ - (i + 1);
					h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
				} // if
			} // for
		}else if(hType_.equalsIgnoreCase("disconnect")){
			h[numberOfObjectives_-1] = 0;
			for (int i = 0; i < numberOfObjectives_-1; i++) {
				if(genType_.equalsIgnoreCase("multiplication"))
					h[i] = xI[i]/(1+g);
				else if(genType_.equalsIgnoreCase("addition"))
					h[i] = xI[i] - g;
				else{
					System.out.println("Error: Generation type " + genType_ + " invalid");
					System.exit(0);
				}	
				h[numberOfObjectives_-1] += (xI[i]*(1+Math.sin(3*Math.PI*xI[i])))/(1+g);
			} // for
			h[numberOfObjectives_-1] = numberOfObjectives_ - h[numberOfObjectives_-1];
		}else {
			System.out.println("Error: H function type " + hType_ + " invalid");
			System.exit(0);
		}
		return h;
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
		else if (gType_.equalsIgnoreCase("F17"))
			return GFunctions.getF17(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F18"))
			return GFunctions.getF18(xII, shiftValues_, rotationMatrix_);
		else if (gType_.equalsIgnoreCase("F19"))
			return GFunctions.getF19(xII, shiftValues_, rotationMatrix_);
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
		else if (gType_.equalsIgnoreCase("HF6"))
			return GFunctions.getHF6(xII, shiftValues_, rotationMatrix_);		
		else if (gType_.equalsIgnoreCase("HF7"))
			return GFunctions.getHF7(xII, shiftValues_, rotationMatrix_);			
		else if (gType_.equalsIgnoreCase("HF8"))
			return GFunctions.getHF8(xII, shiftValues_, rotationMatrix_);
		else {
			System.out.println("Error: g function type " + gType_ + " invalid");
			return Double.NaN;
		}
	}
	
	protected static void degenerateEvalPV(double[] xI, double g){
		int I = 2;
		for(int i=I-1;i<xI.length;i++){
			xI[i] = 0.5*(1.0+2*g*xI[i])/(1.0+g);
		}
	}

	public void setGType(String gType) {
		gType_ = gType;
	}
	
	public void setGenType(String genType) {
		genType_ = genType;
	}

	public String getHType() {
		return hType_;
	}
	
	protected static void shiftEvalPV(double[] xI){
		for(int i=0; i<xI.length; i++){
			xI[i] = 0.5*xI[i] + 0.25;
		}
	}

	@Override
	public void dynamicEvaluate(Solution solution, int currentGeneration) throws JMException {
		
	}

}
