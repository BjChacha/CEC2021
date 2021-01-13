package etmo.problems.base.staticBase;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.util.JMException;

public class MMLZ extends Problem {
	String gType_;
	String genType_;

	public MMLZ() {
		numberOfObjectives_ = 2;
	}

	public MMLZ(int numberOfObjectives, int numberOfVariables, double lg, double ug) {
		numberOfObjectives_ = numberOfObjectives;
		numberOfVariables_ = numberOfVariables;
	
		gType_ = "LF1";
		hType_ = "convex";
		genType_ = "addition";
		
		int num = numberOfVariables_ - numberOfObjectives_ + 1;
		
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
	public void evaluate(Solution solution) throws JMException {
		double vars[] = scaleVariables(solution);

		double[] xI = new double[numberOfObjectives_ - 1];
		double[] xII = new double[numberOfVariables_ - numberOfObjectives_ + 1];

		for (int i = 0; i < numberOfObjectives_ - 1; i++)
			xI[i] = vars[i];
		
		for (int i = numberOfObjectives_ - 1; i < numberOfVariables_; i++)
			xII[i - numberOfObjectives_ + 1] = vars[i];
		
		xII = transformVariables(xII);
		
		int quotient = (numberOfVariables_ - numberOfObjectives_ + 1)/numberOfObjectives_;
		int remainder = (numberOfVariables_ - numberOfObjectives_ + 1)%numberOfObjectives_;

		double[][] xIII = new double[numberOfObjectives_][];
		int[][] index = new int[numberOfObjectives_][];
		for(int i=0; i<numberOfObjectives_; i++){
			if(remainder > i){
				xIII[i] = new double[quotient+1];
				index[i] = new int[quotient+1];
			}else{
				xIII[i] = new double[quotient];
				index[i] = new int[quotient];
			}
			int t = 0;
			for (int j = i; j < xII.length; ){
				xIII[i][t] = xII[j];
				index[i][t] = j + numberOfObjectives_;
				t++;
				j += numberOfObjectives_;
			}	
		}
		double[] h = evalH(xI);
		double[] g = evalG(xI,xIII,index);
		
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

	double[] evalG(double[] xI, double[][] xIII, int[][] index) throws JMException {
		double[] g = null;
		if (gType_.equalsIgnoreCase("LF1"))
			g = GFunctions.getLF1(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("LF2"))
			g = GFunctions.getLF2(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("LF3"))
			g = GFunctions.getLF3(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("LF4_5"))
			g = GFunctions.getLF4_5(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("LF6"))
			g = GFunctions.getLF6(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("LF7"))
			g = GFunctions.getLF7(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("LF8"))
			g = GFunctions.getLF8(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("LF9"))
			g = GFunctions.getLF9(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("LF10"))
			g = GFunctions.getLF10(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("F14"))
			g = GFunctions.getF14(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("F15"))
			g = GFunctions.getF15(xI, xIII, index, numberOfVariables_);
		else if (gType_.equalsIgnoreCase("F16"))
			g = GFunctions.getF16(xI, xIII, index, numberOfVariables_);
		else {
			System.out.println("Error: g function type " + gType_ + " invalid");
			System.exit(0);
		}
		return g;
	}

	double[] evalH(double[] xI) {
		double[] h = new double[numberOfObjectives_];
		if(numberOfObjectives_ == 2){
			if (hType_.equalsIgnoreCase("convex")){
				h[0] = xI[0];
				h[1] = 1 - Math.pow(xI[0], 0.5);
			}else if(hType_.equalsIgnoreCase("concave")){
				h[0] = xI[0];
				h[1] = 1 - Math.pow(xI[0], 2.0);
			}else if(hType_.equalsIgnoreCase("lineoid")){
				h[0] = xI[0];
				h[1] = 1 - xI[0];
			}else {
				System.out.println("Error: H function type " + hType_ + " invalid");
				System.exit(0);
			}
		}else if(numberOfObjectives_ >= 3){
			if (hType_.equalsIgnoreCase("sphere")){
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
			}else if(hType_.equalsIgnoreCase("lineoid")){
				for (int i = 0; i < numberOfObjectives_; i++) {
					h[i] = 1.0;
					for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
						h[i] *= xI[j];
					if (i != 0) {
						int aux = numberOfObjectives_ - (i + 1);
						h[i] *= (1 - xI[aux]);
					} // if
				} // for
			}
				
		}else{
			System.out.println("Error: numberOfObjectives =  " + numberOfObjectives_ + " is invalid");
			System.exit(0);
		}	
		return h;
	}

	public void setGType(String gType) {
		gType_ = gType;
	}
	public void setGenType(String genType) {
		genType_ = genType;
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
