package etmo.problems.base.staticBase;

public class GFunctions {

	protected static double[] shiftValues_;
	protected static double[][] rotationMatrix_;	

	/*F1: Shifted and Rotated Sphere Function*/
	public static double getF1(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);				
		return BaseFunctions.getSphere(x)+Fstar;
	}
	
	/*F2: Shifted and Rotated High Conditioned Elliptic Function*/
	public static double getF2(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);				
		return BaseFunctions.getElliptic(x)+Fstar;
	}
	
	/*F3: Shifted and Rotated Bent Cigar Function*/
	public static double getF3(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);				
		return BaseFunctions.getCigar(x)+Fstar;
	}
	
	/*F4: Shifted and Rotated Discus Function*/
	public static double getF4(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);				
		return BaseFunctions.getDiscus(x)+Fstar;
	}
	
	/*F5: Shifted and Rotated Rosenbrock Function*/	
	public static double getF5(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		
		shiftVariables(x);
		x=rotateVariables(x);			
		return BaseFunctions.getRosenbrock(x)+Fstar;
	}
	
	/*F6: Shifted and Rotated Ackley Function*/
	public static double getF6(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);				
		return BaseFunctions.getAckley(x)+Fstar;
	}
	
	/*F7: Shifted and Rotated Weierstrass Function*/
	public static double getF7(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);				
		return BaseFunctions.getWeierstrass(x)+Fstar;
	}
	
	/*F8: Shifted and Rotated Griewank Function*/
	public static double getF8(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);				
		return BaseFunctions.getGriewank(x)+Fstar;
	}
			
	/*F9: Shifted and Rotated Rastrigin Function*/	
	public static double getF9(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);		
		return BaseFunctions.getRastrigin(x)+Fstar;
	}		
	
	/*F10: Shifted and Rotated Schwefel Function	*/	
	public static double getF10(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);
		
		return BaseFunctions.getMSchwefel(x)+Fstar;
	}	
	
	/*F11: Shifted and Rotated Katsuura Function*/	
	public static double getF11(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);
		
		return BaseFunctions.getKatsuura(x)+Fstar;
	}	

	/*F12: Shifted and Rotated HappyCat Function*/	
	public static double getF12(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);
		
		return BaseFunctions.getHappyCat(x)+Fstar;
	}		
	
	/*F13: Shifted and Rotated Expanded Griewank plus Rosenbrock Function*/
	public static double getF13(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);		
		return BaseFunctions.getExGriewRosen(x)+Fstar;
	}	
	
	/*F14: Mean Function*/
	public static double[] getF14(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				g[i] += 9*Math.abs(xIII[i][j]);
			}
			g[i] = g[i]/xIII[i].length;
		}
		return g;
	}
	
	/*F15: GFunction Used in*/
	public static double[] getF15(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		for(int i=0; i<numObj; i++){
			g[i] = BaseFunctions.getRastrigin(xIII[i]);
		}
		return g;
	}
	
	/*F16: Mean Function*/
	public static double[] getF16(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				g[i] += (0.5*xIII[i][j]-0.25)*(0.5*xIII[i][j]-0.25);
			}
		}
		return g;
	}
	
	
	public static double getF17(double[] x, double[] shiftValues, double[][] rotationMatrix) {
		double g=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);
		for(int j=0; j<x.length; j++){
			g += 9*Math.abs(x[j]);
		}
		g = 1 + g/x.length;
		return g;
	}
	
	public static double getF18(double[] x, double[] shiftValues, double[][] rotationMatrix) {
		double g=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);
		return 1 + BaseFunctions.getSphere(x);
	}
	
	/*F19: GFunction Used in*/
	public static double getF19(double[] x, double[] shiftValues, double[][] rotationMatrix) {
		double g=0.0;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);
		return 1 + BaseFunctions.getRastrigin(x);
	}
	
	
	
	/*HF1: Hybrid Function 1*/
	public static double getHF1(double x[], double[] shiftValues, double[][] rotationMatrix) {
		double Fstar=0.0;
		int dim = x.length;
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);	
		double[] p = {0.3,0.3,0.4};
		int n1 = (int)Math.ceil(p[0]*dim);			
		int n2 = (int)Math.ceil(p[1]*dim);
		int n3 = dim-n1-n2;
		
		double[] xI = new double[n1];
		double[] xII = new double[n2];
		double[] xIII = new double[n3];
		
		for (int i = 0; i < n1; i++)
			xI[i]=x[i];
		for (int i = n1; i < n1+n2; i++)
			xII[i-n1]=x[i];		
		for (int i = n1+n2; i < dim; i++)
			xIII[i-n1-n2]=x[i];		
		return BaseFunctions.getCigar(xI) + BaseFunctions.getRastrigin(xII) + BaseFunctions.getMSchwefel(xIII) + Fstar;
	}	

	/*HF2: Hybrid Function 2*/
	public static double getHF2(double x[], double[] shiftValues, double[][] rotationMatrix) {
		
		double Fstar=0.0;
		int dim = x.length;
		
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);		
		
		double[] p = {0.3,0.3,0.4};
		int n1 = (int)Math.ceil(p[0]*dim);			
		int n2 = (int)Math.ceil(p[1]*dim);
		int n3 = dim-n1-n2;
		
		double[] xI = new double[n1];
		double[] xII = new double[n2];
		double[] xIII = new double[n3];
		
		for (int i = 0; i < n1; i++)
			xI[i]=x[i];
		for (int i = n1; i < n1+n2; i++)
			xII[i-n1]=x[i];		
		for (int i = n1+n2; i < dim; i++)
			xIII[i-n1-n2]=x[i];		
			
		return BaseFunctions.getRosenbrock(xI)+BaseFunctions.getWeierstrass(xII)
				+BaseFunctions.getGriewank(xIII)+Fstar;
	}	
	
	/*HF3: Hybrid Function 3*/
	public static double getHF3(double x[], double[] shiftValues, double[][] rotationMatrix) {
		
		double Fstar=0.0;
		int dim = x.length;
		
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);		
		
		double[] p = {0.3,0.3,0.4};
		int n1 = (int)Math.ceil(p[0]*dim);			
		int n2 = (int)Math.ceil(p[1]*dim);
		int n3 = dim-n1-n2;
		
		double[] xI = new double[n1];
		double[] xII = new double[n2];
		double[] xIII = new double[n3];
		
		for (int i = 0; i < n1; i++)
			xI[i]=x[i];
		for (int i = n1; i < n1+n2; i++)
			xII[i-n1]=x[i];		
		for (int i = n1+n2; i < dim; i++)
			xIII[i-n1-n2]=x[i];		
			
		return BaseFunctions.getDiscus(xI)+BaseFunctions.getAckley(xII)
				+BaseFunctions.getGriewank(xIII)+Fstar;
	}	

	
	/*HF4: Hybrid Function 4*/
	public static double getHF4(double x[], double[] shiftValues, double[][] rotationMatrix) {
		
		double Fstar=0.0;
		int dim = x.length;
		
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);		
		
		double[] p = {0.3,0.3,0.4};
		int n1 = (int)Math.ceil(p[0]*dim);			
		int n2 = (int)Math.ceil(p[1]*dim);
		int n3 = dim-n1-n2;
		
		double[] xI = new double[n1];
		double[] xII = new double[n2];
		double[] xIII = new double[n3];
		
		for (int i = 0; i < n1; i++)
			xI[i]=x[i];
		for (int i = n1; i < n1+n2; i++)
			xII[i-n1]=x[i];		
		for (int i = n1+n2; i < dim; i++)
			xIII[i-n1-n2]=x[i];		
			
		return BaseFunctions.getSphere(xI)+BaseFunctions.getAckley(xII)
				+BaseFunctions.getRastrigin(xIII)+Fstar;
	}		
	
	/*HF5: Hybrid Function 5*/
	public static double getHF5(double x[], double[] shiftValues, double[][] rotationMatrix) {
		
		double Fstar=0.0;
		int dim = x.length;
		
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);		
		
		double[] p = {0.4,0.6};
		int n1 = (int)Math.ceil(p[0]*dim);			
		int n2 = dim-n1;
		
		double[] xI = new double[n1];
		double[] xII = new double[n2];
		
		for (int i = 0; i < n1; i++)
			xI[i]=x[i];		
		for (int i = n1; i < dim; i++)
			xII[i-n1]=x[i];		
			
		return BaseFunctions.getKatsuura(xI)+BaseFunctions.getHappyCat(xII)+Fstar;
	}	
	
	/*HF6: Hybrid Function 6*/
	public static double getHF6(double x[], double[] shiftValues, double[][] rotationMatrix) {
		
		double Fstar=0.0;
		int dim = x.length;
		
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);		
		
		double[] p = {0.4,0.6};
		int n1 = (int)Math.ceil(p[0]*dim);			
		int n2 = dim-n1;
		
		double[] xI = new double[n1];
		double[] xII = new double[n2];
		
		for (int i = 0; i < n1; i++)
			xI[i]=x[i];		
		for (int i = n1; i < dim; i++)
			xII[i-n1]=x[i];		
			
		return 1.0 + BaseFunctions.getSphere(xI)+BaseFunctions.getRosenbrock(xII)+Fstar;
	}	
	
	/*HF7: Hybrid Function 7*/
	public static double getHF7(double x[], double[] shiftValues, double[][] rotationMatrix) {
		
		double Fstar=0.0;
		int dim = x.length;
		
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);		
		
		double[] p = {0.4,0.6};
		int n1 = (int)Math.ceil(p[0]*dim);			
		int n2 = dim-n1;
		
		double[] xI = new double[n1];
		double[] xII = new double[n2];
		
		for (int i = 0; i < n1; i++)
			xI[i]=x[i];		
		for (int i = n1; i < dim; i++)
			xII[i-n1]=x[i];		
			
		return 1.0 + BaseFunctions.getMean(xI)+BaseFunctions.getRastrigin(xII)+Fstar;
	}	
	
	/*HF8: Hybrid Function 8*/
	public static double getHF8(double x[], double[] shiftValues, double[][] rotationMatrix) {
		
		double Fstar=0.0;
		int dim = x.length;
		
		shiftValues_ = shiftValues;
		rotationMatrix_ = rotationMatrix;
		shiftVariables(x);
		x=rotateVariables(x);		
		
		double[] p = {0.4,0.6};
		int n1 = (int)Math.ceil(p[0]*dim);			
		int n2 = dim-n1;
		
		double[] xI = new double[n1];
		double[] xII = new double[n2];
		
		for (int i = 0; i < n1; i++)
			xI[i]=x[i];		
		for (int i = n1; i < dim; i++)
			xII[i-n1]=x[i];		
			
		return 1.0 + BaseFunctions.getSphere(xI) + BaseFunctions.getMean(xII)+Fstar;
	}	

	/*LF1: Linkage Function 1*/
	public static double[] getLF1(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		double linkages;
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				linkages = (0.3*xI[0]*xI[0]*Math.cos(24*Math.PI*xI[0] + 4*(index[i][j]*Math.PI)/D) + 0.6*xI[0])
						*Math.sin(6*Math.PI*xI[0] + (index[i][j]*Math.PI)/D);
				g[i] += Math.pow(xIII[i][j]-linkages, 2);
			}
			g[i] = 2*g[i]/xIII[i].length;
		}
		return g;
	}
	
	/*LF2: Linkage Function 2*/
	public static double[] getLF2(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		double linkage;
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				linkage = xIII[i][j] - Math.pow(xI[0],0.5*(1.0+(3.0*(index[i][j]-2)/(D-2))));
				g[i] += Math.pow(linkage, 2);
			}
			g[i] = 2*g[i]/xIII[i].length;
		}
		return g;
	}
	
	/*LF3: Linkage Function 3*/
	public static double[] getLF3(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		double linkage;
		double part1, part2;
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			part1 = 0.0;
			part2 = 1.0;
			for(int j=0; j<xIII[i].length; j++){
				linkage = xIII[i][j] - Math.pow(xI[0],0.5*(1.0+(3.0*(index[i][j]-2)/(D-2))));
				part1 += linkage*linkage;
				part2 *= Math.cos(20*Math.PI*linkage/Math.sqrt(index[i][j]));
			}
			g[i] = 4*part1 - 2*part2 + 2;
			g[i] = 2*g[i]/xIII[i].length;
		}
		return g;
	}
	
	/*LF4-5: Linkage Function 4 and 5*/
	public static double[] getLF4_5(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		double linkages;
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				if(i % 2 == 0){
					linkages = 0.8*xI[0]*Math.cos(6*Math.PI*xI[0] + (index[i][j]*Math.PI)/D);
				}else {
					linkages = 0.8*xI[0]*Math.sin(6*Math.PI*xI[0] + (index[i][j]*Math.PI)/D);
				}
				g[i] += Math.pow(xIII[i][j]-linkages, 2);
			}
			g[i] = 2*g[i]/xIII[i].length;
		}
		return g;
	}
	
	/*LF6: Linkage Function 6*/
	public static double[] getLF6(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		double linkages;
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				linkages = 2.0*xI[1]*Math.sin(2*Math.PI*xI[0] + (index[i][j]*Math.PI)/D);
				g[i] += Math.pow(xIII[i][j]-linkages, 2);
			}
			g[i] = 2*g[i]/xIII[i].length;
		}
		return g;
	}
	
	/*LF7: Linkage Function 6*/
	public static double[] getLF7(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		double linkages;
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				linkages = Math.sin(6*Math.PI*xI[0] + (index[i][j]*Math.PI)/D);
				g[i] += Math.pow(xIII[i][j]-linkages, 2);
			}
			g[i] = 2*g[i]/xIII[i].length;
		}
		return g;
	}
	
	/*LF8: Linkage Function 8*/
	public static double[] getLF8(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		double linkages;
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				linkages = 0.8*xI[0]*Math.cos(6*Math.PI*xI[0] + (index[i][j]*Math.PI)/D);
				g[i] += Math.pow(xIII[i][j]-linkages, 2);
			}
			g[i] = 2*g[i]/xIII[i].length;
		}
		return g;
	}
	
	/*LF9: Linkage Function 9*/
	public static double[] getLF9(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		double linkages;
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				linkages = 0.8*xI[0]*Math.sin(6*Math.PI*xI[0] + (index[i][j]*Math.PI)/D);
				g[i] += Math.pow(xIII[i][j]-linkages, 2);
			}
			g[i] = 2*g[i]/xIII[i].length;
		}
		return g;
	}
	
	
	/*LF10: Linkage Function 10*/
	public static double[] getLF10(double[] xI, double xIII[][], int[][] index, int D) {
		int numObj = xIII.length;
		double[] g = new double[numObj];
		double linkage;
		for(int i=0; i<numObj; i++){
			g[i] = 0.0;
			for(int j=0; j<xIII[i].length; j++){
				linkage = xIII[i][j] - Math.pow(xI[0],0.5*(1.0+(3.0*(index[i][j]-2)/(D-2))));
				g[i] += 4*Math.pow(linkage, 2) - Math.cos(8*Math.PI*linkage) + 1.0;
			}
			g[i] = 2*g[i]/xIII[i].length;
		}
		return g;
	}
	
	/*DF1: Deep grouping-based Function 1*/
	public static double[] getDF1(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getSphere(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getSphere(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF2: Deep grouping-based Function 2*/
	public static double[] getDF2(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getMean(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getMean(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF3: Deep grouping-based Function 3*/
	public static double[] getDF3(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getSphere(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRosenbrock(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF4: Deep grouping-based Function 4*/
	public static double[] getDF4(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getSphere(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRastrigin(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF5: Deep grouping-based Function 5*/
	public static double[] getDF5(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getSphere(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getMean(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF6: Deep grouping-based Function 6*/
	public static double[] getDF6(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%3 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getSphere(xII[i][j]);
				}
			}else if(i%3 == 1){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRosenbrock(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getAckley(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF7: Deep grouping-based Function 7*/
	public static double[] getDF7(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%3 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getMean(xII[i][j]);
				}
			}else if(i%3 == 1){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRosenbrock(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getGriewank(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF8: Deep grouping-based Function 8*/
	public static double[] getDF8(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%3 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getSphere(xII[i][j]);
				}
			}else if(i%3 == 1){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRastrigin(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getAckley(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF9: Deep grouping-based Function 9*/
	public static double[] getDF9(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%3 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRosenbrock(xII[i][j]);
				}
			}else if(i%3 == 1){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getAckley(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRastrigin(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF10: Deep grouping-based Function 10*/
	public static double[] getDF10(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getSphere(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getAckley(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF11: Deep grouping-based Function 11*/
	public static double[] getDF11(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getSphere(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getGriewank(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF12: Deep grouping-based Function 12*/
	public static double[] getDF12(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRosenbrock(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getGriewank(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF13: Deep grouping-based Function 13*/
	public static double[] getDF13(double x1, double xII[][][],int[][][] index, int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			if(i%2 == 0){
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRosenbrock(xII[i][j]);
				}
			}else{
				for(int j=0;j<xII[i].length;j++){
					g[i] += BaseFunctions.getRastrigin(xII[i][j]);
				}
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF14: Deep grouping-based Function 14*/
	public static double[] getDF14(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getRosenbrock(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF15: Deep grouping-based Function 15*/
	public static double[] getDF15(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getAckley(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF16: Deep grouping-based Function 16*/
	public static double[] getDF16(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getWeierstrass(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF17: Deep grouping-based Function 17*/
	public static double[] getDF17(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getGriewank(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF18: Deep grouping-based Function 18*/
	public static double[] getDF18(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getRastrigin(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF19: Deep grouping-based Function 19*/
	public static double[] getDF19(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getElliptic(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF20: Deep grouping-based Function 20*/
	public static double[] getDF20(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getCigar(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF21: Deep grouping-based Function 21*/
	public static double[] getDF21(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getDiscus(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF22: Deep grouping-based Function 22*/
	public static double[] getDF22(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getMSchwefel(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF23: Deep grouping-based Function 23*/
	public static double[] getDF23(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getKatsuura(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF24: Deep grouping-based Function 24*/
	public static double[] getDF24(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getHappyCat(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	/*DF25: Deep grouping-based Function 25*/
	public static double[] getDF25(double x1, double xII[][][], int[][][] index,int dim, String linkageType) {
		if(linkageType.equalsIgnoreCase("linear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + (index[i][j][k]+1.0)/(dim))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}else if(linkageType.equalsIgnoreCase("nonlinear")){
			for(int i=0;i<xII.length;i++){
				for(int j=0;j<xII[i].length;j++){
					for(int k=0;k<xII[i][j].length;k++){
						xII[i][j][k] = (1.0 + Math.cos(0.5*Math.PI*((index[i][j][k]+1.0)/(dim))))*xII[i][j][k] - 10.0*x1;
					}
				}
			}
		}	
		int numObj = xII.length;
		double[] g = new double[numObj];
		for (int i = 0; i < numObj; i++) {
			g[i] = 0.0;
			for(int j=0;j<xII[i].length;j++){
				g[i] += BaseFunctions.getExGriewRosen(xII[i][j]);
			}
			g[i] = g[i]/xII[i].length;
		}
		
		return g;
	}
	
	
	protected static void shiftVariables(double x[]) {
		for (int i = 0; i < x.length; i++) {
			x[i] -= shiftValues_[i];
		}
	}

	protected static void shrinkVariables(double x[], double mul) {
		for (int i = 0; i < x.length; i++) 
			x[i] *= mul;
	}	
	
	protected static double[] rotateVariables(double x[]) {
		int len = x.length;
		double res[] = new double[len];

		for (int i = 0; i < len; i++) {
			double[] y = rotationMatrix_[i];

			double sum = 0;
			for (int j = 0; j < len; j++)
				sum += x[j] * y[j];
			res[i] = sum;
//			res[i] = x[i];
		}

		return res;
	}	
	
}
