package etmo.problems.base.dynamicBase;

public class DGFunctions {
	
	/*Dynamic Sphere Function without linkages used in DF1, DF2, DF5*/
	public static double getSphere(double x[], double gt) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++)
			sum += (x[i]-gt)*(x[i]-gt);		
		return sum;
	}
	/*Dynamic Sphere Function with linkages used in DF3*/
	public static double getSphereLinkage1(double x[], double gt, double f1, double ht) {
		double sum = 0.0;
		double linkage = gt + Math.pow(f1, ht);
		for (int i = 0; i < x.length; i++)
			sum += (x[i] - linkage)*(x[i] - linkage);		
		return sum;
	}
	
	/*Dynamic Sphere Function with linkages used in DF10*/
	public static double getSphereLinkage2(double x[], double[] pv, double gt) {
		double sum = 0.0;
		double sumPV = 0.0;
		
		for(int i = 0; i < pv.length; i++){
			sumPV += pv[i];
		}
		double linkage = Math.sin(2*Math.PI*sumPV)/(1.0+Math.abs(gt));
		for (int i = 0; i < x.length; i++)
			sum += (x[i]-linkage)*(x[i]-linkage);		
		return sum;
	}
	
	/*Dynamic Sphere Function with linkages used in DF11*/
	public static double getSphereLinkage3(double x[], double[] pv, double gt) {
		double sum = 0.0;
		double linkage = 0.5*gt*pv[0];
		for (int i = 0; i < x.length; i++)
			sum += (x[i]-linkage)*(x[i]-linkage);		
		return sum + gt;
	}
	
	/*Dynamic Sphere Function with linkages used in DF12*/
	public static double getSphereLinkage4(double x[], double[] pv, double t) {
		double sum = 0.0;
		double linkage = Math.sin(t*pv[0]);
		double prod = 1.0;
		for (int i = 0; i < x.length; i++)
			sum += (x[i]-linkage)*(x[i]-linkage);
		for (int i = 0; i < pv.length; i++)
			prod *= Math.sin(Math.floor(10*Math.sin(Math.PI*t)*(2*pv[i]-1))*Math.PI/2);
		return sum + Math.abs(prod);
	}
	
	/*Dynamic Sphere Function with linkages used in DF8*/
	public static double getSphereLinkage5(double x[], double gt, double x1) {
		double sum = 0.0;
		double bt = 100*gt*gt;
		double linkage = gt*Math.sin(4*Math.PI*Math.pow(x1, bt))/(1+Math.abs(gt));
		for (int i = 0; i < x.length; i++)
			sum += (x[i]-linkage)*(x[i]-linkage);
		return sum;
	}
	
	/*Dynamic Rastrigin Function without linkages used in DF6*/
	public static double getRastrigin(double x[], double gt) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++){
			double y = x[i] - gt;
			sum += Math.abs(gt)*y*y - 10*Math.cos(2*Math.PI*y) + 10;	
		}		
		return sum;
	}
	
	/*Dynamic Rosenbrock Function without linkages*/
	public static double getRosenbrock(double x[], double gt) {
		double sum = 0.0;
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i++){
			y[i] = x[i] - gt;	
		}	
		for (int i = 0; i < x.length - 1; i++) {
			double t = 100 * (y[i] * y[i] - y[i + 1]) * (y[i] * y[i] - y[i + 1]) + (1 - y[i]) * (1 - y[i]);
			sum += t;
		}

		return sum;
	}
	
	/*Dynamic Ackley Function without linkages*/
	public static double getAckley(double x[], double gt) {
		double sum1 = 0;
		double sum2 = 0;
		for (int i = 0; i < x.length; i++) {
			double y = x[i] - gt;
			sum1 += ((y * y) / x.length);
			sum2 += (Math.cos(2 * Math.PI * y) / x.length);
		}
		return -20 * Math.exp(-0.2 * Math.sqrt(sum1)) - Math.exp(sum2) + 20 + Math.E;
	}
	
	/*Dynamic Griewank Function without linkages*/
	public static double getGriewank(double x[], double gt) {
		int k = 1;
		double sum = 0;
		double prod = 1; 
		for (int i = 0; i < x.length; i++) {
			double y = x[i] - gt;
			sum += (y * y);
			prod *= (k * Math.cos(y / Math.sqrt(i + 1)));
		}
		return k + sum / 4000 - prod;
	}
}
