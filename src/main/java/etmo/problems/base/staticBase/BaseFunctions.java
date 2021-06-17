package etmo.problems.base.staticBase;

public class BaseFunctions {
	
	/*f-1: Mean Function*/
	public static double getMean(double x[]) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++)
			sum += 9*Math.abs(x[i]);		
		return sum/x.length;
	}	
	
	/*f0: Sphere Function*/
	public static double getSphere(double x[]) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++)
			sum += x[i]*x[i];		
		return sum;
	}	
	
	/*f1: High Conditioned Elliptic Function*/
	public static double getElliptic(double x[]) {
		double sum = 0;
		int dim = x.length;
		double a = Math.pow(10.0, 6);
		
		for (int i = 0; i < dim; i++)
			sum += Math.pow(a, i/(dim-1.0))*(x[i] * x[i]);
		return sum;
	}	
	
	/*f2: Bent Cigar Function*/
	public static double getCigar(double x[]) {
		
		double sum = 0;
		for (int i = 1; i < x.length; i++)
			sum += (x[i] * x[i]);
		
		return Math.pow(x[0], 2)+Math.pow(10, 6)*sum;
	}	
	
	/*f3: Discus Function*/
	public static double getDiscus(double x[]) {
		double sum = 0;
		for (int i = 1; i < x.length; i++)
			sum += (x[i] * x[i]);
		
		return Math.pow(10, 6)*Math.pow(x[0], 2)+sum;
	}
	
	/*f4: Rosenbrock Function*/
	public static double getRosenbrock(double x[]) {
		double sum = 0;
		
		//x[0] += 1.0;//shift to origin
		for (int i = 0; i < x.length - 1; i++) {
			
			//x[i+1] += 1.0;//shift to origin
			
			double t = 100 * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1]) + (1 - x[i]) * (1 - x[i]);
			sum += t;
		}

		return sum;
	}

	/*f5: Ackley Function*/
	public static double getAckley(double x[]) {
		double sum1 = 0;
		double sum2 = 0;

		for (int i = 0; i < x.length; i++) {
			sum1 += ((x[i] * x[i]) / x.length);
			sum2 += (Math.cos(2 * Math.PI * x[i]) / x.length);
		}

		return -20 * Math.exp(-0.2 * Math.sqrt(sum1)) - Math.exp(sum2) + 20 + Math.E;

	}	
	
	/*f6: Weierstrass Function*/
	public static double getWeierstrass(double x[]) {

		double result= 0.0, part1= 0.0, part2 = 0.0;
		double a = 0.5;
		double b = 3.0;	
		double k = 20.0;
		double w = 2 * Math.PI;
		

		for (int i = 0; i < x.length; i++) {
			for(int j = 0; j <= k; j++)
				part1 += Math.pow(a, j)* Math.cos(w*Math.pow(b, j)*(x[i]+0.5));
		}
		
		for(int j = 0; j <= k; j++)
			part2 +=  Math.pow(a, j)* Math.cos(w*Math.pow(b, j)*0.5);		
		
		result=part1-x.length*part2;

		return result;
	}	
	
	/*f7: Griewank Function*/
	public static double getGriewank(double x[]) {
		int k = 1;

		double sum = 0;
		double prod = 1; 

		for (int i = 0; i < x.length; i++) {
			sum += (x[i] * x[i]);
			prod *= (k * Math.cos(x[i] / Math.sqrt(i + 1)));
		}

		return k + sum / 4000 - prod;

	}
	public static double getGriewank1D(double x) {
		
		return Math.pow(x, 2)/4000-Math.cos(x)+1;

	}
	
	/*f8: Rastrigin Function*/
	public static double getRastrigin(double x[]) {

		double result = 0.0;
		double a = 10.0;
		double w = 2 * Math.PI;

		for (int i = 0; i < x.length; i++) {
			result += x[i] * x[i] - a * Math.cos(w * x[i]);
		}
		result += a * x.length;

		return result;
	}
	
	/*f9: Modified Schwefel Function*/
	public static double getMSchwefel(double x[]) {
		double sum = 0;
		int dim = x.length;
		double[] z = new double[dim];
		double[] gz = new double[dim];
		for (int i = 0; i < dim; i++) {
			z[i]=x[i]+4.209687462275036e+002;
			if(Math.abs(z[i])<=500)
				gz[i]=z[i]*Math.sin(Math.pow(Math.abs(z[i]), 0.5));
			else if(z[i]>500)
				gz[i]=(500-z[i]%500)*Math.sin(Math.pow(Math.abs(500-z[i]%500), 0.5))   -Math.pow(z[i]-500,2)/(10000*dim);
			else
				gz[i]=(Math.abs(z[i])%500-500)*Math.sin(Math.pow(500-Math.abs(z[i])%500, 0.5))   -Math.pow(z[i]+500,2)/(10000*dim);
			
			sum += gz[i];
		}
		
		sum=418.9829*dim-sum;
		
		return sum;
	}	

	/*f10: Katsuura Function*/
	public static double getKatsuura(double x[]) {
		int index=32;
		int dim = x.length;
		
		double sum = 0;
		double prod = 1;
		for (int i = 0; i < dim; i++) {
			for (int j = 1; j <= index; j++)
				sum += Math.abs(Math.pow(2, j)*x[i]-Math.round(Math.pow(2, j)*x[i]))/Math.pow(2, j);
			prod *= Math.pow(1+i*sum , 10/Math.pow(dim, 1.2));
		}
		
		prod = 10/Math.pow(dim, 1.2)*prod-10/Math.pow(dim, 1.2);
		
		return prod;
	}
	
	/*f11: HappyCat Function*/
	public static double getHappyCat(double x[]) {
		double sum, sum1 = 0, sum2=0;
		int dim = x.length;
		
		for (int i = 0; i < dim; i++) {
			x[i] -= 1.0;//shift to origin
			sum1 += (x[i] * x[i]);
			sum2 += x[i];
		}
		
		sum = Math.pow(Math.abs(sum1-dim), 0.25)+(0.5*sum1+sum2)/dim+0.5;
		return sum;
	}
	
	/*f12: HGBat Function*/
	public static double getHGBat(double x[]) {
		
		double sum, sum1 = 0, sum2=0;
		int dim = x.length;
		
		for (int i = 0; i < dim; i++) {
			x[i] -= 1.0;//shift to origin
			sum1 += (x[i] * x[i]);
			sum2 += x[i];
		}
		
		sum = Math.pow(Math.abs(Math.pow(sum1, 2)-Math.pow(sum2, 2)), 0.5)+(0.5*sum1+sum2)/dim+0.5;
		
		return sum;
	}	
	
	/*f13: Expanded Griewank plus Rosenbrock Function*/
	public static double getExGriewRosen(double x[]) {
		double sum = 0;
		int dim = x.length;
		double tmp1,tmp2,temp;
		x[0] += 1.0;//shift to origin
		
	    for (int i=0; i<dim-1; i++)
	    {
			x[i+1] += 1.0;//shift to origin
			tmp1 = x[i]*x[i]-x[i+1];
			tmp2 = x[i]-1.0;
	        temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
	         sum += (temp*temp)/4000.0 - Math.cos(temp) + 1.0;
	    }
	    
		tmp1 = x[dim-1]*x[dim-1]-x[0];
		tmp2 = x[dim-1]-1.0;
	    temp = 100.0*tmp1*tmp1 + tmp2*tmp2;;
	    		
	    sum += (temp*temp)/4000.0 - Math.cos(temp) + 1.0 ;
		
		return sum;
	}	
	
	/*f14: Expanded Scaffer F6 Function*/
	public static double getScafferF6(double x[]) {
		double sum = 0;
		for (int i = 0; i < x.length; i++)
			sum += ScafferBase(x[i], x[(i+1)% x.length]);
		return sum;
	}
		
	private static double ScafferBase(double x, double y) {
		double result = 0.5;		
		double pSum=Math.pow(x, 2)+Math.pow(y, 2);
		
		result +=(Math.pow(Math.sin(Math.pow(pSum, 0.5)),2)-0.5)/Math.pow(1+0.001*(pSum),2);
				
		return result;
	}
		
}
