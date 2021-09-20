package etmo.util.math;

public class NonLinear {
    public static double[] softmax(double[] arr){
		double[] res = new double[arr.length];
		double sum = 0;

		for (int i = 0; i < res.length; i++){
			res[i] = Math.exp(arr[i]);
			sum += res[i];
		}

		for (int i = 0; i < res.length; i++) {
			res[i] /= sum;
		}
		return res;
	}
}
