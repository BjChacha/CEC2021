package etmo.util.math;

import etmo.util.PseudoRandom;

public class Random {
    public static int rouletteWheel(double[] arr){
		double[] sm = NonLinear.softmax(arr);

		double s = 0;
		double p = PseudoRandom.randDouble();
		int idx;

		for (idx = 0; idx < sm.length - 1; idx++) {
			s += sm[idx];
			if (s >= p) break;
		}
		return idx;
	}

	public static int rouletteWheel(double[] arr, int masked){
		double[] newArr = new double[arr.length - 1];
		for (int i = 0; i < newArr.length; i ++) {
			if (i < masked) {
				newArr[i] = arr[i];
			} else {
				newArr[i] = arr[i+1];
			}
		} 
		
		double[] sm = NonLinear.softmax(newArr);

		double s = 0;
		double p = PseudoRandom.randDouble();
		int idx;

		for (idx = 0; idx < sm.length - 1; idx++) {
			s += sm[idx];
			if (s >= p) break;
		}

		idx = idx < masked ? idx : idx + 1; 
		return idx;
	}
}
