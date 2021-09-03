package etmo.util.sorting;

import java.util.Arrays;

public class SortingIdx {
    static public int[] sort(double[] a, boolean inversed){
        int[] resIdx = new int[a.length];
        int sign = inversed ? -1 : 1;
        Number sorted[] = new Number[a.length];
        for (int i = 0; i < a.length; ++i) {
            sorted[i] = new Number(sign * a[i], i);
        }
        Arrays.sort(sorted);

        for (int i = 0; i < sorted.length; i++){
            resIdx[i] = sorted[i].index;
        }

        return resIdx;
    }
}