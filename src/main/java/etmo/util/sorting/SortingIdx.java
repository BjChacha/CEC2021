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

class Number implements Comparable<Number>{
    Double data;
    int index;

    Number(double d, int i){
        this.data = d;
        this.index = i;
    }

    @Override
    public int compareTo(Number o){
        return this.data.compareTo(o.data);
    }

    public double getData(){
        return this.data;
    }

    public int getIndex(){
        return this.index;
    }
}