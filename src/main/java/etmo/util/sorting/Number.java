package etmo.util.sorting;

public class Number implements Comparable<Number>{
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
