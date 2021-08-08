package etmo.test;

public class Test_WeightVector {
    public static void main(String[] args){
        int H_= 0;
        int populationSize_ = 100;
        int objectNumber = 3;
        double[][] lambda_ = new double[populationSize_][objectNumber];

        while ((H_+2) * (H_+1) / 2 < populationSize_)
            H_++;

        int nw = 0;
        for (int i = 0; i <= H_; i++) {
            for (int j = 0; j <= H_ - i; j++) {
                lambda_[nw][0] = (1.0 * i) / H_;
                lambda_[nw][1] = (1.0 * j) / H_;
                lambda_[nw][2] = 1.0 * (H_ - i - j) / H_;
                nw++;
            } // for
        }
        System.out.println("Break point");
    }
}
