package etmo.test;
import etmo.test.Test_ojb;

public class Test_2 {
    public static void main(String[] args) throws InterruptedException {
        while (true){
            Test_ojb o = new Test_ojb();
            System.out.println(o.test1());
            Thread.sleep(1000);
        }
    }
}
