package etmo.metaheuristics.matmy2;

import jdk.jshell.execution.Util;

public class test {
    public static void main(String[] args) {
        double[] vector1 = new double[]{1, 0};
        double[] vector2 = new double[]{1, 1};
        double[] vector3 = new double[]{0, 1};
        double[] vector4 = new double[]{-1, 0};

        // 0
        System.out.println(Utils.calVectorAngle(vector1, vector1));
        // 0.78
        System.out.println(Utils.calVectorAngle(vector1, vector2));
        // 1.57
        System.out.println(Utils.calVectorAngle(vector1, vector3));
        // 0
        System.out.println(Utils.calVectorAngle(vector1, vector4));
    }
}
