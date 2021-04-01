package etmo.metaheuristics.matmy2;

import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.Ranking;
import org.apache.commons.math3.analysis.function.Inverse;
import org.datavec.api.transform.analysis.columns.NDArrayAnalysis;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.api.ops.impl.transforms.custom.MatrixInverse;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.eigen.Eigen;
import static org.nd4j.linalg.factory.Nd4j.linspace;
import org.nd4j.linalg.schedule.InverseSchedule;
import org.nd4j.linalg.inverse.InvertMatrix;

import java.util.Arrays;

public class LTR {
    public static double[][][] LTR(SolutionSet[] archives, int D) throws JMException {
        // assume population are evaluated
        classify(archives);
        double[][] Ls = calLs(archives);
        double[][] Ld = calLd(archives);
        double[][] Lk = calLk(archives);
        double[][] Z = getZ(archives);
        double[][] M = calM(Ls, Ld, Lk, Z, D);
        double[][][] Ms = sliceM(M, archives.length);
        return Ms;
    }

    private static double[][] calLs(SolutionSet[] archives) throws JMException {
        int taskNum = archives.length;
        // 1. Get Ws
        int SN = 0;
        for (int k = 0; k < taskNum; k++){
            SN += archives[k].size();
        }
        double[][] Ws = new double[SN][SN];

        int rowAccu = 0;
        for (int k1 = 0; k1 < taskNum; k1++){
            int colAccu = 0;
            for (int k2 = 0; k2 < taskNum; k2++){
                for (int i = 0; i < archives[k1].size(); i++){
                    for (int j = 0; j < archives[k2].size(); j++){
                        Ws[rowAccu+i][colAccu+j] = isSimilar(archives[k1].get(i), archives[k2].get(j));
                    }
                }
                colAccu += archives[k2].size();
            }
            rowAccu += archives[k1].size();
        }

        // 2. Get Ds
        double[][] Ds = new double[SN][SN];
        for (int i = 0; i < Ds.length; i++){
            Ds[i][i] = Arrays.stream(Ws[i]).sum();
        }

        // 3. Get Ls
        double[][] Ls = new double[SN][SN];
        for (int i = 0; i < Ls.length; i++){
            int finalI = i;
            Arrays.setAll(Ls[i], j -> Ds[finalI][j] - Ws[finalI][j]);
        }

        return Ls;
    }


    private static double[][] calLd(SolutionSet[] archives) throws JMException {
        int taskNum = archives.length;
        // 1. Get Wd
        int SN = 0;
        for (int k = 0; k < taskNum; k++){
            SN += archives[k].size();
        }
        double[][] Wd = new double[SN][SN];

        int rowAccu = 0;
        for (int k1 = 0; k1 < taskNum; k1++){
            int colAccu = 0;
            for (int k2 = 0; k2 < taskNum; k2++){
                for (int i = 0; i < archives[k1].size(); i++){
                    for (int j = 0; j < archives[k2].size(); j++){
                        Wd[rowAccu+i][colAccu+j] = 1 - isSimilar(archives[k1].get(i), archives[k2].get(j));
                    }
                }
                colAccu += archives[k2].size();
            }
            rowAccu += archives[k1].size();
        }

        // 2. Get Dd
        double[][] Dd = new double[SN][SN];
        for (int i = 0; i < Dd.length; i++){
            Dd[i][i] = Arrays.stream(Wd[i]).sum();
        }

        // 3. Get Ld
        double[][] Ld = new double[SN][SN];
        for (int i = 0; i < Ld.length; i++){
            int finalI = i;
            Arrays.setAll(Ld[i], j -> Dd[finalI][j] - Wd[finalI][j]);
        }

        return Ld;
    }


    private static double[][] calLk(SolutionSet[] archives) throws JMException {
        int taskNum = archives.length;
        // 1. Get Wk
        int SN = 0;
        for (int k = 0; k < taskNum; k++){
            SN += archives[k].size();
        }
        double[][] Wk = new double[SN][SN];

        int accu = 0;
        for (int k = 0; k < taskNum; k++){
            for (int i = 0; i < archives[k].size(); i++){
                for (int j = 0; j < archives[k].size(); j++){
                    Wk[accu+i][accu+j] = getSimilarity(archives[k].get(i), archives[k].get(j));
                }
            }
            accu += archives[k].size();
        }

        // 2. Get Dk
        double[][] Dk = new double[SN][SN];
        for (int i = 0; i < Dk.length; i++){
            Dk[i][i] = Arrays.stream(Wk[i]).sum();
        }

        // 3. Get Lk
        double[][] Lk = new double[SN][SN];
        for (int i = 0; i < Lk.length; i++){
            int finalI = i;
            Arrays.setAll(Lk[i], j -> Dk[finalI][j] - Wk[finalI][j]);
        }

        return Lk;
    }

    private static double[][] getZ(SolutionSet[] archives) throws JMException {
        int taskNum = archives.length;
        // Get Z
        int rows = 0;
        int cols = 0;
        for (int k = 0; k < taskNum; k++){
            rows += archives[k].size();
            cols += archives[k].get(0).getDecisionVariables().length;
        }

        double[][] Z = new double[rows][cols];
        int rowAccu = 0;
        int colAccu = 0;
        for (int k = 0; k < taskNum; k++){
            int iLen = archives[k].size();
            int jLen = archives[k].get(0).getDecisionVariables().length;
            for (int i = 0; i < iLen; i++){
                for (int j = 0; j < jLen; j++){
                    Z[rowAccu+i][colAccu+j] = archives[k].get(i).getDecisionVariables(j);
                }
            }
            rowAccu += iLen;
            colAccu += jLen;
        }
        return Z;
    }

    private static double[][] calM(double[][] ls, double[][] ld, double[][] lk, double[][] z, int D) {
        NDArray Ls = new NDArray(ls);
        NDArray Ld = new NDArray(ld);
        NDArray Lk = new NDArray(lk);
        NDArray Z = new NDArray(z);

        // Left part
        INDArray tmp1 = Lk.add(Ls);
        INDArray left = Z.transpose().mmul(tmp1);
        left = left.mmul(Z);

        // Right part
        INDArray right = Z.transpose().mmul(Ld);
        right = right.mmul(Z);

        INDArray res = left;

        // 特征分解
        INDArray eigens = Eigen.symmetricGeneralizedEigenvalues(res, right, true);
        int[] cols = linspace(0, D-1, D).toIntVector();
        INDArray r = res.getColumns(cols);
        return r.toDoubleMatrix();
    }


    private static double[][][] sliceM(double[][] M, int K){
        int l = M.length / K;
        double[][][] Ms = new double[K][][];
        NDArray NM = new NDArray(M);
        for (int i = 0; i < K; i++){
            int[] rows = linspace(i*l, (i+1)*l-1, l).toIntVector();
            Ms[i] = NM.getRows(rows).toDoubleMatrix();
        }
        return Ms;
    }


    private static int isSimilar(Solution s1, Solution s2) throws JMException {
        int res = s1.getLabel() == s2.getLabel() ? 1 : 0;
        return res;
    }


    private static double getSimilarity(Solution s1, Solution s2) throws JMException {
        double e = 0;
        for (int i = 0; i < s1.numberOfVariables(); i++){
            e += Math.pow((s1.getDecisionVariables(i) - s2.getDecisionVariables(i)), 2);
        }
        e = Math.sqrt(e);
        double res = Math.exp(e);
        return res;
    }


    public static void classify(SolutionSet[] archives){
        for (int k = 0; k < archives.length; k++) {
            Distance distance = new Distance();
            Ranking ranking = new Ranking(archives[k]);
            int remain = archives[k].size();
            int index = 0;
            SolutionSet front = null;

            int nObjs = archives[k].get(0).getProblemSet().get(0).getNumberOfObjectives();

            archives[k].clear();

            int tied1Cnt = 0;
            int tied3Cnt = 0;

            while (remain > 0) {
                front = ranking.getSubfront(index);
                if (index == 0){
                    tied1Cnt = front.size();
                    tied3Cnt = remain - tied1Cnt;
                    tied1Cnt = tied1Cnt / 2;
                    tied3Cnt = tied3Cnt / 2;
                }
                distance.crowdingDistanceAssignment(front, nObjs);
                for (int i = 0; i < front.size(); i++) {
                    if (index == 0){
                        if (tied1Cnt > 0){
                            front.get(i).setLabel(1);
                            tied1Cnt --;
                        }else
                            front.get(i).setLabel(2);
                    }else{
                        if (tied3Cnt > 0){
                            front.get(i).setLabel(3);
                            tied3Cnt --;
                        }else
                            front.get(i).setLabel(4);
                    }
                    archives[k].add(front.get(i));
                }
                remain -= front.size();
                index++;
            }
        }
    }


    public static Solution Mapping(Solution individual, double[][] M1, double[][] M2) throws JMException {
        double[][] variables = new double[1][];
        variables[0] = individual.getDecisionVariablesInDouble();
        
        NDArray Nv1 = new NDArray(variables);
        NDArray NM1 = new NDArray(M1);
        NDArray NM2 = new NDArray(M2);
        NDArray Nv2 = (NDArray) (Nv1.mmul(NM1)).mmul(InvertMatrix.pinvert(NM2, false));

        double[] newVariables = Nv2.toDoubleMatrix()[0];

        for (int i = 0; i < newVariables.length; i++){
            if (newVariables[i] > 1.0 || newVariables[i] < 0.0)
                System.out.println("Error! Border");
            if (Double.isNaN(newVariables[i]))
                System.out.println("Error NaN");
        }

        Solution newIndividual = new Solution(individual);
        newIndividual.setDecisionVariables(newVariables);
        return newIndividual;
    }


    public static void test(double[][] input){
        NDArray a = new NDArray(input);
        INDArray es = Eigen.symmetricGeneralizedEigenvalues(a);
        System.out.println("Demo");
    }
}