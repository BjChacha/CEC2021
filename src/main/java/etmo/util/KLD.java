package etmo.util;

import etmo.core.ProblemSet;
import etmo.core.SolutionSet;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.correlation.Covariance;

public class KLD {
    private int maxVarNum;
    private int varNum;
    private final int taskNum;

    private double cov[][][];
    private double covInv[][][];
    private double covDet[];
    private double kld[];

    private SolutionSet[] archives;
    private ProblemSet problemSet;

    public KLD(ProblemSet problemSet, SolutionSet[] solutionSet){
        archives = solutionSet;
        this.problemSet = problemSet;

        varNum = Integer.MAX_VALUE;
        maxVarNum = 0;
        for (int i = 0; i < archives.length; i++){
            varNum = Math.min(varNum, archives[i].get(0).numberOfVariables());
            maxVarNum = Math.max(maxVarNum, archives[i].get(0).numberOfVariables());
        }
        taskNum = problemSet.size();
        cov = new double[taskNum][maxVarNum][maxVarNum];
        covInv = new double[taskNum][maxVarNum][maxVarNum];
        covDet = new double[taskNum];
        kld = new double[taskNum];
    }

    public KLD(int taskNum, SolutionSet[] solutionSet){
        archives = solutionSet;

        varNum = Integer.MAX_VALUE;
        maxVarNum = 0;
        for (int i = 0; i < archives.length; i++){
            varNum = Math.min(varNum, archives[i].get(0).numberOfVariables());
            maxVarNum = Math.max(maxVarNum, archives[i].get(0).numberOfVariables());
        }
        this.taskNum = taskNum;
        cov = new double[taskNum][maxVarNum][maxVarNum];
        covInv = new double[taskNum][maxVarNum][maxVarNum];
        covDet = new double[taskNum];
        kld = new double[taskNum];
    }

    public double[] getKDL(int task) throws JMException {
        double[] kld = new double[taskNum];
        double tr, u;
        double s1, s2;

        cov[task] = getCov(task).getData();
        covDet[task] = getCovDet(task);
        covInv[task] = getCovInv(task).getData();

        int varNum = archives[task].get(0).getDecisionVariables().length;
        for (int i = 0; i < taskNum; i++){
            if (i == task)
                continue;

            cov[i] = getCov(i).getData();
            covDet[i] = getCovDet(i);
            covInv[i] = getCovInv(i).getData();

            tr = getTrace(task, i);
            u = getMul(task, i);
            s1 = Math.abs(0.5 * (tr + u - varNum + Math.log(covDet[task] / covDet[i])));

            tr = getTrace(i, task);
            u = getMul(i, task);
            s2 = Math.abs(0.5 * (tr + u - varNum + Math.log(covDet[i] / covDet[task])));
            
            kld[i] = 0.5 * (s1 + s2);
        }
        return kld;
    }

    private double getMul(int t1, int t2) throws JMException {
        double a[] = new double[maxVarNum];
        double sum;
        int i, j;

        for (i = 0; i < varNum; i++){
            sum = 0;
            for (j = 0; j < varNum; j++){
                sum += (archives[t1].getMeanOfIdx(j) - archives[t2].getMeanOfIdx(j)) * covInv[t1][j][i];
            }
            a[i] = sum;
        }

        sum = 0;
        for (i = 0; i < varNum; i++){
            sum += (archives[t1].getMeanOfIdx(i) - archives[t2].getMeanOfIdx(i)) * a[i];
        }

        return sum;
    }

    // 计算t1 cov_inv和t2 cov的乘积的迹
    private double getTrace(int t1, int t2) {
//        double sum;
//        double[][] result = new double[maxVarNum][maxVarNum];
//        for (int i = 0; i < varNum; i++){
//            for (int j = 0; j < varNum; j++){
//                sum = 0;
//                for (int l = 0; l < varNum; l++)
//                    sum += covInv[t1][i][l] * cov[t2][l][j];
//                result[i][j] = sum;
//            }
//        }
//
//        sum = 0;
//        for (int i = 0; i < varNum; i++)
//            sum += result[i][i];
//
//        return sum;
        RealMatrix m1 = new Array2DRowRealMatrix(covInv[t1]);
        RealMatrix m2 = new Array2DRowRealMatrix(cov[t2]);
        RealMatrix m = m1.multiply(m2);
        return m.getTrace();
    }

    // 计算协方差矩阵的逆
    private RealMatrix getCovInv(int task) throws JMException {
//        int[] is = new int[maxVarNum];
//        int[] js = new int[maxVarNum];
//        double d, p;
//        double[][] inv = new double[maxVarNum][maxVarNum];
//        for (int i = 0; i < varNum; i++) {
//            for (int j = 0; j < varNum; j++)
//                inv[i][j] = cov[task][i][j];
//        }
//
//        for (int k = 0; k < varNum; k++){
//            d = 0.0;
//            for (int i = k; i < varNum; i++){
//                for (int j = k; j < varNum; j++){
//                    p = Math.abs(inv[i][j]);
//                    if (p > d) { d = p; is[k] = i; js[k] = j; }
//                }
//            }
//
//            if (d + 1.0 == 1.0) {
//                System.out.println("error! (line 134");
//                System.exit(0);
//            }
//
//            if (is[k] != k){
//                for (int j = 0; j < varNum; j++){
//                    p = inv[k][j];
//                    inv[k][j] = inv[is[k]][j];
//                    inv[is[k]][j] = p;
//                }
//            }
//
//            if (js[k] != k){
//                for (int i = 0; i < varNum; i++){
//                    p = inv[i][k];
//                    inv[i][k] = inv[i][js[k]];
//                    inv[i][js[k]] = p;
//                }
//            }
//
//            inv[k][k] = 1.0 / inv[k][k];
//
//            for (int j = 0; j < varNum; j++)
//                if (j != k) inv[k][j] *= inv[k][k];
//            for (int i = 0; i < varNum; i++){
//                if (i != k){
//                    for (int j = 0; j < varNum; j++) {
//                        if (j != k)
//                            inv[i][j] = inv[i][j] - inv[i][k] * inv[k][j];
//                    }
//                }
//            }
//
//            for (int i = 0; i < varNum; i++){
//                if (i != k){
//                    inv[i][k] = -inv[i][k] * inv[k][k];
//                }
//            }
//
//        }
//
//        for (int k = varNum - 1; k >= 0; k--){
//            if (js[k] != k){
//                for (int j = 0; j < varNum; j++){
//                    p = inv[k][j];
//                    inv[k][j] = inv[js[k]][j];
//                    inv[js[k]][j] = p;
//                }
//            }
//
//            if (is[k] != k){
//                for (int i = 0; i < varNum; i++){
//                    p = inv[i][k];
//                    inv[i][k] = inv[i][is[k]];
//                    inv[i][is[k]] = p;
//                }
//            }
//        }
//        return inv;
        RealMatrix m = new Array2DRowRealMatrix(cov[task]);
        SingularValueDecomposition svd = new SingularValueDecomposition(m);
        DecompositionSolver solver = svd.getSolver();
        return solver.getInverse();
    }


    // 计算协方差矩阵的行列式
    private double getCovDet(int task) throws JMException {
//        int is, js;
//        double f, det, q, d;
//        double[][] a = new double[maxVarNum][maxVarNum];
//        for (int i = 0; i < varNum; i++)
//            for (int j = 0; j < varNum; j++)
//                a[i][j] = cov[task][i][j];
//
//        f = 1.0; det = 1.0;
//        for (int k = 0; k < varNum - 1; k++){
//            q = 0.0;
//            is = js = 0;
//            for (int i = k; i < varNum; i++){
//                for (int j = k; j < varNum; j++){
//                    d = Math.abs(a[i][j]);
//                    if (d > q) {
//                        q = d;
//                        is = i;
//                        js = j;
//                    }
//                }
//            }
//            if (q + 1.0 == 1.0) {
//                det = 0.0;
//                System.out.println("error! (line 221)");
//                System.exit(0);
//            }
//
//            if (is != k){
//                f = -f;
//                for (int j = k; j < varNum; j++){
//                    d = a[k][j];
//                    a[k][j] = a[is][j];
//                    a[is][j] = d;
//                }
//            }
//
//            if (js != k){
//                f = -f;
//                for (int i = k; i < varNum; i++){
//                    d = a[i][js];
//                    a[i][js] = a[i][k];
//                    a[i][k] = d;
//                }
//            }
//
//            det = det * a[k][k];
//            for (int i = k + 1; i < varNum; i++){
//                d = a[i][k] / a[k][k];
//                for (int j = k + 1; j < varNum; j++){
//                    a[i][j] = a[i][j] - d * a[k][j];
//                }
//            }
//        }
//
//        det = f * det * a[varNum - 1][varNum - 1];
//        if (det < 1e-3)
//            det = 1e-3;
//        return det;
        RealMatrix m = new Array2DRowRealMatrix(cov[task]);
        LUDecomposition L = new LUDecomposition(m);
        double det = L.getDeterminant();
        if (det < 1e-3)
            det = 1e-3;
        return det;
    }

    // 计算种群archives[task]中变量间的相关性（所有个体之和，取平均）
    private RealMatrix getCov(int task) throws JMException {
        Covariance c = new Covariance(archives[task].getMat());
        return c.getCovarianceMatrix();
    }
}
