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

    private SolutionSet[] archives;

    public KLD(ProblemSet problemSet, SolutionSet[] solutionSet) {
        archives = solutionSet;
        taskNum = problemSet.size();
        varNum = Integer.MAX_VALUE;
        maxVarNum = 0;
        
        for (int i = 0; i < archives.length; i++) {
            varNum = Math.min(varNum, archives[i].get(0).numberOfVariables());
            maxVarNum = Math.max(maxVarNum, archives[i].get(0).numberOfVariables());
        }
        cov = new double[taskNum][maxVarNum][maxVarNum];
        covInv = new double[taskNum][maxVarNum][maxVarNum];
        covDet = new double[taskNum];
    }

    public KLD(int taskNum, SolutionSet[] solutionSet) {
        archives = solutionSet;

        varNum = Integer.MAX_VALUE;
        maxVarNum = 0;
        for (int i = 0; i < archives.length; i++) {
            varNum = Math.min(varNum, archives[i].get(0).numberOfVariables());
            maxVarNum = Math.max(maxVarNum, archives[i].get(0).numberOfVariables());
        }
        this.taskNum = taskNum;
        cov = new double[taskNum][maxVarNum][maxVarNum];
        covInv = new double[taskNum][maxVarNum][maxVarNum];
        covDet = new double[taskNum];
    }

    public double[] getKDL(int task) throws JMException {
        double[] kld = new double[taskNum];
        double tr, u;
        double s1, s2;

        cov[task] = getCov(task).getData();
        covDet[task] = getCovDet(task);
        covInv[task] = getCovInv(task).getData();

        int varNum = archives[task].get(0).getDecisionVariables().length;
        for (int i = 0; i < taskNum; i++) {
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

        for (i = 0; i < varNum; i++) {
            sum = 0;
            for (j = 0; j < varNum; j++) {
                sum += (archives[t1].getMeanOfIdx(j) - archives[t2].getMeanOfIdx(j)) * covInv[t1][j][i];
            }
            a[i] = sum;
        }

        sum = 0;
        for (i = 0; i < varNum; i++) {
            sum += (archives[t1].getMeanOfIdx(i) - archives[t2].getMeanOfIdx(i)) * a[i];
        }

        return sum;
    }

    // 计算t1 cov_inv和t2 cov的乘积的迹
    private double getTrace(int t1, int t2) {
        RealMatrix m1 = new Array2DRowRealMatrix(covInv[t1]);
        RealMatrix m2 = new Array2DRowRealMatrix(cov[t2]);
        RealMatrix m = m1.multiply(m2);
        return m.getTrace();
    }

    // 计算协方差矩阵的逆
    private RealMatrix getCovInv(int task) throws JMException {
        RealMatrix m = new Array2DRowRealMatrix(cov[task]);
        SingularValueDecomposition svd = new SingularValueDecomposition(m);
        DecompositionSolver solver = svd.getSolver();
        return solver.getInverse();
    }

    // 计算协方差矩阵的行列式
    private double getCovDet(int task) throws JMException {
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
