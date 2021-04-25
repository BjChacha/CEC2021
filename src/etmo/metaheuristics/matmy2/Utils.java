//  Utils.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package etmo.metaheuristics.matmy2;

import etmo.core.Operator;
import etmo.core.Solution;
import etmo.encodings.variable.Real;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.sorting.SortingIdx;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.bytedeco.opencv.opencv_dnn.TanHLayer;
import org.deeplearning4j.datasets.iterator.DataSetIteratorSplitter;
import org.deeplearning4j.datasets.iterator.impl.ListDataSetIterator;
import org.deeplearning4j.datasets.iterator.impl.MnistDataSetIterator;
import org.deeplearning4j.datasets.iterator.loader.DataSetLoaderIterator;
import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.GradientNormalization;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.inputs.InputType;
import org.deeplearning4j.nn.conf.layers.BatchNormalization;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.listeners.PerformanceListener;
import org.nd4j.evaluation.classification.Evaluation;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.api.buffer.DataType;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.api.ops.impl.transforms.strict.Tan;
import org.nd4j.linalg.api.rng.distribution.Distribution;
import org.nd4j.linalg.api.rng.distribution.impl.NormalDistribution;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.dimensionalityreduction.PCA;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.learning.config.AdaGrad;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.learning.config.IUpdater;
import org.nd4j.linalg.learning.config.Sgd;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.nd4j.nativeblas.Nd4jCpu;

import javax.swing.*;
import java.io.IOException;
import java.util.Arrays;
import java.util.function.Supplier;


public class Utils {
    public static double distVector(double[] vector1, double[] vector2) {
        int dim = vector1.length;
        double sum = 0;
        for (int n = 0; n < dim; n++) {
            sum += (vector1[n] - vector2[n]) * (vector1[n] - vector2[n]);
        }
        return Math.sqrt(sum);
    } // distVector

    public static void minFastSort(double x[], int idx[], int n, int m) {
        for (int i = 0; i < m; i++) {
            for (int j = i + 1; j < n; j++) {
                if (x[i] > x[j]) {
                    double temp = x[i];
                    x[i] = x[j];
                    x[j] = temp;
                    int id = idx[i];
                    idx[i] = idx[j];
                    idx[j] = id;
                } // if
            }
        } // for

    } // minFastSort

    public static void maxFastSort(double x[], int idx[], int n, int m) {
        for (int i = 0; i < m; i++) {
            for (int j = i + 1; j < n; j++) {
                if (x[i] < x[j]) {
                    double temp = x[i];
                    x[i] = x[j];
                    x[j] = temp;
                    int id = idx[i];
                    idx[i] = idx[j];
                    idx[j] = id;
                } // if
            }
        } // for

    } // minFastSort

    public static void randomPermutation(int[] perm, int size) {
        int[] index = new int[size];
        boolean[] flag = new boolean[size];

        for (int n = 0; n < size; n++) {
            index[n] = n;
            flag[n] = true;
        }

        int num = 0;
        while (num < size) {
            int start = etmo.util.PseudoRandom.randInt(0, size - 1);
            // int start = int(size*nd_uni(&rnd_uni_init));
            while (true) {
                if (flag[start]) {
                    perm[num] = index[start];
                    flag[start] = false;
                    num++;
                    break;
                }
                if (start == (size - 1)) {
                    start = 0;
                } else {
                    start++;
                }
            }
        } // while
    } // randomPermutation

    public static void randomPermutation(int[] perm, int size, int range){
        if (size > range){
            System.out.println("randomPermutation error: size should be less than range.");
            return;
        }
        int[] index = new int[range];
        randomPermutation(index, range);

        for (int i = 0; i < perm.length; i++){
            perm[i] = index[i];
        }
    }

    static void QuickSort(double[] array, int[] idx, int from, int to) {
        if (from < to) {
            double temp = array[to];
            int tempIdx = idx[to];
            int i = from - 1;
            for (int j = from; j < to; j++) {
                if (array[j] <= temp) {
                    i++;
                    double tempValue = array[j];
                    array[j] = array[i];
                    array[i] = tempValue;
                    int tempIndex = idx[j];
                    idx[j] = idx[i];
                    idx[i] = tempIndex;
                }
            }
            array[to] = array[i + 1];
            array[i + 1] = temp;
            idx[to] = idx[i + 1];
            idx[i + 1] = tempIdx;
            QuickSort(array, idx, from, i);
            QuickSort(array, idx, i + 1, to);
        }
    }

    public static double[] softMaxOnlyPositive(double[] arr){
        double[] res = new double[arr.length];
        double sum = 0;
        double max = 0;
        for (int i = 0; i < arr.length; i++) {
            if (max < arr[i])
                max = arr[i];
        }
        for (int i = 0; i < res.length; i++){
            if (arr[i] > 0){
                res[i] = Math.exp(arr[i] / max);
                sum += res[i];
            }
            else{
                res[i] = 0;
            }
        }

        for (int i = 0; i < res.length; i++) {
            res[i] /= sum;
        }
        return res;
    }

    public static int roulette(double[] arr){
        double s = 0;
        double p = PseudoRandom.randDouble();
        int idx;
        // 轮盘赌算法
        for (idx = 0; idx < arr.length; idx++) {
            s += arr[idx];
            if (s >= p)
                break;
        }
        if (idx >= arr.length)
            idx = arr.length - 1;
        return idx;
    }

    public static double[] vectorMinus(double[] vec1, double[] vec2){
        if (vec1.length != vec2.length){
            System.out.println("Error: only can minus vectors with identical length.");
            return vec1;
        }

        double[] res = new double[vec1.length];
        for (int i = 0; i < res.length; i++)
            res[i] = vec1[i] - vec2[i];
        return res;
    }

    public static double calVectorAngle(double[] vec1, double[] vec2){
        if (vec1.length != vec2.length){
            System.out.println("Error: only can cal vectors with identical length.");
            return 1;
        }

        double v1 = 0;
        double v2 = 0;
        double dotP = 0;
        for (int i = 0; i < vec1.length; i++){
            v1 += vec1[i] * vec1[i];
            v2 += vec2[i] * vec2[i];
            dotP += vec1[i] * vec2[i];
        }
        v1 = Math.sqrt(v1);
        v2 = Math.sqrt(v2);

        double res = Math.abs(dotP / (v1 * v2));
        res = Math.acos(res);

        if (Double.isNaN(res) && v1 == 0)
            res = Double.MAX_VALUE;

        return res;
    }

    // PCA
    public static RealMatrix PCAGetBaseVectors(double[][] mat){
        // 归一化+缩放，使均值为零，且最大最小值之差为1.
        MatNormalize(mat);
        // 获得协方差矩阵
        RealMatrix covarianceMatrix = new Covariance(MatrixUtils.createRealMatrix(mat)).getCovarianceMatrix();

        // 奇异值分解
        SingularValueDecomposition svd = new SingularValueDecomposition(covarianceMatrix);
        // 取前95%信息的奇异值数量
        double[] singleValues = svd.getSingularValues();
        double sum = 0;
        for (double value: singleValues)
            sum += value;
        int d = 0;
        double curSum = 0;
        for (double value: singleValues){
            curSum += value;
            d ++;
            if (curSum >= sum * 0.95)
                break;
        }

        RealMatrix U = svd.getU();
        RealMatrix baseVectors = U.getSubMatrix(0,U.getRowDimension() - 1,0,d - 1);

        // 输入：n x D   输出：d x D   降维：D -> d
        return baseVectors;
    }

    public static RealMatrix PCAGetBaseVectors(double[][] mat, int d){
        // 归一化+缩放，使均值为零，且最大最小值之差为1.
        MatNormalize(mat);
        // 获得协方差矩阵
        RealMatrix covarianceMatrix = new Covariance(MatrixUtils.createRealMatrix(mat)).getCovarianceMatrix();

        // 奇异值分解
        SingularValueDecomposition svd = new SingularValueDecomposition(covarianceMatrix);
        RealMatrix U = svd.getU();
        RealMatrix baseVectors = U.getSubMatrix(0,U.getRowDimension() - 1,0,d - 1);

        // 输入：n x D   输出：d x D   降维：D -> d
        return baseVectors;
    }

    static void MatNormalize(double[][] mat){
        int row = mat.length;
        int col = mat[0].length;

        double[] means = new double[col];
        double[] max = new double[col];
        double[] min = new double[col];
        for (int i = 0; i < col; i++) {
            float sum = 0;
            max[i] = Double.MIN_VALUE;
            min[i] = Double.MAX_VALUE;
            // 获得每一列的均值、最大值、最小值
            for (int j = 0; j < row; j++) {
                if (max[i] < mat[j][i])
                    max[i] = mat[j][i];
                if (min[i] > mat[j][i])
                    min[i] = mat[j][i];
                sum += mat[j][i];
            }
            means[i] = sum / row;
            // 数据预处理
            for (int j = 0; j < row; j++) {
                if (max[i] == min[i])
                    mat[j][i] = 0;
                else
                    mat[j][i] = (mat[j][i] - means[i]) / (max[i] - min[i]);
            }
        }
    }

    static void MatReform(double[][] mat, double minimun, double maximun){
        int col = mat[0].length;
        int row = mat.length;

        double[] max = new double[col];
        double[] min = new double[col];
        for (int i = 0; i < col; i++){
            max[i] = Double.MIN_VALUE;
            min[i] = Double.MAX_VALUE;
            for (int j = 0; j < row; j++){
                if (max[i] < mat[j][i])
                    max[i] = mat[j][i];
                if (min[i] > mat[j][i])
                    min[i] = mat[j][i];
            }
            for (int j = 0; j < row; j++){
                if ((max[i] - min[i]) == 0)
                    mat[j][i] = 0;
                else
                    mat[j][i] = (mat[j][i] - min[i]) / (max[i] - min[i]) * (maximun - minimun) + minimun;
                if (mat[j][i] < 0)
                    System.out.println("No!!!");
            }
        }
    }

    public static double[][] MappingViaPCA(double[][] mat1, double[][] mat2){
        // 迁移方向: mat1->mat2
        RealMatrix realMatrix1 = MatrixUtils.createRealMatrix(mat1);

        RealMatrix baseVectors1 = PCAGetBaseVectors(mat1);
        RealMatrix baseVectors2 = PCAGetBaseVectors(mat2, baseVectors1.getColumnDimension());

        RealMatrix M = baseVectors1.transpose().multiply(baseVectors2);

        RealMatrix resMatrix = realMatrix1.multiply(baseVectors1);
        resMatrix = resMatrix.multiply(M);
        resMatrix = resMatrix.multiply(baseVectors2.transpose());

        double[][] res = resMatrix.getData();
        MatReform(res, 0, 1);
        return res;
    }

    public static double[][] MappingViaPCA1(double[][] mat1, double[][] mat2){
        // 迁移方向: mat1->mat2
        RealMatrix realMatrix1 = MatrixUtils.createRealMatrix(mat1);
        RealMatrix realMatrix2 = MatrixUtils.createRealMatrix(mat2);
        RealMatrix baseVectors1 = PCAGetBaseVectors(mat1);
        RealMatrix baseVectors2 = PCAGetBaseVectors(mat2, baseVectors1.getColumnDimension());

        NDArray input = new NDArray(realMatrix1.multiply(baseVectors1).getData());
        NDArray output = new NDArray(realMatrix2.multiply(baseVectors2).getData());

        AutoEncoder ae = new AutoEncoder(input, output);
        ae.train();
        double[][] res = ae.predict(input);

        MatReform(res, 0, 1);
        return res;
    }

    public static double[][] MappingViaAE(double[][] mat1, double[][] mat2, int size){
        AutoEncoder ae = new AutoEncoder(mat1, mat2);
        ae.train();
        double[][] goodPops = Arrays.copyOfRange(mat1, 0, 25);
        double[][] res = ae.predict(Sample(goodPops, size));

        MatReform(res, 0, 1);
        return res;
    }

    public static double[][] MappingViaPE(double[][] mat1, double[][] mat2, int size){
        // 迁移方向: mat1->mat2
        RealMatrix realMatrix1 = MatrixUtils.createRealMatrix(mat1);
        RealMatrix realMatrix2 = MatrixUtils.createRealMatrix(mat2);
        RealMatrix baseVectors1 = PCAGetBaseVectors(mat1);
        RealMatrix baseVectors2 = PCAGetBaseVectors(mat2, baseVectors1.getColumnDimension());

        // 降维后的输入输出
        NDArray input = new NDArray(realMatrix1.multiply(baseVectors1).getData());
        NDArray output = new NDArray(realMatrix2.multiply(baseVectors2).getData());

        // 训练自编码器
        AutoEncoder ae = new AutoEncoder(input, output);
        ae.train();

        // 取输入样本前25个个体，建立正态分布，并采样。
        double[][] goodPops = Arrays.copyOfRange(mat1, 0, 25);
        RealMatrix goodPopsMat = MatrixUtils.createRealMatrix(goodPops);
        // 预测前需要先降维
        goodPopsMat = goodPopsMat.multiply(baseVectors1);
        double[][] predicts = ae.predict(Sample(goodPopsMat.getData(), size));
        // 输出个体要升维复原
        RealMatrix recoveredPopsMat = MatrixUtils.createRealMatrix(predicts).multiply(baseVectors2.transpose());
        double [][] res = recoveredPopsMat.getData();
        MatReform(res, 0, 1);
        return res;
    }


    public static double[][] MappingViaGAN(double[][] mat1, double[][] mat2, int size){
        // mat1 -> mat2
        int d1 = mat1[0].length;
        int d2 = mat2[0].length;
        int hG = (d1 + d2);
        int hD = (int) (d2 * 1);

        int SEED = 42;
        int EPOCHS = 10;
        int LATENT_DIM = d1;
        double LR = 5e-2;

        IUpdater UPDATER = Adam.builder().learningRate(LR).beta1(0.5).build();

        Supplier<MultiLayerNetwork> genSupplier = () -> {
            return new MultiLayerNetwork(new NeuralNetConfiguration.Builder()
                    .weightInit(WeightInit.XAVIER)
                    .activation(Activation.IDENTITY)
                    .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                    .gradientNormalizationThreshold(100)
                    .list()
                    .layer(new DenseLayer.Builder().nIn(d1).nOut(hG).weightInit(WeightInit.NORMAL).build())
                    .layer(new DenseLayer.Builder().nIn(hG).nOut(hG).weightInit(WeightInit.NORMAL).build())
                    .layer(new DenseLayer.Builder().nIn(hG).nOut(d2).activation(Activation.TANH).build())
                    .build());
        };

        GAN.DiscriminatorProvider discriminatorProvider = (updater) -> {
            return new MultiLayerNetwork(new NeuralNetConfiguration.Builder()
                    .seed(SEED)
                    .updater(updater)
                    .weightInit(WeightInit.XAVIER)
                    .activation(Activation.IDENTITY)
                    .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                    .gradientNormalizationThreshold(100)
                    .list()
                    .layer(new DenseLayer.Builder().nIn(d2).nOut(hD).build())
                    .layer(new DenseLayer.Builder().nIn(hD).nOut(hD).build())
//                    .layer(new DenseLayer.Builder().nIn(hD).nOut(hD).build())
                    .layer(new OutputLayer.Builder().lossFunction(LossFunctions.LossFunction.XENT).nIn(hD).nOut(1).
                            activation(Activation.SIGMOID).build())
                    .build());
        };

        GAN gan = new GAN.Builder()
                .generator(genSupplier)
                .discriminator(discriminatorProvider)
                .latentDimension(LATENT_DIM)
                .seed(SEED)
                .updater(UPDATER)
                .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                .gradientNormalizationThreshold(100)
                .build();

//        Nd4j.getMemoryManager().setAutoGcWindow(15 * 1000);

        gan.getGenerator().setListeners(new PerformanceListener(1, true));
        gan.getDiscriminator().setListeners(new PerformanceListener(1, true));

        INDArray features = new NDArray((new NDArray(mat2)).toFloatMatrix());
        INDArray labels = Nd4j.zeros(mat2.length, 1);

        DataSet trainData = new DataSet(features, labels);

        for (int e = 0; e < EPOCHS; e++){
            gan.fit(trainData);
        }

//        // DEBUG: evaluate
//        // existed target individuals
//        INDArray fake = Nd4j.rand(DataType.DOUBLE, new int[]{mat2.length, mat2[0].length});
//        INDArray fakeLabels = Nd4j.ones(mat2.length, 1);
//        DataSet fakeSet = new DataSet(fake, fakeLabels);
//        INDArray real = new NDArray(mat2);
//        real = real.muli(2).subi(1);
//        INDArray realLabels = Nd4j.zeros(mat2.length, 1);
//        DataSet realSet = new DataSet(real, realLabels);
//        DataSet testSet = DataSet.merge(Arrays.asList(fakeSet, realSet));
//        DataSetIterator testIterator = new ListDataSetIterator(testSet.asList());
//        Evaluation eval = gan.discriminator.evaluate(testIterator);
//        System.out.println("existed: \n"+eval.stats());

//        // generated individuals
//        INDArray noise = Nd4j.rand(new int[]{mat2.length, mat2[0].length});
//        INDArray test3 = gan.generator.output(noise);
//        INDArray testLabels3 = Nd4j.ones(mat2.length, 1);
//        DataSet testDataset3 = new DataSet(test3, testLabels3);
//        DataSetIterator testIterator3 = new ListDataSetIterator(testDataset3.asList());
//        Evaluation eval3 = gan.discriminator.evaluate(testIterator3);
//        System.out.println("fake: \n"+eval3.stats());

        INDArray samples = new NDArray(mat1);
        INDArray predicts = gan.getGenerator().output(samples);

        gan = null;
        Nd4j.getMemoryManager().invokeGc();

        INDArray res = predicts.addi(1).divi(2);
        return res.toDoubleMatrix();
    }


    public static double[][] MappingViaGANMixCrossover(double[][] mat1, double[][] mat2, boolean[] isGAN, int size){
        // mat1 -> mat2x
        int d1 = mat1[0].length;
        int d2 = mat2[0].length;
        int hG = (d1+d2)/2;
        int hD = (d1+d2)/2;

        int SEED = 42;
        int EPOCHS = 3;
        int LATENT_DIM = d1;
        double LR = 1e-1;

        IUpdater UPDATER = Adam.builder().learningRate(LR).beta1(0.5).build();

        Supplier<MultiLayerNetwork> genSupplier = () -> {
            return new MultiLayerNetwork(new NeuralNetConfiguration.Builder()
                    .weightInit(WeightInit.XAVIER)
                    .activation(Activation.IDENTITY)
                    .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                    .gradientNormalizationThreshold(100)
                    .list()
                    .layer(new DenseLayer.Builder().nIn(d1).nOut(hG).weightInit(WeightInit.NORMAL).build())
//                    .layer(new DenseLayer.Builder().nIn(hG).nOut(hG).weightInit(WeightInit.NORMAL).build())
                    .layer(new DenseLayer.Builder().nIn(hG).nOut(d2).activation(Activation.TANH).build())
                    .build());
        };

        GAN.DiscriminatorProvider discriminatorProvider = (updater) -> {
            return new MultiLayerNetwork(new NeuralNetConfiguration.Builder()
                    .seed(SEED)
                    .updater(updater)
                    .weightInit(WeightInit.XAVIER)
                    .activation(Activation.IDENTITY)
                    .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                    .gradientNormalizationThreshold(100)
                    .list()
                    .layer(new DenseLayer.Builder().nIn(d2).nOut(hD).build())
//                    .layer(new DenseLayer.Builder().nIn(hD).nOut(hD).build())
//                    .layer(new DenseLayer.Builder().nIn(hD).nOut(hD).build())
                    .layer(new OutputLayer.Builder().lossFunction(LossFunctions.LossFunction.XENT).nIn(hD).nOut(1).
                            activation(Activation.SIGMOID).build())
                    .build());
        };

        GAN gan = new GAN.Builder()
                .generator(genSupplier)
                .discriminator(discriminatorProvider)
                .latentDimension(LATENT_DIM)
                .seed(SEED)
                .updater(UPDATER)
                .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                .gradientNormalizationThreshold(100)
                .build();

//        Nd4j.getMemoryManager().setAutoGcWindow(15 * 1000);

        gan.getGenerator().setListeners(new PerformanceListener(1, true));
        gan.getDiscriminator().setListeners(new PerformanceListener(1, true));

        INDArray features = new NDArray((new NDArray(mat2)).toFloatMatrix());
        INDArray labels = Nd4j.zeros(mat2.length, 1);

        DataSet trainData = new DataSet(features, labels);

        for (int e = 0; e < EPOCHS; e++){
            gan.fit(trainData);
        }

        // DEBUG: evaluate
        // existed target individuals
        INDArray fake = Nd4j.rand(DataType.DOUBLE, new int[]{mat2.length, mat2[0].length});
        INDArray fakeLabels = Nd4j.ones(mat2.length, 1);
        DataSet fakeSet = new DataSet(fake, fakeLabels);
        INDArray real = new NDArray(mat2);
        real = real.muli(2).subi(1);
        INDArray realLabels = Nd4j.zeros(mat2.length, 1);
        DataSet realSet = new DataSet(real, realLabels);
        DataSet testSet = DataSet.merge(Arrays.asList(fakeSet, realSet));
        DataSetIterator testIterator = new ListDataSetIterator(testSet.asList());
        Evaluation eval = gan.discriminator.evaluate(testIterator);
//        System.out.println("existed: \n"+eval.stats());

//        System.out.println(eval.accuracy() +"\t"+ eval.recall() +"\t"+ eval.precision() +"\t"+ eval.f1());
        if (eval.f1() > 0.8)
            isGAN[0] = true;

        INDArray samples = new NDArray(mat1);
        INDArray predicts = gan.getGenerator().output(samples);
        INDArray res = predicts.addi(1).divi(2);

        // clear memory
        gan = null;
        Nd4j.getMemoryManager().invokeGc();

        return res.toDoubleMatrix();
    }



    public static double[] MappingViaGAN(Solution srcIndividual, double[][] targetDomain) throws JMException {
        // mat1 -> mat2
        int shape = targetDomain.length;
        int d1 = srcIndividual.getDecisionVariables().length;
        int d2 = targetDomain[0].length;

        int EPOCHS = 200;
        int LATENT_DIM = d1;
        double LR = 1e-4;

        GAN gan = GANBuilder(d2, LATENT_DIM, LR);

        INDArray features = new NDArray((new NDArray(targetDomain)).toFloatMatrix());
        INDArray labels = Nd4j.zeros(shape, 1);

        DataSet trainData = new DataSet(features, labels);

        for (int e = 0; e < EPOCHS; e++){
            gan.fit(trainData);
        }

        double[][] srcI = new double[1][];
        srcI[0] = srcIndividual.getDecisionVariablesInDouble();
        INDArray samples = new NDArray(srcI);
        INDArray predicts = gan.getGenerator().output(samples);

        gan = null;
        Nd4j.getMemoryManager().invokeGc();

        INDArray res = predicts.addi(1).divi(2);

        return res.toDoubleMatrix()[0];
    }


    public static double[][] GANTest(double[][] mat1, double[][] mat2, int HG, int HD, double LRG, double LRD, int E) throws IOException {
        // mat1 -> mat2x
        int d1 = mat1[0].length;
        int d2 = mat2[0].length;
        int hG = HG;
        int hD = HD;

        int SEED = 42;
        int EPOCHS = E;
        int LATENT_DIM = d1;
        double gLR = LRG;
        double dLR = LRD;

        IUpdater dUPDATER = Adam.builder().learningRate(dLR).beta1(0.5).build();
        IUpdater gUPDATER = Adam.builder().learningRate(gLR).beta1(0.5).build();

        Supplier<MultiLayerNetwork> genSupplier = () -> {
            return new MultiLayerNetwork(new NeuralNetConfiguration.Builder()
                    .weightInit(WeightInit.XAVIER)
                    .activation(Activation.IDENTITY)
                    .updater(gUPDATER)
//                    .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
//                    .gradientNormalizationThreshold(100)
                    .list()
                    .layer(new DenseLayer.Builder().nIn(d1).nOut(hG).weightInit(WeightInit.NORMAL).build())
                    .layer(new DenseLayer.Builder().nIn(hG).nOut(hG).weightInit(WeightInit.NORMAL).build())
                    .layer(new DenseLayer.Builder().nIn(hG).nOut(d2).activation(Activation.TANH).build())
                    .build());
        };

        GAN.DiscriminatorProvider discriminatorProvider = (updater) -> {
            return new MultiLayerNetwork(new NeuralNetConfiguration.Builder()
                    .seed(SEED)
                    .updater(updater)
                    .weightInit(WeightInit.XAVIER)
                    .activation(Activation.IDENTITY)
//                    .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
//                    .gradientNormalizationThreshold(100)
                    .list()
                    .layer(new DenseLayer.Builder().nIn(d2).nOut(hD).build())
                    .layer(new DenseLayer.Builder().nIn(hD).nOut(hD).build())
                    .layer(new OutputLayer.Builder().lossFunction(LossFunctions.LossFunction.XENT).nIn(hD).nOut(1).
                            activation(Activation.SIGMOID).build())
                    .build());
        };

        GAN gan = new GAN.Builder()
                .generator(genSupplier)
                .discriminator(discriminatorProvider)
                .latentDimension(LATENT_DIM)
                .seed(SEED)
                .updater(dUPDATER)
                .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                .gradientNormalizationThreshold(100)
                .build();

//        Nd4j.getMemoryManager().setAutoGcWindow(15 * 1000);

        gan.getGenerator().setListeners(new PerformanceListener(1, true));
        gan.getDiscriminator().setListeners(new PerformanceListener(1, true));

        INDArray trueFeature = new NDArray((new NDArray(mat2)).toFloatMatrix());
        INDArray trueLabel = Nd4j.zeros(mat2.length, 1);
        DataSet trueData = new DataSet(trueFeature, trueLabel);

        INDArray fakeFeature = new NDArray((new NDArray(mat1)).toFloatMatrix());
        INDArray fakeLabel = Nd4j.zeros(mat1.length, 1);
        DataSet fakeData = new DataSet(fakeFeature, fakeLabel);

        for (int e = 0; e < EPOCHS; e++){
//            gan.fit2(trueData, fakeData);
            gan.fit(trueData);
        }


        // DEBUG: evaluate
        // existed target individuals
        INDArray fake = Nd4j.rand(DataType.DOUBLE, new int[]{mat2.length, mat2[0].length});
        fake = fake.muli(2).subi(1);
        INDArray fakeLabels = Nd4j.ones(mat2.length, 1);
        DataSet fakeSet = new DataSet(fake, fakeLabels);
        INDArray real = new NDArray(mat2);
        real = real.muli(2).subi(1);
        INDArray realLabels = Nd4j.zeros(mat2.length, 1);
        DataSet realSet = new DataSet(real, realLabels);
        DataSet testSet = DataSet.merge(Arrays.asList(fakeSet, realSet));
        DataSetIterator testIterator = new ListDataSetIterator(testSet.asList());
        Evaluation eval = gan.discriminator.evaluate(testIterator);
//        System.out.println(eval.stats());
//        System.out.println(eval.accuracy() +"\t"+ eval.recall() +"\t"+ eval.precision() +"\t"+ eval.f1());

        INDArray samples = new NDArray(mat1);
        INDArray predicts = gan.getGenerator().output(samples);
        INDArray res = predicts.addi(1).divi(2);

        // clear memory
        gan = null;
        Nd4j.getMemoryManager().invokeGc();

        return res.toDoubleMatrix();
    }

    static GAN GANBuilder(int d, int dNoise, double lr){
        int SEED = 42;

        IUpdater GenUPDATER = Adam.builder().learningRate(lr*4).beta1(0.5).beta2(0.999).build();
        IUpdater DisUPDATER = Adam.builder().learningRate(lr).beta1(0.5).beta2(0.999).build();

        Supplier<MultiLayerNetwork> genSupplier = () -> {
            return new MultiLayerNetwork(new NeuralNetConfiguration.Builder()
                    .weightInit(WeightInit.XAVIER)
                    .updater(GenUPDATER)
                    .activation(Activation.IDENTITY)
                    .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                    .gradientNormalizationThreshold(100)
                    .list()
                    .layer(new DenseLayer.Builder().nIn(dNoise).nOut(d).hasBias(true).weightInit(WeightInit.UNIFORM).build())
                    .layer(new BatchNormalization.Builder().nIn(d).nOut(d).build())
                    .layer(new DenseLayer.Builder().nIn(d).nOut(d).hasBias(true).weightInit(WeightInit.UNIFORM).build())
                    .layer(new BatchNormalization.Builder().nIn(d).nOut(d).build())
                    .layer(new DenseLayer.Builder().nIn(d).nOut(d).hasBias(true).weightInit(WeightInit.UNIFORM).build())
                    .layer(new BatchNormalization.Builder().nIn(d).nOut(d).build())
                    .layer(new DenseLayer.Builder().nIn(d).nOut(d).activation(Activation.SIGMOID).build())
                    .build());
        };

        GAN.DiscriminatorProvider discriminatorProvider = (updater) -> {
            return new MultiLayerNetwork(new NeuralNetConfiguration.Builder()
                    .seed(SEED)
                    .updater(updater)
                    .weightInit(WeightInit.XAVIER)
                    .activation(Activation.IDENTITY)
                    .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                    .gradientNormalizationThreshold(100)
                    .list()
                    .layer(new DenseLayer.Builder().nIn(d).nOut(d).hasBias(true).weightInit(WeightInit.UNIFORM).build())
                    .layer(new OutputLayer.Builder().nIn(d).nOut(1).hasBias(true)
                            .lossFunction(LossFunctions.LossFunction.XENT).activation(Activation.SIGMOID).build())
                    .build());
        };

        GAN gan = new GAN.Builder()
                .generator(genSupplier)
                .discriminator(discriminatorProvider)
                .latentDimension(dNoise)
                .seed(SEED)
                .updater(DisUPDATER)
                .gradientNormalization(GradientNormalization.RenormalizeL2PerLayer)
                .gradientNormalizationThreshold(100)
                .build();

        return gan;
    }


    public static double[][] Sample(double[][] mat, int size){
        NDArray input = new NDArray(mat);
        double[] means = input.mean(0).toDoubleVector();
        double[] stds = input.var(0).toDoubleVector();

        double[][] res = new double[size][input.columns()];
        for (int j = 0; j < input.columns(); j ++){
            if (stds[j] < 0.1)
                stds[j] = 0.1;
            Distribution d = new NormalDistribution(means[j], stds[j]);
            for (int i = 0; i < size; i ++){
                res[i][j] = d.sample();
                if (res[i][j] < 0)
                    res[i][j] = 0;
                else if (res[i][j] > 1)
                    res[i][j] = 1;
            }
        }
        return res;
    }
}