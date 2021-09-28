package etmo.metaheuristics.matmy3.models;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.deeplearning4j.datasets.iterator.impl.ListDataSetIterator;
import org.nd4j.autodiff.listeners.impl.ScoreListener;
import org.nd4j.autodiff.listeners.impl.UIListener;
import org.nd4j.autodiff.listeners.records.History;
import org.nd4j.autodiff.samediff.SDIndex;
import org.nd4j.autodiff.samediff.SDVariable;
import org.nd4j.autodiff.samediff.SameDiff;
import org.nd4j.autodiff.samediff.TrainingConfig;
import org.nd4j.evaluation.classification.Evaluation;
import org.nd4j.evaluation.classification.Evaluation.Metric;
import org.nd4j.linalg.api.buffer.DataType;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.weightinit.impl.XavierInitScheme;

import shapeless.DataT;

public class CoralClassifier {

    SameDiff sd;
    int epochs;
    double lr;
    double[][] trgDataMat;

    public CoralClassifier(int epochs, double lr, double[][] trgData) {
        this.epochs = epochs;
        this.lr = lr;
        this.trgDataMat = trgData;
    }

    public void train(double[][] features, double[][] labels) {
        int srcN = features.length;
        int trgN = trgDataMat.length;
        int inD = features[0].length;
        int hD1 = inD / 2;
        int hD2 = hD1 / 2;
        int outD = 2;

        double lambda1 = 1e-2;
        double lambda2 = 1e-2;
        double lambda3 = 1e-2;

        this.sd = SameDiff.create();

        // INDArray tmp = new NDArray(trgDataMat);
        // INDArray trgData = tmp.castTo(DataType.DOUBLE);
        INDArray srcData = new NDArray(features);
        // INDArray srcLabel = new NDArray(new double[][] {labels});
        INDArray srcLabel = new NDArray(labels);
        // srcLabel.transposei();

        SDVariable srcInput = sd.placeHolder("srcInput", DataType.DOUBLE, -1, inD);
        // SDVariable trgInput = sd.placeHolder("trgInput", DataType.DOUBLE, -1, inD);
        SDVariable label = sd.placeHolder("label", DataType.DOUBLE, -1, outD);

        SDVariable Wsrc1 = sd.var("ws1", new XavierInitScheme('c', inD, hD1), DataType.DOUBLE, inD, hD1);
        SDVariable Bsrc1 = sd.zero("bs1", hD1);
        SDVariable Wsrc2 = sd.var("ws2", new XavierInitScheme('c', hD1, hD2), DataType.DOUBLE, hD1, hD2);
        SDVariable Bsrc2 = sd.zero("bs2", hD2);
        SDVariable Wsrc3 = sd.var("ws3", new XavierInitScheme('c', hD2, outD), DataType.DOUBLE, hD2, outD);
        SDVariable Bsrc3 = sd.zero("bs3", outD);

        SDVariable Wtrg1 = sd.var("wt1", new XavierInitScheme('c', inD, hD1), DataType.DOUBLE, inD, hD1);
        SDVariable Btrg1 = sd.zero("bt1", hD1);
        SDVariable Wtrg2 = sd.var("wt2", new XavierInitScheme('c', hD1, hD2), DataType.DOUBLE, hD1, hD2);
        SDVariable Btrg2 = sd.zero("bt2", hD2);
        SDVariable Wtrg3 = sd.var("wt3", new XavierInitScheme('c', hD2, outD), DataType.DOUBLE, hD2, outD);
        SDVariable Btrg3 = sd.zero("bt3", outD);

        // SDVariable sn = sd.constant("n", DataType.INT64, n);
        // SDVariable sind = sd.constant("inD", DataType.INT64, inD);
        // SDVariable shd1 = sd.constant("hD1", DataType.INT64, hD1);
        // SDVariable shd2 = sd.constant("hD2", DataType.INT64, hD2);

        // SDVariable trgDomainData = sd.var(trgData);

        // compute flow
        SDVariable ws1 = sd.nn().tanh(srcInput.mmul(Wsrc1).add(Bsrc1));
        SDVariable ws2 = sd.nn().tanh(ws1.mmul(Wsrc2).add(Bsrc2));
        // SDVariable srcOutput = sd.nn().sigmoid("srcOutput", ws2.mmul(Wsrc3).add(Bsrc3));
        SDVariable srcOutput = ws2.mmul(Wsrc3).add("srcOutput", Bsrc3);
        
        // SDVariable wtd1 = sd.nn().tanh(trgDomainData.mmul(Wtrg1).add(Btrg1));
        // SDVariable wtd2 = sd.nn().tanh(wtd1.mmul(Wtrg2).add(Btrg2));
        // SDVariable wtd3 = wtd2.mmul(Wtrg3).add(Btrg3);

        // SDVariable wt1 = sd.nn().tanh(trgInput.mmul(Wtrg1).add(Btrg1));
        // SDVariable wt2 = sd.nn().tanh(wt1.mmul(Wtrg2).add(Btrg2));
        // // SDVariable trgOutput = sd.nn().sigmoid("trgOutput", wt2.mmul(Wsrc3).add(Bsrc3));
        // SDVariable trgOutput = wt2.mmul(Wtrg3).add("trgOutput", Btrg3);

        // SDVariable coralLoss1 = coralLoss("coralLoss1", sd, ws1, wtd1, srcN, trgN, hD1);
        // SDVariable coralLoss2 = coralLoss("coralLoss2", sd, ws2, wtd2, srcN, trgN, hD2);
        
        // SDVariable oneVector1 = sd.one("ones1", DataType.DOUBLE, srcN, 1);
        // SDVariable oneVector2 = sd.one("ones2", DataType.DOUBLE, trgN, 1);

        // SDVariable srcMean1 = 
        //     sd.transpose(sd.transpose(oneVector1).mmul(ws1))
        //     .mmul(sd.transpose(oneVector1).mmul(ws1)).div(srcN);
        // SDVariable trgMean1 = 
        //     sd.transpose(sd.transpose(oneVector2).mmul(wtd1))
        //     .mmul(sd.transpose(oneVector2).mmul(wtd1)).div(trgN);
        // SDVariable srcCov1 = 
        //     sd.transpose(ws1).mmul(ws1).sub(srcMean1).div(srcN - 1);
        // SDVariable trgCov1 = 
        //     sd.transpose(wtd1).mmul(wtd1).sub(trgMean1).div(trgN - 1);
        // SDVariable coralLoss1 = srcCov1.sub(trgCov1).norm2().div(hD1 * 4);
        
        // SDVariable srcMean2 = 
        //     sd.transpose(sd.transpose(oneVector1).mmul(ws2))
        //     .mmul(sd.transpose(oneVector1).mmul(ws2)).div(srcN);
        // SDVariable trgMean2 = 
        //     sd.transpose(sd.transpose(oneVector2).mmul(wtd2))
        //     .mmul(sd.transpose(oneVector2).mmul(wtd2)).div(trgN);
        // SDVariable srcCov2 = 
        //     sd.transpose(ws2).mmul(ws2).sub(srcMean2).div(srcN - 1);
        // SDVariable trgCov2 = 
        //     sd.transpose(wtd2).mmul(wtd2).sub(trgMean2).div(trgN - 1);
        // SDVariable coralLoss2 = srcCov2.sub(trgCov2).norm2().div(hD2 * 4);

        // SDVariable srcMean3 = 
        //     sd.transpose(sd.transpose(oneVector1).mmul(srcOutput))
        //     .mmul(sd.transpose(oneVector1).mmul(srcOutput)).div(srcN);
        // SDVariable trgMean3 = 
        //     sd.transpose(sd.transpose(oneVector2).mmul(wtd3))
        //     .mmul(sd.transpose(oneVector2).mmul(wtd3)).div(trgN);
        // SDVariable srcCov3 = 
        //     sd.transpose(srcOutput).mmul(srcOutput).sub(srcMean3).div(srcN - 1);
        // SDVariable trgCov3 = 
        //     sd.transpose(wtd3).mmul(wtd3).sub(trgMean3).div(trgN - 1);
        // SDVariable coralLoss3 = srcCov3.sub(trgCov3).norm2().div(outD * 4);
        
        SDVariable XENTLoss = sd.loss().softmaxCrossEntropy("loss", label, srcOutput, null);
        // SDVariable XENTLoss = sd.mean("xentLoss",
        //     label.mul(-1.0).mul(sd.math.log(sd.constant(1.0).sub(srcOutput)))
        //     .sub(sd.constant(1.0).sub(label).mul(sd.math.log(srcOutput))));

        // SDVariable totalLoss = XENTLoss.add(coralLoss1.mul(lambda1))
        //     .add(coralLoss2.mul(lambda2))
        //     .add(coralLoss3.mul(lambda3));
            
            sd.setLossVariables(XENTLoss);
        // sd.setLossVariables(totalLoss);

        Evaluation evaluation = new Evaluation();

        TrainingConfig config = TrainingConfig.builder()
            // .l2(2e-5)
            .updater(new Adam(lr))
            .dataSetFeatureMapping("srcInput")
            // .dataSetFeatureMapping("trgInput")
            .dataSetLabelMapping("label")
            // .trainEvaluation("out", 0, evaluation)
            .build();

        sd.setTrainingConfig(config);
        sd.setListeners(new ScoreListener(epochs / 10));

        // UIListener listener = UIListener.builder(new File("ui_demo.bin"))
        // .plotLosses(1)
        // .updateRatios(10)
        // .build();
        // sd.setListeners(listener);

        DataSet trainSet = new DataSet(srcData, srcLabel);
        trainSet.shuffle();
        DataSetIterator iterator = new ListDataSetIterator<DataSet>(trainSet.asList(), 200);
        History hist = sd.fit().train(iterator, epochs).exec();
    }

    public boolean predictTargetDomain(double[] features) {
        double[][] tmp = new double[][] {features};
        INDArray input = new NDArray(tmp);

        HashMap<String, INDArray> map = new HashMap<>();
        map.put("trgInput", input);
        Map<String, INDArray> o = sd.output(map, "trgOutput");
        INDArray result = o.getOrDefault("trgOutput", null);

        return result.getDouble(0) > result.getDouble(1);
    }

    public boolean predict(double[] features) {
        double[][] tmp = new double[][] {features};
        INDArray input = new NDArray(tmp);

        HashMap<String, INDArray> map = new HashMap<>();
        map.put("srcInput", input);
        Map<String, INDArray> o = sd.output(map, "srcOutput");
        INDArray result = o.getOrDefault("srcOutput", null);

        return result.getDouble(0) > result.getDouble(1);
    }

    SDVariable coralLoss(String name, SameDiff sd, SDVariable src, SDVariable trg, int n1, int n2, int d) {
        SDVariable oneVector1 = sd.one("ones1_" + name, DataType.DOUBLE, n1, 1);
        SDVariable oneVector2 = sd.one("ones2_" + name, DataType.DOUBLE, n2, 1);
        SDVariable srcMean = 
            sd.transpose(sd.transpose(oneVector1).mmul(src))
            .mmul(sd.transpose(oneVector1).mmul(src)).div(n1);
        SDVariable trgMean = 
            sd.transpose(sd.transpose(oneVector2).mmul(trg))
            .mmul(sd.transpose(oneVector2).mmul(trg)).div(n2);
        SDVariable srcCov = 
            sd.transpose(src).mmul(src).sub(srcMean).div(n1 - 1);
        SDVariable trgCov = 
            sd.transpose(trg).mmul(trg).sub(trgMean).div(n2 - 1);
        
        SDVariable loss = srcCov.sub(trgCov).norm2().div(d * 4);
        return loss;
    }
}
