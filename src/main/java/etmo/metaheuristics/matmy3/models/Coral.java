package etmo.metaheuristics.matmy3.models;

import org.nd4j.autodiff.samediff.SDVariable;
import org.nd4j.autodiff.samediff.SameDiff;
import org.nd4j.linalg.api.buffer.DataType;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.weightinit.impl.XavierInitScheme;

public class Coral {
    public Coral() {
    
    }

    public void train(double[][] sourceDomainData, double[][] targetDomainData) {
        int inD = sourceDomainData[0].length;
        int hD1 = inD / 2;
        int hD2 = hD1 / 2;
        int outD = 1;

        SameDiff sd = SameDiff.create();

        INDArray srcData = new NDArray(sourceDomainData);
        INDArray trgData = new NDArray(targetDomainData);

        SDVariable srcInput = sd.placeHolder("srcInput", DataType.DOUBLE, -1, inD);
        SDVariable trgInput = sd.placeHolder("trgInput", DataType.DOUBLE, -1, inD);
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

        // compute flow
        SDVariable ws1 = sd.nn().tanh(srcInput.mmul(Wsrc1).add(Bsrc1));
        SDVariable ws2 = sd.nn().tanh(ws1.mmul(Wsrc2).add(Bsrc2));
        SDVariable output = sd.nn().sigmoid(ws2.mmul(Wsrc3).add(Bsrc3));

        SDVariable wt1 = sd.nn().tanh(trgInput.mmul(Wtrg1).add(Btrg1));
        SDVariable wt2 = sd.nn().tanh(wt1.mmul(Wtrg2).add(Btrg2));

        // TODO: compute loss
    }

    // TODO: Coral loss function
}
