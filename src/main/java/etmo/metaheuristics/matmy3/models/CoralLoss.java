package etmo.metaheuristics.matmy3.models;

import org.nd4j.autodiff.samediff.SDVariable;
import org.nd4j.autodiff.samediff.SameDiff;
import org.nd4j.common.primitives.Pair;
import org.nd4j.linalg.activations.IActivation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.lossfunctions.ILossFunction;
import org.nd4j.linalg.lossfunctions.SameDiffLoss;

public class CoralLoss extends SameDiffLoss {

    @Override
    public SDVariable defineLoss(SameDiff sd, SDVariable layerInput, SDVariable labels) {
        // TODO Auto-generated method stub
        return null;
    }
}
