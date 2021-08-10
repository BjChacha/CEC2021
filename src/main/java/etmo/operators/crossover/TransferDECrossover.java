package etmo.operators.crossover;

import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.encodings.solutionType.RealSolutionType;
import etmo.util.Configuration;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.wrapper.XReal;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class TransferDECrossover extends Crossover{
    /**
     * DEFAULT_CR defines a default CR (crossover operation control) value
     */
    private static final double DEFAULT_CR = 0.9;

    /**
     * DEFAULT_F defines the default F (Scaling factor for mutation) value
     */
    private static final double DEFAULT_F1_ = 0.5;
    private static final double DEFAULT_F2_ = 0.5;
    private static final double DEFAULT_F3_ = 0.5;

    /**
     * DEFAULT_K defines a default K value used in variants current-to-rand/1
     * and current-to-best/1
     */
    private static final double DEFAULT_K = 0.5;

    /**
     * DEFAULT_VARIANT defines the default DE variant
     */

    private static final String DEFAULT_DE_VARIANT = "rand/1/bin";

    /**
     * Valid solution types to apply this operator
     */
    private static final List VALID_TYPES = Arrays.asList(RealSolutionType.class);

    private double K_;
    private String DE_Variant_; // DE variant (rand/1/bin, rand/1/exp, etc.)

    private double[] CR_;
    private double[] F1_;
    private double[] F2_;
    private double[] F3_;

    /**
     * Constructor
     */
    public TransferDECrossover(HashMap<String, Object> parameters) {
        super(parameters);

        K_ = DEFAULT_K;
        DE_Variant_ = DEFAULT_DE_VARIANT;

        if (parameters.get("CR") != null)
            Arrays.fill(CR_, (Double) parameters.get("CR"));
        if (parameters.get("K") != null)
            K_ = (Double) parameters.get("K");
        if (parameters.get("DE_VARIANT") != null)
            DE_Variant_ = (String) parameters.get("DE_VARIANT");

    } // Constructor

    /**
     * Constructor
     */
    // public RandomDECrossover(Properties properties) {
    // this();
    // CR_ = (new Double((String)properties.getProperty("CR_")));
    // F_ = (new Double((String)properties.getProperty("F_")));
    // K_ = (new Double((String)properties.getProperty("K_")));
    // DE_Variant_ = properties.getProperty("DE_Variant_") ;
    // } // Constructor

    public void adaptive(SolutionSet target, SolutionSet source){
        // assume target <- source
        double[] stdTarget = target.getStd();
        double[] stdSource = source.getStd();
        double[] meanTarget = target.getMean();
        double[] meanSource = source.getMean();
        int len = meanSource.length;
        double[] deltaMean = new double[len];
        for (int i = 0; i < len; i++){
            deltaMean[i] = Math.abs(meanTarget[i] - meanSource[i]);
        }

        // adaptive
        CR_ = new double[len];
        F1_ = new double[len];
        F2_ = new double[len];
        for (int i = 0; i < len; i++) {

            F1_[i] = 2 * Math.sqrt(deltaMean[i]);
            F2_[i] = stdSource[i];
        }
    }

    /**
     * Executes the operation
     *
     * @param object
     *            An object containing an array of three parents
     * @return An object containing the offSprings
     */
    public Object execute(Object object) throws JMException {
        Solution[] parent = (Solution[]) object;
        Solution current = parent[0];

        Solution child;
        int jRand;

        child = new Solution(current);
        XReal x_target = new XReal(parent[0]);      // x_target
        XReal x_source = new XReal(parent[1]);      // x_source
        XReal x_source_r1 = new XReal(parent[2]);   // x_source_r1
        XReal x_source_r2 = new XReal(parent[3]);   // x_source_r2
        XReal x_target_r1 = new XReal(parent[4]);   // x_target_r1
        XReal x_target_r2 = new XReal(parent[5]);   // x_target_r2
        XReal xCurrent = new XReal(current);
        XReal xChild = new XReal(child);

        int numberOfVariables = x_target.getNumberOfDecisionVariables();
        jRand = PseudoRandom.randInt(0, numberOfVariables - 1);

        for (int j = 0; j < numberOfVariables; j++) {
            if (PseudoRandom.randDouble(0, 1) < CR_[j] || j == jRand) {
                double value;
                value = x_source_r1.getValue(j) + F1_[j] * (x_source.getValue(j) - x_target.getValue(j))
                                                + F2_[j] * (x_source_r1.getValue(j) - x_source_r2.getValue(j))
                                                + F3_[j] * (x_target_r1.getValue(j) - x_target_r2.getValue(j));
                if (value < xChild.getLowerBound(j))
                    value = xChild.getLowerBound(j);
                if (value > xChild.getUpperBound(j))
                    value = xChild.getUpperBound(j);
                xChild.setValue(j, value);
            } else {
                double value;
                value = xCurrent.getValue(j);
                xChild.setValue(j, value);
            } // else
        } // for

        return child;
    }
}
