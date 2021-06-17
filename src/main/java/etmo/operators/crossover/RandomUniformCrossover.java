package etmo.operators.crossover;

import etmo.core.Solution;

import etmo.encodings.solutionType.RealSolutionType;
import etmo.util.Configuration;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.wrapper.XReal;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class RandomUniformCrossover extends Crossover{

    private static final double DEFAULT_CR_UB = 0.9;
    private static final double DEFAULT_CR_LB = 0.2;

    private static final String DEFAULT_DE_VARIANT = "rand/1/bin";

    /**
     * Valid solution types to apply this operator
     */
    private static final List VALID_TYPES = Arrays.asList(RealSolutionType.class);

    private double CR_UB_;
    private double CR_LB_;

    /**
     * Constructor
     */
    public RandomUniformCrossover(HashMap<String, Object> parameters) {
        super(parameters);

        CR_UB_ = DEFAULT_CR_UB;
        CR_LB_ = DEFAULT_CR_LB;

        if (parameters.get("CR_UB") != null)
            CR_UB_ = (Double) parameters.get("CR_UB");
        if (parameters.get("CR_LB") != null)
            CR_LB_ = (Double) parameters.get("CR_LB");

    } // Constructor

    /**
     * Constructor
     */
    // public RandomUniformCrossover(Properties properties) {
    // this();
    // CR_ = (new Double((String)properties.getProperty("CR_")));
    // F_ = (new Double((String)properties.getProperty("F_")));
    // K_ = (new Double((String)properties.getProperty("K_")));
    // DE_Variant_ = properties.getProperty("DE_Variant_") ;
    // } // Constructor

    /**
     * Executes the operation
     *
     * @param object
     *            An object containing an array of three parents
     * @return An object containing the offSprings
     */
    public Object execute(Object object) throws JMException {
        Solution[] parents = (Solution[]) object;

        if (parents.length != 2) {
            Configuration.logger_.severe("RandomUniformCrossover.execute: operator needs two " + "parents");
            Class cls = java.lang.String.class;
            String name = cls.getName();
            throw new JMException("Exception in " + name + ".execute()");
        } // if

        if (!(VALID_TYPES.contains(parents[0].getType().getClass())
                && VALID_TYPES.contains(parents[1].getType().getClass()))) {
            Configuration.logger_.severe("RandomUniformCrossover.execute: the solutions " + "type " + parents[0].getType()
                    + " is not allowed with this operator");

            Class cls = java.lang.String.class;
            String name = cls.getName();
            throw new JMException("Exception in " + name + ".execute()");
        } // if

        Solution[] offSpring;
        double crossoverProbability_ = PseudoRandom.randDouble(CR_LB_, CR_UB_);
        offSpring = doCrossover(crossoverProbability_, parents[0], parents[1]);
        return offSpring;
    }

    private Solution[] doCrossover(double probability, Solution parent1, Solution parent2) throws JMException {
        Solution[] offSpring = new Solution[2];

        offSpring[0] = new Solution(parent1);
        offSpring[1] = new Solution(parent2);
        XReal x1 = new XReal(parent1);
        XReal x2 = new XReal(parent2);
        XReal offs1 = new XReal(offSpring[0]);
        XReal offs2 = new XReal(offSpring[1]);

        double valueX1, valueX2;

        int maxVar = x1.getNumberOfDecisionVariables();

        int k = PseudoRandom.randInt(0, maxVar - 1);

        for (int i = 0; i < maxVar; i++) {
            valueX1 = x1.getValue(i);
            valueX2 = x2.getValue(i);
            if (i == k || PseudoRandom.randDouble() < probability){
                offs1.setValue(i, valueX2);
                offs2.setValue(i, valueX1);
            }
            else{
                offs1.setValue(i, valueX1);
                offs2.setValue(i, valueX2);
            }
        }
        return offSpring;
    }
}
