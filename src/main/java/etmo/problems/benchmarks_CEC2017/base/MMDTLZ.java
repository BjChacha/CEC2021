package etmo.problems.benchmarks_CEC2017.base;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.problems.base.staticBase.GFunctions;
import etmo.util.JMException;
import org.uma.jmetal.solution.doublesolution.DoubleSolution;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.SocketHandler;

/**
 * @Author: Zhi-Ming Dong, dzm.neu@gmail.com
 * @Date: created in 19-1-15 15:11
 * @Version: v
 * @Descriptiom: #
 * 1#
 * @Modified by:
 */
public class MMDTLZ extends Problem {
    String gType_;
    String genType_;
    private int alpha;

    public MMDTLZ(int numberOfObjectives, int numberOfVariables, int alpha, double lg, double ug) {
        numberOfObjectives_ = numberOfObjectives;
        numberOfVariables_ = numberOfVariables;
        
        if (numberOfObjectives == 2)
            hType_ = "circle";
        else
            hType_ = "sphere";
        gType_ = "sphere";
        genType_ = "multiplication";
        
        int num = numberOfVariables - numberOfObjectives + 1;

        shiftValues_ = new double[num];
        rotationMatrix_ = new double[num][num];

        upperLimit_ = new double[numberOfVariables_];
        lowerLimit_ = new double[numberOfVariables_];

        for (int var = 0; var < numberOfObjectives_ - 1; var++) {
            lowerLimit_[var] = 0.0;
            upperLimit_[var] = 1.0;
        } // for

        for (int var = numberOfObjectives_ - 1; var < numberOfVariables; var++) {
            lowerLimit_[var] = lg;
            upperLimit_[var] = ug;
        }

        for (int i = 0; i < num; i++)
            shiftValues_[i] = 0;

        for (int i = 0; i < num; i++) {
            for (int j = 0; j < num; j++) {
                if (i != j)
                    rotationMatrix_[i][j] = 0;
                else
                    rotationMatrix_[i][j] = 1;
            }
        }
    }
    
    public void evaluate(Solution solution) throws JMException {
        double vars[] = scaleVariables(solution);

        double[] xI = new double[numberOfObjectives_ - 1];
        double[] xII = new double[numberOfVariables_ - numberOfObjectives_ + 1];

        for (int i = 0; i < numberOfObjectives_ - 1; i++)
            xI[i] = vars[i];

        for (int i = numberOfObjectives_ - 1; i < numberOfVariables_; i++)
            xII[i - numberOfObjectives_ + 1] = vars[i];

        double[] f = new double[numberOfObjectives_];

        double g = evalG(xII);

        for (int i = 0; i < numberOfObjectives_; i++) {
            f[i] = 1 + g;
        }

//        solution.setGFunValue(1 + g);

        for (int i = 0; i < numberOfObjectives_; i++) {
            for (int j = 0; j < numberOfObjectives_ - (i + 1); j++) {
                f[i] *= Math.cos(Math.pow(xI[j], alpha) * 0.5 * Math.PI);
            }
            if (i != 0) {
                int aux = numberOfObjectives_ - (i + 1);
                f[i] *= Math.sin(Math.pow(xI[aux], alpha) * 0.5 * Math.PI);
            } // if
        } // for

        for (int i = 0; i < getNumberOfObjectives(); i++) {
            solution.setObjective(startObjPos_ + i, f[i]);
        }
    }

    double[] evalH(double[] xI, double g) {
        double[] h = new double[numberOfObjectives_];
        if (hType_.equalsIgnoreCase("lineoid")){
            for (int i = 0; i < numberOfObjectives_; i++) {
                h[i] = 1.0;
                for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
                    h[i] *= xI[j];
                if (i != 0) {
                    int aux = numberOfObjectives_ - (i + 1);
                    h[i] *= (1 - xI[aux]);
                } // if
            } // for
        }else if(hType_.equalsIgnoreCase("circle" ) || hType_.equalsIgnoreCase("sphere")){
            for (int i = 0; i < numberOfObjectives_; i++) {
                h[i] = 1.0;
                for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
                    h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
                if (i != 0) {
                    int aux = numberOfObjectives_ - (i + 1);
                    h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
                } // if
            } // for
        }else if(hType_.equalsIgnoreCase("convex")){
            for (int i = 0; i < numberOfObjectives_; i++) {
                h[i] = 1.0;
                for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
                    h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
                if (i != 0) {
                    int aux = numberOfObjectives_ - (i + 1);
                    h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
                } // if
                if(i != numberOfObjectives_-1)
                    h[i] = Math.pow(h[i], 4);
                else
                    h[i] = Math.pow(h[i], 2);
            } // for
        }else if(hType_.equalsIgnoreCase("degenerate")){
            degenerateEvalPV(xI, g);
            for (int i = 0; i < numberOfObjectives_; i++) {
                h[i] = 1.0;
                for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
                    h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
                if (i != 0) {
                    int aux = numberOfObjectives_ - (i + 1);
                    h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
                } // if
            } // for
        }else if(hType_.equalsIgnoreCase("irconcave")){
            shiftEvalPV(xI);
            for (int i = 0; i < numberOfObjectives_; i++) {
                h[i] = 1.0;
                for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
                    h[i] *= Math.cos(xI[j] * 0.5 * Math.PI);
                if (i != 0) {
                    int aux = numberOfObjectives_ - (i + 1);
                    h[i] *= Math.sin(xI[aux] * 0.5 * Math.PI);
                } // if
            } // for
        }else if(hType_.equalsIgnoreCase("disconnect")){
            h[numberOfObjectives_-1] = 0;
            for (int i = 0; i < numberOfObjectives_-1; i++) {
                if(genType_.equalsIgnoreCase("multiplication"))
                    h[i] = xI[i]/(1+g);
                else if(genType_.equalsIgnoreCase("addition"))
                    h[i] = xI[i] - g;
                else{
                    System.out.println("Error: Generation type " + genType_ + " invalid");
                    System.exit(0);
                }
                h[numberOfObjectives_-1] += (xI[i]*(1+Math.sin(3*Math.PI*xI[i])))/(1+g);
            } // for
            h[numberOfObjectives_-1] = numberOfObjectives_ - h[numberOfObjectives_-1];
        }else {
            System.out.println("Error: H function type " + hType_ + " invalid");
            System.exit(0);
        }
        return h;
    }

    double evalG(double[] xII) throws JMException {
        if(gType_.equalsIgnoreCase("F1") || gType_.equalsIgnoreCase("sphere"))
            return etmo.problems.base.staticBase.GFunctions.getF1(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F2"))
            return etmo.problems.base.staticBase.GFunctions.getF2(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F3"))
            return etmo.problems.base.staticBase.GFunctions.getF3(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F4"))
            return etmo.problems.base.staticBase.GFunctions.getF4(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F5") || gType_.equalsIgnoreCase("rosenbrock"))
            return etmo.problems.base.staticBase.GFunctions.getF5(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F6") || gType_.equalsIgnoreCase("ackley"))
            return etmo.problems.base.staticBase.GFunctions.getF6(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F7"))
            return etmo.problems.base.staticBase.GFunctions.getF7(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F8") || gType_.equalsIgnoreCase("griewank"))
            return etmo.problems.base.staticBase.GFunctions.getF8(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F9") || gType_.equalsIgnoreCase("rastrigin"))
            return etmo.problems.base.staticBase.GFunctions.getF9(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F10"))
            return etmo.problems.base.staticBase.GFunctions.getF10(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F11"))
            return etmo.problems.base.staticBase.GFunctions.getF11(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F12"))
            return etmo.problems.base.staticBase.GFunctions.getF12(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F13"))
            return etmo.problems.base.staticBase.GFunctions.getF13(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("mean"))
            return etmo.problems.base.staticBase.GFunctions.getMean(xII);
        else if (gType_.equalsIgnoreCase("F17"))
            return etmo.problems.base.staticBase.GFunctions.getF17(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F18"))
            return etmo.problems.base.staticBase.GFunctions.getF18(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("F19"))
            return etmo.problems.base.staticBase.GFunctions.getF19(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("HF1"))
            return etmo.problems.base.staticBase.GFunctions.getHF1(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("HF2"))
            return etmo.problems.base.staticBase.GFunctions.getHF2(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("HF3"))
            return etmo.problems.base.staticBase.GFunctions.getHF3(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("HF4"))
            return etmo.problems.base.staticBase.GFunctions.getHF4(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("HF5"))
            return etmo.problems.base.staticBase.GFunctions.getHF5(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("HF6"))
            return etmo.problems.base.staticBase.GFunctions.getHF6(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("HF7"))
            return etmo.problems.base.staticBase.GFunctions.getHF7(xII, shiftValues_, rotationMatrix_);
        else if (gType_.equalsIgnoreCase("HF8"))
            return GFunctions.getHF8(xII, shiftValues_, rotationMatrix_);
        else {
            System.out.println("Error: g function type " + gType_ + " invalid");
            return Double.NaN;
        }
    }

    protected static void degenerateEvalPV(double[] xI, double g){
        int I = 2;
        for(int i=I-1;i<xI.length;i++){
            xI[i] = 0.5*(1.0+2*g*xI[i])/(1.0+g);
        }
    }

    public void setGType(String gType) {
        gType_ = gType;
    }

    public void setGenType(String genType) {
        genType_ = genType;
    }

    public String getHType() {
        return hType_;
    }

    protected static void shiftEvalPV(double[] xI){
        for(int i=0; i<xI.length; i++){
            xI[i] = 0.5*xI[i] + 0.25;
        }
    }

    @Override
    public void dynamicEvaluate(Solution solution, int currentGeneration) throws JMException {

    }
}
