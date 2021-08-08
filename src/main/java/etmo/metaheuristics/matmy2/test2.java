package etmo.metaheuristics.matmy2;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.problems.benchmarks_CEC2021.*;
import etmo.util.JMException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class test2 {
    public static void main(String[] args) throws IOException, JMException, ClassNotFoundException {
        System.out.println("Populations optimizing...");

        ProblemSet problemSet = ETMOF28.getProblem();
//        MtoAlgorithm algorithm = getOptimizer(problemSet);
//        SolutionSet[] population = algorithm.execute();

        System.out.println("Populations optimization completed!");
        // -------------------------EXPERIMENTS--------------------
        System.out.println("Start my experiments...");


        double[][] mat1 = readMatFromFile("src/etmo/metaheuristics/matmy2/mat/mat1.txt");
        double[][] mat2 = readMatFromFile("src/etmo/metaheuristics/matmy2/mat/mat2.txt");;

        int[][] taskPairs = {
                {0, 3},  // good -> good
                {0, 1},  // good -> bad
                {1, 0},  // bad -> good
        };

        for (int[] pair: taskPairs) {
            int srcTask = pair[0];
            int trgTask = pair[1];

//            int d1 = population[srcTask].getMat()[0].length;
//            int d2 = population[trgTask].getMat()[0].length;

            int d1 = mat1[0].length;
            int d2 = mat2[0].length;

            int[] HGs = {(d1+d2)/5, (d1+d2)/2, (d1+d2)};
            int[] HDs = {(d1+d2)/5, (d1+d2)/2, (d1+d2)};
            double[] LRGs = {1e-1, 5e-2, 1e-2, 5e-3, 1e-3};
            double[] LRDs = {1e-1, 5e-2, 1e-2, 5e-3, 1e-3};
            int[] EPOCHs = {3, 5, 8};

            for (int HG : HGs) {
                for (int HD : HDs) {
                    for (double LRG: LRGs) {
                        for (double LRD : LRDs) {
                            for (int e: EPOCHs) {
                                for (int c = 0; c < 3; c++) {
//                                System.out.println("---------------Schedue------------------");
//                                System.out.println("Transfer: " + srcTask + " -> " + trgTask);
//                                System.out.println("hG: " + HG + "\thD: " + HD + "\tLRG: "+LRG + "\tLRD: "+LRD+ "\trun time: " + (c + 1));

//                                double[][] srcMat = population[srcTask].getMat();
//                                double[][] trgMat = population[trgTask].getMat();

                                    double[][] srcMat = mat1;
                                    double[][] trgMat = mat2;

//                                for (int i = 0; i < srcMat.length; i++) {
//                                    System.out.println(Arrays.toString(srcMat[i]));
//                                }
//                                System.out.println(" ");
//                                for (int i = 0; i < srcMat.length; i++) {
//                                    System.out.println(Arrays.toString(trgMat[i]));
//                                }

                                    double[][] transferVariables = Utils.GANTest(
                                            srcMat,
                                            trgMat,
                                            HG, HD, LRG, LRD, e
                                    );

                                    int betterCnt = 0;
                                    double offset = 0;
                                    for (int i = 0; i < trgMat.length; i++) {
                                        for (int j = 0; j < trgMat[0].length; j++) {
                                            offset += Math.pow((trgMat[i][j] - transferVariables[i][j]), 2);
                                        }
//                                    Solution transferIndividual = new Solution(problemSet);
//                                    transferIndividual.setDecisionVariables(transferVariables[i]);
//                                    problemSet.get(trgTask).evaluate(transferIndividual);
//                                    System.out.println("Pause");
                                    }
                                    offset = offset / trgMat.length;
                                    if (offset < 0.35) {
                                        System.out.println("hG: " + HG + "\thD: " + HD + "\tLRG: " + LRG + "\tLRD: " + LRD + "\tEpoch: "+e+ "\trun time: " + (c + 1));
                                        System.out.println("offset: " + offset);
                                    }
                                }
                            }
                        }
//                    // compare
//                    double[] trgAveObjs = population[trgTask].getAverageObjectiveVector();
//                    double[] trObjs = transferIndividual.getObjectives();
//                    int objStart = problemSet.get(trgTask).getStartObjPos();
//                    int objCnt = problemSet.get(trgTask).getNumberOfObjectives();
//                    double[] tmp1 = new double[objCnt];
//                    double[] tmp2 = new double[objCnt];
//                    for (int k = 0; k < objCnt; k++) {
//                        tmp1[k] = trgAveObjs[k+objStart];
//                        tmp2[k] = trObjs[k+objStart];
//                    }
//                    trgAveObjs = tmp1;
//                    trObjs = tmp2;
//                    System.out.println("Target population average objectives: "
//                            +Arrays.toString(trgAveObjs));
//                    System.out.println("Transfer individual objectives: "
//                            + Arrays.toString(trObjs));

//                    boolean better = true;
//                    for (int j = 0; j < trgAveObjs.length; j++){
//                        if (trgAveObjs[j] < trObjs[j]){
//                            better = false;
//                            break;
//                        }
//                    }
//                    if (better){
//                        betterCnt++;
//                    }
//                        System.out.println("Better Transfer Individuals Count: "+betterCnt);
                    }
                }
            }
        }
    }

    static MtoAlgorithm getOptimizer(ProblemSet problemSet) throws JMException {
        MtoAlgorithm algorithm = new MaTMY2(problemSet);

        algorithm.setInputParameter("populationSize", 100);
        algorithm.setInputParameter("maxEvaluations", 1000 * 100 * problemSet.size());
        algorithm.setInputParameter("transferVolume", 10);
        algorithm.setInputParameter("baseRunTime", 3);

        algorithm.setInputParameter("forceTransferRate", 0.2);

        algorithm.setInputParameter("scoreIncrement", 0.5);
        algorithm.setInputParameter("scoreDecreaseRate", 0.2);
        algorithm.setInputParameter("isDRA", false);
        algorithm.setInputParameter("algoName", "MaOEAC");

        HashMap parameters;
        parameters = new HashMap();
        parameters.put("probability", 0.9);
        parameters.put("distributionIndex", 20.0);
        Operator crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

        Operator selection= SelectionFactory.getSelectionOperator("BinaryTournament2", null);

        algorithm.addOperator("crossover",crossover);
        algorithm.addOperator("selection", selection);

        return algorithm;
    }

    static double[][] readMatFromFile(String path) throws IOException {
        FileReader fr = new FileReader(path);
        BufferedReader br = new BufferedReader(fr);
        String line = null;
        double[][] array = new double[99][];
        double[] lineArray;
        String[] sp;
        int ind = 0;
        while ((line = br.readLine()) != null){
            sp = line.split(",");
            lineArray = new double[sp.length];
            for (int j = 0; j < sp.length; j++){
                lineArray[j] = Double.parseDouble(sp[j]);
            }
            array[ind++] = lineArray;
        }
        return array;
    }
}
