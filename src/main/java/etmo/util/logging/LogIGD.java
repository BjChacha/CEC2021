package etmo.util.logging;

import etmo.util.Configuration;

import java.io.*;
import java.util.Arrays;

public class LogIGD {
    public static void LogIGD(String algoName, int problemNum, double[][] igds){
        String folderPath = "./data/IGDs/_collections/" + algoName;
        File folder = new File(folderPath);
        if (!folder.exists()){
            folder.mkdirs();
        }

        int taskNum = igds.length;
        int times = igds[0].length;
        String filePath = folderPath + "/" + algoName + "_ETMOF" + problemNum + "_" + times + ".txt";
        try {
            /* Open the file */
            FileOutputStream fos = new FileOutputStream(filePath);
            OutputStreamWriter osw = new OutputStreamWriter(fos);
            BufferedWriter bw = new BufferedWriter(osw);

            for (double[] igd : igds) {
                // if (this.vector[i].getFitness()<1.0) {
                bw.write(Arrays.toString(igd));
                bw.newLine();
                // }
            }

            /* Close the file */
            bw.close();
        } catch (IOException e) {
            Configuration.logger_.severe("Error acceding to the file");
            e.printStackTrace();
        }
    }

    public static void MarkLog(String AlgoName, String problemName, double[][] igds){
        String folderPath = "./data/IGDs/_study/" + AlgoName;
        File folder = new File(folderPath);
        if (!folder.exists()){
            folder.mkdirs();
        }

        int taskNum = igds.length;
        String filePath = folderPath + "/" + AlgoName + "_" + problemName + ".txt";
        try {
            /* Open the file */
            FileOutputStream fos = new FileOutputStream(filePath);
            OutputStreamWriter osw = new OutputStreamWriter(fos);
            BufferedWriter bw = new BufferedWriter(osw);

            for (double[] igd : igds) {
                // if (this.vector[i].getFitness()<1.0) {
                bw.write(Arrays.toString(igd));
                bw.newLine();
                // }
            }

            /* Close the file */
            bw.close();
        } catch (IOException e) {
            Configuration.logger_.severe("Error acceding to the file");
            e.printStackTrace();
        }
    }
}
