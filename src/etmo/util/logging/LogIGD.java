package etmo.util.logging;

import etmo.core.Solution;
import etmo.util.Configuration;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

public class LogIGD {
    public static void LogIGD(String fileName, double[] igds){
        try {
            /* Open the file */
            fileName = "data//" + fileName;
            FileOutputStream fos = new FileOutputStream(fileName);
            OutputStreamWriter osw = new OutputStreamWriter(fos);
            BufferedWriter bw = new BufferedWriter(osw);

            for (double igd : igds) {
                // if (this.vector[i].getFitness()<1.0) {
                bw.write(Double.toString(igd));
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
