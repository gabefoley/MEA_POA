package Alignment;

/**
 * Created by gabe on 6/09/2016.
 */

import org.rosuda.JRI.*;

public class RTest {
    public static void main(String a[]) {

//        // Create an R vector in the form of a string.
//        String javaVector = "c(1,2,3,4,5)";
//
//        // Start Rengine.
//        Rengine engine = new Rengine(new String[] { "--no-save" }, false, null);
//
//        // The vector that was created in JAVA context is stored in 'rVector' which is a variable in R context.
//        engine.eval("rVector=" + javaVector);
//
//        //Calculate MEAN of vector using R syntax.
//        engine.eval("meanVal=mean(rVector)");
//
//        //Retrieve MEAN value
//        double mean = engine.eval("meanVal").asDouble();
//
//        //Print output values
//        System.out.println("Mean of given vector is=" + mean);

        Rengine engine = new Rengine(new String[] { "--no-save"}, false, null);

        engine.eval("require(ggplot2)");
        engine.eval("Seqs <- c(100, 200, 300, 400)");
        engine.eval("Time <- c(2, 4, 16, 25)");
        engine.eval("joined <- data.frame(Seqs, Time)");
        engine.eval("plot <- ggplot(data=joined, aes(x=Seqs, y=Time, group=1)) + geom_line()");
        engine.eval("ggsave('/Users/gabe/Dropbox/java.png', plot)");
        System.out.println("Done");



    }
}
