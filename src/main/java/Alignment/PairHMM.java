package Alignment;

import SubstitutionModels.ExampleModel;
import com.sun.xml.internal.bind.v2.TODO;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by gabe on 3/05/2016.
 */
public class PairHMM {

    private double transitionMM, transitionMX, transitionMY, transitionXX, transitionYY, transitionXM, transitionYM;
    private String seq1, seq2;
    private double[][] vM;
    private double[][] vX;
    private double[][] vY;
    private String[][] tracebackM, tracebackX, tracebackY;

    private HashMap<Integer, Integer> alignedPairs;

    private int longerSide;
    private int shorterSide;
    private String longerSeq;
    private String shorterSeq;

    private double tau;


    public PairHMM(String seq1, String seq2, double tau, double epsilon, double delta) {

        this.seq1 = seq1;
        this.seq2 = seq2;

        this.transitionMM = 1 - (2 * delta) - tau;
        this.transitionMX = this.transitionMY = delta;
        this.transitionXX = this.transitionYY = epsilon;
        this.transitionXM = this.transitionYM = 1 - epsilon - tau;

        this.tau = tau;

        // Probability of aligning A,G,C,T with a gap


        this.vM = new double[seq1.length() + 1][seq2.length() + 1];
        this.vX = new double[seq1.length() + 1][seq2.length() + 1];
        this.vY = new double[seq1.length() + 1][seq2.length() + 1];

        this.tracebackM = new String[seq1.length() + 1][seq2.length() + 1];
        this.tracebackX = new String[seq1.length() + 1][seq2.length() + 1];
        this.tracebackY = new String[seq1.length() + 1][seq2.length() + 1];

        this.alignedPairs = new HashMap<Integer, Integer>();


        vM[0][0] = 1;

        double[][] mat1 = new double[2][2];
        double[][] mat2 = new double[2][2];





//        double[][] mat3 = mat1 * mat2


//        if (vM.length > vM[0].length) {
//
//            this.longerSide = vM.length;
//            this.shorterSide = vM[0].length;
//            this.longerSeq = seq1;
//            this.shorterSeq = seq2;
//
//        } else {
//            this.longerSide = vM[0].length;
//            this.shorterSide = vM.length;
//            this.longerSeq = seq2;
//            this.shorterSeq = seq1;
//        }

        for (int i = 1; i <= seq1.length() ; i++) {
            for (int j = 1; j <= seq2.length(); j ++) {
                fillVM(i, j, vM, vX, vY, tracebackM);
                fillVX(i, j, vM, vX, tracebackX);
                fillVY(i, j, vM, vY, tracebackY);
            }

        }

//        printMatrix(tracebackM);
//        printMatrix(tracebackX);
//        printMatrix(tracebackY);

        traceback();

    }


        // Fill out the match matrix

    public void fillVM(int i, int j, double[][] vM, double[][] vX, double[][] vY, String[][] tracebackM) {

//        printMatrix(tracebackM);
//        for (int k = 1; k < j; k++) {


            double emissionM = ExampleModel.getDistance(seq1.charAt(i-1), seq2.charAt(j - 1));

            // Get the actual costs for transitioning for each of the states
            double currentTransitionMM = transitionMM * vM[i - 1][j - 1];
            double currentTransitionXM = transitionXM * vX[i - 1][j - 1];
            double currentTransitionYM = transitionYM * vY[i - 1][j - 1];

//            System.out.println(i + " " + k);
//            System.out.println(currentTransitionMM * emissionM + " " + currentTransitionXM * emissionM + " " + currentTransitionYM * emissionM);



            // Work out the optimal cost and set the cell
            if (currentTransitionMM > currentTransitionXM && currentTransitionMM > currentTransitionYM) {
                vM[i][j] = currentTransitionMM * emissionM;
                tracebackM[i][j] = "M";
            } else if (currentTransitionXM > currentTransitionMM && currentTransitionXM > currentTransitionYM) {
                vM[i][j] = currentTransitionXM * emissionM;
                tracebackM[i][j] = "X";
            } else {
                vM[i][j] = currentTransitionYM * emissionM;
                tracebackM[i][j] = "Y";

            }


        }


//    }


    // Fill out the gap in X matrix
    public void fillVX(int i, int j, double[][] vM, double[][] vX, String[][] tracebackX) {
        double emissionX = 0.25;


//        for (int k = 1; k < j; k++) {


            double currentTransitionXX = transitionXX * vX[i][j - 1];
            double currentTransitionMX = transitionMX * vM[i][j - 1];


            if (currentTransitionXX > currentTransitionMX) {
                vX[i][j] = currentTransitionXX * emissionX;
                tracebackX[i][j] = "X";

            } else {
                vX[i][j] = currentTransitionMX * emissionX;
                tracebackX[i][j] = "M";

            }

//        }

    }


    // Fill out the gap in Y matrix
    public  void fillVY(int i, int j, double[][] vM, double[][] vY, String[][] tracebackY) {

        double emissionY = 0.25;
//        for (int k = 1; k < j; k++) {

            double currentTransitionYY = transitionYY * vY[i - 1][j];
            double currentTransitionMY = transitionMY * vM[i - 1][j];

            if (currentTransitionYY > currentTransitionMY) {
                vY[i][j] = currentTransitionYY * emissionY;
                tracebackY[i][j] = "Y";

            } else {
                vY[i][j] = currentTransitionMY * emissionY;
                tracebackY[i][j] = "M";

            }

//        }

    }

    public void traceback(){

//        printAllMatrices();

        String seq1Output = "";
        String seq2Output = "";
        String lastState;
        int i = seq1.length();
        int j = seq2.length();
        int n = longerSide - 1;


        while ((i > 0) && (j > 0)) {
            if ((vM[i][j] > vX[i][j]) && (vM[i][j] > vY[i][j])) {
                seq1Output = seq1.substring(i-1, i) + seq1Output;
                seq2Output = seq2.substring(j-1, j) + seq2Output;
                alignedPairs.put(i,j);
                lastState = tracebackM[i][j];
                i--;
                j--;

            }

            else if((vX[i][j]) > vM[i][j] && (vX[i][j]) > vY[i][j]){
                seq1Output = "-" + seq1Output;
                seq2Output = seq2.substring(j-1, j) + seq2Output;

                lastState = tracebackX[i][j];
                j--;
            }

            else {
                seq1Output = seq1.substring(i-1, i) + seq1Output;
                seq2Output = "-" + seq2Output;
                lastState = tracebackY[i][j];
                i--;
            }

            n--;

            while ((i > 0) && (j > 0)) {
                if (lastState == "M"){
                    seq1Output = seq1.charAt(i-1) + seq1Output;
                    seq2Output = seq2.charAt(j-1) + seq2Output;
                    lastState = tracebackM[i][j];
                    alignedPairs.put(i,j);
                    i--;
                    j--;
                }
                else if (lastState == "Y"){
                    seq1Output = seq1.charAt(i-1) + seq1Output;
                    seq2Output = "-" + seq2Output;
                    lastState = tracebackY[i][j];
                    i--;
                }

                else {
                    seq1Output = "-" + seq1Output;
                    seq2Output = seq2.charAt(j-1) + seq2Output;
                    lastState = tracebackX[i][j];
                    j--;
                }
                n--;
            }
            n++;
//            System.out.println(n);


            System.out.println(seq1Output);
            System.out.println(seq2Output);
        }

    }

    public double[][] forwardAlgorithm() {

        double[][] fM = new double[seq1.length() + 1][seq2.length() + 1];
        double[][] fX = new double[seq1.length() + 1][seq2.length() + 1];
        double[][] fY = new double[seq1.length() + 1][seq2.length() + 1];

        fM[0][0] = 1;

        for (int i = 0; i <= seq1.length(); i++) {
            for (int j = 0; j <= seq2.length(); j++) {
                if (i!= 0 || j!=0) {
                    sumfM(i, j, fM, fX, fY);
                    sumfX(i, j, fM, fX);
                    sumfY(i, j, fM, fY);
                }
            }

        }

        double forwardProb = tau * (fM[seq1.length()][seq2.length()] + fX[seq1.length()][seq2.length()] +
                fY[seq1.length()][seq2.length()]);

        return fM;

//        System.out.println("Forward probability is " + forwardProb);



    }

    public void sumfM(int i, int j, double[][] fM, double[][]fX, double[][] fY){

            if((i - 1 < 0) || (j - 1 < 0)){
                fM[i][j] = 0;
            }

            else {

                double emissionM = ExampleModel.getDistance(seq1.charAt(i - 1), seq2.charAt(j - 1));
                double forwardMM = transitionMM * fM[i - 1][j - 1];
                double forwardXM = transitionXM * fX[i - 1][j - 1];
                double forwardYM = transitionYM * fY[i - 1][j - 1];

                fM[i][j] = emissionM * (forwardMM + forwardXM + forwardYM);
            }



    }

    public void sumfX(int i, int j, double[][] fM, double[][]fX){
        double emissionX = 0.25;


            // If we're in the first column
            if (j - 1 < 0) {
                fX[i][j] = 0;
            } else {


                double forwardXX = transitionXX * fX[i][j - 1];
                double forwardMX = transitionMX * fM[i][j - 1];

                fX[i][j] = emissionX * (forwardMX + forwardXX);

            }

    }

    public void sumfY(int i, int j, double[][] fM, double[][] fY){

        double emissionY = 0.25;



            if (i - 1 < 0) {
                fY[i][j] = 0;
            } else {


                double forwardYY = transitionYY * fY[i - 1][j];
                double forwardMY = transitionMY * fM[i - 1][j];

                fY[i][j] = emissionY * (forwardMY + forwardYY);

            }


    }

    public  double[][] backwardAlgorithm() {

        double[][] bM = new double[seq1.length() + 2][seq2.length() + 2];
        double[][] bX = new double[seq1.length() + 2][seq2.length() + 2];
        double[][] bY = new double[seq1.length() + 2][seq2.length() + 2];

        bM[seq1.length()][seq2.length()] = tau;
        bX[seq1.length()][seq2.length()] = tau;
        bY[seq1.length()][seq2.length()] = tau;


        for (int i = seq1.length(); i > 0; i--) {
            for (int j = seq2.length(); j > 0 ; j--) {
                if (i!= seq1.length() || j!= seq2.length()) {
                    sumbM(i, j, bM, bX, bY);
                    sumbX(i, j, bM, bX);
                    sumbY(i, j, bM, bY);
                }
            }

        }

//        printMatrix(bM);
//        printMatrix(bX);
//        printMatrix(bY);

        double backwardProb = bM[1][1] + bX[1][1] + bY[1][1];

        System.out.println("Backward probability is " + backwardProb);

        return bM;


    }

    public void sumbM(int i, int j, double[][] bM, double[][] bX, double[][] bY){
//        System.out.println("i = " + i + " j = " + j );
//        System.out.println("seq1length = " + seq1.length() + " seq2length = " + seq2.length());
//        System.out.println("seq1 = " + seq1 + " seq2= " + seq2);
//        System.out.println(seq1.charAt(i-1));
//        System.out.println(seq1.charAt(i));

//        System.out.println(" BM: Cell is " + i + " " + j + " seq1char is " + seq1.charAt(i - 1) + " seq2char is " + seq2.charAt(j -1));

        double emissionM;
        if (i > seq1.length() - 1 || j > seq2.length() - 1){
            emissionM = 0;
        }
        else {
            emissionM = ExampleModel.getDistance(seq1.charAt(i), seq2.charAt(j));
        }        double emissionX = 0.25;
        double emissionY = 0.25;

//        if((i + 1 > seq1.length()) || (j + 1 > seq2.length())){
//            bM[i][j] = 0;
//        }
//        else {


            double backwardMM = emissionM * transitionMM * bM[i + 1][j + 1];
            double backwardXM = emissionX * transitionMX * bX[i][j + 1];
            double backwardYM = emissionY * transitionMY * bY[i + 1][j];

            bM[i][j] = backwardMM + backwardXM + backwardYM;
//        }

    }

    public void sumbX(int i, int j, double[][] bM, double[][] bX){

//        System.out.println(" BX: Cell is " + i + " " + j + " seq1char is " + seq1.charAt(i - 1) + " seq2char is " + seq2.charAt(j-1));


//        if (j + 1 > seq2.length()) {
//            bX[i][j] = 0;
//        } else {


        double emissionM;
        if (i > seq1.length() - 1 || j > seq2.length() - 1){
            emissionM = 0;
        }
        else {
            emissionM = ExampleModel.getDistance(seq1.charAt(i), seq2.charAt(j));
        }            double emissionX = 0.25;


            double backwardMX = emissionM * transitionXM * bM[i + 1][j + 1];
            double backwardXX = emissionX * transitionXX * bX[i][j + 1];

            bX[i][j] = backwardMX + backwardXX;

//        }


    }

    public void sumbY(int i, int j, double[][] bM, double[][] bY) {

//        System.out.println(" BY: Cell is " + i + " " + j + " seq1char is " + seq1.charAt(i - 1) + " seq2char is " + seq2.charAt(j-1));


//        if (i + 1 > seq1.length()) {
//            bY[i][j] = 0;
//        } else {
            double emissionM;
            if (i > seq1.length() - 1 || j > seq2.length() - 1){
                 emissionM = 0;
            }
            else {
                 emissionM = ExampleModel.getDistance(seq1.charAt(i), seq2.charAt(j));
            }
            double emissionY = 0.25;

            double backwardMY = emissionM * transitionYM * bM[i + 1][j + 1];
            double backwardYY = emissionY * transitionYY * bY[i + 1][j];

            bY[i][j] = backwardMY + backwardYY;

//        }
    }

    public void performMEA(){
        double[][] fM = this.forwardAlgorithm();
        double[][] bM = this.backwardAlgorithm();
        double[][] pM = calcPosteriorMatrix(fM, bM);

        POGraphAlignment poGraphAlignment = new POGraphAlignment(seq1, seq2, 0, 0, pM, true);



    }

    public void printAllMatrices(){

        printMatrix(vM);
        printMatrix(vX);
        printMatrix(vY);
    }

    public static void printMatrix(double[][] matrix) {
        NumberFormat formatter = new DecimalFormat();
        formatter = new DecimalFormat("0.#####E0");

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", formatter.format(matrix[i][j])) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }

    public static void printMatrix(int[][] matrix) {
        NumberFormat formatter = new DecimalFormat();
        formatter = new DecimalFormat("0.#####E0");

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", formatter.format(matrix[i][j])) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }

    public static void printMatrix(String[][] matrix) {


        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", matrix[i][j]) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }

    public HashMap<Integer, Integer> getAlignedPairs(){
        return this.alignedPairs;

    }

    public static double[][] calcPosteriorMatrix(double[][] fM, double[][]bM){
        System.out.println("Forward Matrix:");
        printMatrix(fM);
        System.out.println("Backward Matrix:");

        printMatrix(bM);

        double[][] pM = new double[fM.length][fM[0].length];

        for (int i = 0; i < fM.length; i++) {
            for (int j = 0; j < fM[0].length; j++) {
                pM[i][j] = fM[i][j] * bM[i][j];
            }


            }
        System.out.println("Posterior Matrix:");

        printMatrix(pM);

        return pM;


        }

}


