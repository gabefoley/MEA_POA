//package Alignment;
//
//import SubstitutionModels.ExampleModel;
//
///**
// * Created by gabe on 9/08/2016.
// */
//public class PairNW {
//
//    private double transitionMM, transitionMX, transitionMY, transitionXX, transitionYY, transitionXM, transitionYM;
//    private String seq1, seq2;
//    private double[][] vM;
//    private double[][] vX;
//    private double[][] vY;
//    private String[][] tracebackM, tracebackX, tracebackY;
//
//    private int longerSide;
//    private int shorterSide;
//    private String longerSeq;
//    private String shorterSeq;
//
//    private double tau;
//
//
//    public PairNW(String seq1, String seq2, double tau, double epsilon, double delta) {
//
//        this.seq1 = seq1;
//        this.seq2 = seq2;
//
//
//
//        // Probability of aligning A,G,C,T with a gap
//
//
//        this.vM = new double[seq1.length() + 1][seq2.length() + 1];
//        this.vX = new double[seq1.length() + 1][seq2.length() + 1];
//        this.vY = new double[seq1.length() + 1][seq2.length() + 1];
//
//        this.tracebackM = new String[seq1.length() + 1][seq2.length() + 1];
//        this.tracebackX = new String[seq1.length() + 1][seq2.length() + 1];
//        this.tracebackY = new String[seq1.length() + 1][seq2.length() + 1];
//
//
//        vM[0][0] = 1;
//
//
//
//        for (int i = 1; i <= seq1.length() ; i++) {
//            for (int j = 1; j <= seq2.length(); j ++) {
//                fillVM(i, j, vM, vX, vY, tracebackM);
//                fillVX(i, j, vM, vX, tracebackX);
//                fillVY(i, j, vM, vY, tracebackY);
//            }
//
//        }
//
//
//        traceback();
//
//    }
//
//
//    // Fill out the match matrix
//
//    public void fillVM(int i, int j, double[][] vM, double[][] vX, double[][] vY, String[][] tracebackM) {
//
//
//
//        double emissionM = ExampleModel.getDistance(seq1.charAt(i-1), seq2.charAt(j - 1));
//
//        // Get the actual costs for transitioning for each of the states
//        double currentTransitionMM = transitionMM * vM[i - 1][j - 1];
//        double currentTransitionXM = transitionXM * vX[i - 1][j - 1];
//        double currentTransitionYM = transitionYM * vY[i - 1][j - 1];
//
////            System.out.println(i + " " + k);
////            System.out.println(currentTransitionMM * emissionM + " " + currentTransitionXM * emissionM + " " + currentTransitionYM * emissionM);
//
//
//
//        // Work out the optimal cost and set the cell
//        if (currentTransitionMM > currentTransitionXM && currentTransitionMM > currentTransitionYM) {
//            vM[i][j] = currentTransitionMM * emissionM;
//            tracebackM[i][j] = "M";
//        } else if (currentTransitionXM > currentTransitionMM && currentTransitionXM > currentTransitionYM) {
//            vM[i][j] = currentTransitionXM * emissionM;
//            tracebackM[i][j] = "X";
//        } else {
//            vM[i][j] = currentTransitionYM * emissionM;
//            tracebackM[i][j] = "Y";
//
//        }
//
//
//    }
//
//
////    }
//
//
//    // Fill out the gap in X matrix
//    public void fillVX(int i, int j, double[][] vM, double[][] vX, String[][] tracebackX) {
//        double emissionX = 0.25;
//
//
////        for (int k = 1; k < j; k++) {
//
//
//        double currentTransitionXX = transitionXX * vX[i][j - 1];
//        double currentTransitionMX = transitionMX * vM[i][j - 1];
//
//
//        if (currentTransitionXX > currentTransitionMX) {
//            vX[i][j] = currentTransitionXX * emissionX;
//            tracebackX[i][j] = "X";
//
//        } else {
//            vX[i][j] = currentTransitionMX * emissionX;
//            tracebackX[i][j] = "M";
//
//        }
//
////        }
//
//    }
//
//
//    // Fill out the gap in Y matrix
//    public  void fillVY(int i, int j, double[][] vM, double[][] vY, String[][] tracebackY) {
//
//        double emissionY = 0.25;
////        for (int k = 1; k < j; k++) {
//
//        double currentTransitionYY = transitionYY * vY[i - 1][j];
//        double currentTransitionMY = transitionMY * vM[i - 1][j];
//
//        if (currentTransitionYY > currentTransitionMY) {
//            vY[i][j] = currentTransitionYY * emissionY;
//            tracebackY[i][j] = "Y";
//
//        } else {
//            vY[i][j] = currentTransitionMY * emissionY;
//            tracebackY[i][j] = "M";
//
//        }
//
////        }
//
//    }
//
//    public void traceback(){
//
////        printAllMatrices();
//
//        String seq1Output = "";
//        String seq2Output = "";
//        String lastState;
//        int i = seq1.length();
//        int j = seq2.length();
//        int n = longerSide - 1;
//
//
//        while ((i > 0) && (j > 0)) {
//            if ((vM[i][j] > vX[i][j]) && (vM[i][j] > vY[i][j])) {
//                seq1Output = seq1.substring(i-1, i) + seq1Output;
//                seq2Output = seq2.substring(j-1, j) + seq2Output;
//                lastState = tracebackM[i][j];
//                i--;
//                j--;
//
//            }
//
//            else if((vX[i][j]) > vM[i][j] && (vX[i][j]) > vY[i][j]){
//                seq1Output = "-" + seq1Output;
//                seq2Output = seq2.substring(j-1, j) + seq2Output;
//
//                lastState = tracebackX[i][j];
//                j--;
//            }
//
//            else {
//                seq1Output = seq1.substring(i-1, i) + seq1Output;
//                seq2Output = "-" + seq2Output;
//                lastState = tracebackY[i][j];
//                i--;
//            }
//
//            n--;
//
//            while ((i > 0) && (j > 0)) {
//                if (lastState == "M"){
//                    seq1Output = seq1.charAt(i-1) + seq1Output;
//                    seq2Output = seq2.charAt(j-1) + seq2Output;
//                    lastState = tracebackM[i][j];
//                    i--;
//                    j--;
//                }
//                else if (lastState == "Y"){
//                    seq1Output = seq1.charAt(i-1) + seq1Output;
//                    seq2Output = "-" + seq2Output;
//                    lastState = tracebackY[i][j];
//                    i--;
//                }
//
//                else {
//                    seq1Output = "-" + seq1Output;
//                    seq2Output = seq2.charAt(j-1) + seq2Output;
//                    lastState = tracebackX[i][j];
//                    j--;
//                }
//                n--;
//            }
//            n++;
////            System.out.println(n);
//
//
//            System.out.println(seq1Output);
//            System.out.println(seq2Output);
//        }
//
//    }
//}
