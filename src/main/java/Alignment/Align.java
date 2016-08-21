package Alignment;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by gabe on 31/07/2016.
 */
public class Align {

    public static void main(String[] args){

//        // A,G,C,T
//        double[][] emissions = new double[][]{
//                {0.50, 0.05, 0.15, 0.30},
//                {0.05, 0.50, 0.30, 0.15},
//                {0.15, 0.30, 0.50, 0.05},
//                {0.30, 0.15, 0.05, 0.50}};

//        // Begin, Match, Gap X, Gap Y, End
//
//        double[][] transitions = new double[][]{
//                {0.0, 0.5, 0.2, 0.2, 0.1},
//                {0.0, 0.5, 0.2, 0.2, 0.1},
//                {0.0, 0.8, 0.1, 0.0, 0.1},
//                {0.0, 0.8, 0.0, 0.1, 0.1},
//                {0.0, 0.0, 0.0, 0.0, 0.0}};


//        final String seq1 = "GAT";
//        final String seq2 = "CGAT";

//        final String seq1 = "AGCG";
//        final String seq2 = "AGCAG";

        final String seq1 = "TAG";
        final String seq2 = "TA";

//        final String seq1 = "TTACG";
//        final String seq2 = "TAG";


        double tau = 0.1;
        double epsilon = 0.1;
        double delta = 0.2;

        PairHMM pairHMM = new PairHMM(seq1, seq2, tau, epsilon, delta);
        pairHMM.performMEA();



//        pairHMM.printAllMatrices();
//        double[][] forwardProb = pairHMM.forwardAlgorithm();
//        HashMap<Integer, Integer> alignedPairs = pairHMM.getAlignedPairs();
//
//        System.out.println(alignedPairs);
//
//
//
//
//
//
//        System.out.println("Full probability is " + forwardProb);
//        pairHMM.backwardAlgorithm();






    }
}
