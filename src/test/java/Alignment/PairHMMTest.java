package Alignment;

import org.junit.Before;
import org.junit.Test;

/**
 * Created by gabe on 22/08/2016.
 */
public class PairHMMTest {

    @Before


    public void makeSimplePairHMM (){
        double tau = 0.1;
        double epsilon = 0.1;
        double delta = 0.2;
        final String seq1 = "TTG";
        final String seq2 = "CTG";


        PairHMM pairHMM = new PairHMM(seq1, seq2, tau, epsilon, delta);
    }

    @Test
    public void testForwardAlgorithm(){

    }




//        final String seq1 = "TTACG";
//        final String seq2 = "TAG";




}