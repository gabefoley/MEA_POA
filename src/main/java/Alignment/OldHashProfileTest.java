package Alignment;

import SubstitutionModels.Blosum62;

/**
 * Created by gabe on 15/08/2016.
 */
public class OldHashProfileTest {
    public static void main(String[] args) {



        HashProfile hashProfile = new HashProfile("AGCAG");

//        System.out.println("Check sequences: ");
//        System.out.println(hashProfile.getSequences().toString());

        HashProfile combinedSeqs = new HashProfile(hashProfile, "AG-TG");

        HashProfile hashProfile2 = new HashProfile("AG-TG");
        HashProfile hashProfile3 = new HashProfile("TGCAG");

        HashProfile hashProfile4 = new HashProfile("-TCGG");
        HashProfile hashProfile5 = new HashProfile("-GCAG");
//        HashProfile hashProfile6 = new HashProfile("-CAG");
        HashProfile hashProfile6 = new HashProfile("AGC--");


        HashProfile combined1 = new HashProfile(hashProfile, hashProfile2);
        HashProfile three1 = new HashProfile(combined1, hashProfile3);

        HashProfile combined2 = new HashProfile(hashProfile4, hashProfile5);
        HashProfile three2 = new HashProfile(combined2, hashProfile6);

        System.out.println("Check sequences: ");
        System.out.println(combined1.toString());


        three1.printProfileArray();
        three2.printProfileArray();

        HashProfile simpleProfile = new HashProfile("AAAGA");

        HashProfile simplerProfile = new HashProfile("G");

        POGraphAlignment alignment = new POGraphAlignment(three1, three2, -2, -1, Blosum62.getMatrix(), false);
//        POGraphAlignment alignment = new POGraphAlignment(simpleProfile, simplerProfile, -2, -1, Blosum62.getMatrix(), false);


//        SimpleProfile simpleProfile = new SimpleProfile("AGCAG");
//
//        SimpleProfile first = new SimpleProfile("AGCAG", "AG-TG");
//
//        SimpleProfile second = new SimpleProfile("GCAG", "-CAG");
//
//        POGraphAlignment alignment = new POGraphAlignment(first, second, -2, -1, Blosum62.getMatrix(), false);

    }

}
