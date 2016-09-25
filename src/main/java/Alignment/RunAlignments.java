package Alignment;

import SubstitutionModels.SubstitutionMatrix;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.util.ConcurrencyTools;

import java.util.ArrayList;
import java.util.List;

/**
 * Helper method to test running alignments
 */
public class RunAlignments {

    public static void main(String[] args) throws Exception {


        // Setup pairwise inputs
//        String pairwiseQuery = "TTACG";
//        String pairwiseTarget = "TAG";


        String pairwiseQuery = "PPGKRNDTGGTK";
        String pairwiseTarget = "PPGKRTTAAWERAGTKK";

        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);



        int gapOpen = 2;
        int gapExtend = 1;

        //Setup HMM inputs
        double tau = 0.1;
        double epsilon = 0.1;
        double delta = 0.2;

        //Setup MSA inputs
        int multiGapOpen = -10;
        int multiGapExtend = -1;


        String multiSeq1 = "GMKKWPRDCC";
        String multiSeq2 = "WPGMSVTNDC";
        String multiSeq3 = "NAREEREDCA";
        String multiSeq4 = "DDPGAEERED";


//        String multiSeq1 = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGEKKKKKKKKK";
//        String multiSeq2 = "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK";
//        String multiSeq3 = "MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK";
//        String multiSeq4 = "MHIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK";



        SubstitutionMatrix blosum62 = new SubstitutionMatrix("blosum62");
        SubstitutionMatrix blosum62Probs = new SubstitutionMatrix("blosum62Probs");
//        SubstitutionMatrix exampleModel = new SubstitutionMatrix("exampleModel");

        // Run pairwise alignments
        runBioJavaPairwise(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend);
        runPairwise(profile1, profile2, gapOpen, gapExtend, blosum62);
        runViterbiPairwiseOnProfile(profile1, profile2, tau, epsilon, delta, blosum62Probs);
        runMEAPairwiseOnProfile(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, blosum62Probs);


        // Run multiple sequence alignments
        runBioJavaMSA(multiSeq1, multiSeq2, multiSeq3, multiSeq4, multiGapOpen, multiGapExtend);
        runMSA(multiSeq1, multiSeq2, multiSeq3, multiSeq4, multiGapOpen, multiGapExtend, blosum62);
        runViterbiMSA(multiSeq1, multiSeq2, multiSeq3, multiSeq4, tau, epsilon, delta, blosum62Probs);
        runMEAMSA(multiSeq1, multiSeq2, multiSeq3, multiSeq4, tau, epsilon, delta, blosum62Probs);

        // Run Baum Welch
        runBaumWelch();
        runBWDefaults(multiSeq1, multiSeq2, multiSeq3, multiSeq4, blosum62Probs);






    }

    private static void runBaumWelch(){
        double [] start = {.3, .7};

        double [][] transition = {
                {.5, .5},
                {.4, .6},
        };

        double [][] emission = {
                {.2, .8 },
                {.9, .1},
        };

        String [] seqArray = {
                "AA",
                "AA",
                "AA",
                "AA",
                "AT",
                "TT",
                "TA",
                "AA",
                "AA"
        };

        BaumWelch bw = new BaumWelch(seqArray, start, transition, emission, "nucleotide");

        System.out.println(bw);

    }

    private static void runBWDefaults(String multiSeq1, String multiSeq2, String multiSeq3, String multiSeq4, SubstitutionMatrix subMatrix){
        String[] seqArray = new String[4];
        seqArray[0] = multiSeq1;
        seqArray[1] = multiSeq2;
        seqArray[2] = multiSeq3;
        seqArray[3] = multiSeq4;


        BaumWelch bw = new BaumWelch(seqArray, "protein");

        double[] start = bw.getStart();
        double[][] transition = bw.getTransition();
        double[][] emission = bw.getEmission();





        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);

        PairHMM bwPairHMM = new PairHMM(first, second, start, transition, emission, subMatrix);

        System.out.println(bwPairHMM);
    }

    private static void runBioJavaPairwise(String pairwiseQuery, String pairwiseTarget, int gapOpen, int gapExtend)
            throws Exception{
        System.out.println("--------------------------");
        System.out.println("\nBioJava NW Alignment: ");

        GapPenalty penalty = new SimpleGapPenalty(gapOpen, gapExtend);

        // Setup BioJava Pairwise alignment
        PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments.getPairwiseAligner(
                new ProteinSequence(pairwiseQuery, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
                new ProteinSequence(pairwiseTarget, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
                Alignments.PairwiseSequenceAlignerType.GLOBAL,
                penalty, SubstitutionMatrixHelper.getBlosum62());


        SequencePair<ProteinSequence, AminoAcidCompound>
                alignment = aligner.getPair();


        System.out.println(alignment);

    }

    private static void runPairwise(HashProfile pairwiseQuery, HashProfile pairwiseTarget, int gapOpen, int gapExtend, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("NW Alignment: ");


        Alignment alignment = new Alignment(pairwiseQuery, pairwiseTarget,
                -1 * gapOpen, -1 * gapExtend, subMatrix, false);

        System.out.println(alignment.getUpdatedProfile());

//        System.out.println(alignment.getUpdatedProfile());


    }



    private static void runViterbiPairwiseOnProfile(HashProfile pairwiseQuery, HashProfile pairwiseTarget, double tau, double epsilon,
                                           double delta, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("\nViterbi Pairwise using PairHMM: ");



        PairHMM pairHMM = new PairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, subMatrix);
        HashProfile alignment = pairHMM.getViterbiAlignment();
        System.out.println(alignment);

    }

    private static void runMEAPairwiseOnProfile(String pairwiseQuery, String pairwiseTarget, double tau, double epsilon,
                                       double delta, SubstitutionMatrix subMatrix ){
        System.out.println("--------------------------");
        System.out.println("\nMaximum Expected Accuracy Alignment using PairHMM: ");



        PairHMM pairHMM = new PairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, subMatrix);
        HashProfile alignment = pairHMM.getMEAAlignment();

        System.out.println(alignment);

//        Alignment meaAlignment = pairHMM.getMEAAlignment();

//        System.out.println("Not being called correctly!");


//        PairHMM pairHMM = new PairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);

    }

    private static void runBioJavaMSA(String multiSeq1, String multiSeq2, String multiSeq3, String multiSeq4,
                                      int multiGapOpen, int multiGapExtend) throws CompoundNotFoundException {
        System.out.println("--------------------------");
        System.out.println("\nBioJava Multiple Sequences Alignment: ");

        ProteinSequence firstProt = new ProteinSequence(multiSeq1, AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence secondProt = new ProteinSequence(multiSeq2, AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence thirdProt = new ProteinSequence(multiSeq3, AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence fourthProt = new ProteinSequence(multiSeq4, AminoAcidCompoundSet.getAminoAcidCompoundSet());

        GapPenalty penalty = new SimpleGapPenalty(multiGapOpen,multiGapExtend);

//        Alignments.getMultipleSequenceAlignment()
//
//        PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments.getPairwiseAligner(
//                new ProteinSequence(pairwiseQuery, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
//                new ProteinSequence(pairwiseTarget, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
//                Alignments.PairwiseSequenceAlignerType.GLOBAL,
//                penalty, SubstitutionMatrixHelper.getBlosum62());

        List<ProteinSequence> list = new ArrayList<ProteinSequence>();
        list.add(firstProt);
        list.add(secondProt);
        list.add(thirdProt);
        list.add(fourthProt);
        Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(list, penalty,
                SubstitutionMatrixHelper.getBlosum62());

        ConcurrencyTools.shutdown();

        System.out.println(profile);


    }

    private static void runMSA(String multiSeq1, String multiSeq2, String multiSeq3, String multiSeq4,
                               int multiGapOpen, int multiGapExtend, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("\nMultiple Sequences Alignment: ");
        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
        HashProfile fourth = new HashProfile(multiSeq4);




        System.out.println("\nFirst alignment:");
        Alignment firstAlignment = new Alignment(first, second, multiGapOpen, multiGapExtend, subMatrix, false);
        HashProfile firstProfile = firstAlignment.getUpdatedProfile();

        System.out.println(firstProfile);

        System.out.println("\nSecond alignment:");


        Alignment secondAlignment = new Alignment(firstProfile, third, multiGapOpen, multiGapExtend,  subMatrix, false);
        HashProfile secondProfile = secondAlignment.getUpdatedProfile();
        System.out.println(secondProfile);


        System.out.println("\nThird alignment:");

        Alignment thirdAlignment = new Alignment(secondProfile, fourth, multiGapOpen, multiGapExtend,  subMatrix, false);
        HashProfile thirdProfile = thirdAlignment.getUpdatedProfile();

        System.out.println(thirdProfile);

    }


    private static void runViterbiMSA(String multiSeq1, String multiSeq2, String multiSeq3, String multiSeq4, double tau, double epsilon,
                                      double delta, SubstitutionMatrix subMatrix){

        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
        HashProfile fourth = new HashProfile(multiSeq4);


        System.out.println("--------------------------");
        System.out.println("\nViterbi Multiple Sequences Alignment: ");



        System.out.println("\nFirst alignment:");
        PairHMM firstAlignment = new PairHMM(first, second, tau, epsilon, delta, subMatrix);
        HashProfile firstProfile = firstAlignment.getViterbiAlignment();
        System.out.println(firstProfile);

        System.out.println("\nSecond alignment:");


        PairHMM secondAlignment = new PairHMM(firstProfile, third, tau, epsilon, delta, subMatrix);
        HashProfile secondProfile = secondAlignment.getViterbiAlignment();
        System.out.println(secondProfile);


        System.out.println("\nThird alignment:");

        PairHMM thirdAlignment = new PairHMM(secondProfile, fourth, tau, epsilon, delta, subMatrix);
        HashProfile thirdProfile = thirdAlignment.getViterbiAlignment();
        System.out.println(thirdProfile);

    }


    private static void runMEAMSA(String multiSeq1, String multiSeq2, String multiSeq3, String multiSeq4, double tau, double epsilon,
                                  double delta, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("\nMaximum Expected Accuracy Multiple Sequences Alignment: ");

        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
        HashProfile fourth = new HashProfile(multiSeq4);

        System.out.println("\nFirst alignment:");
        PairHMM firstAlignment = new PairHMM(first, second, tau, epsilon, delta, subMatrix);
        HashProfile firstProfile = firstAlignment.getMEAAlignment();
        System.out.println(firstProfile);

        System.out.println("\nSecond alignment:");


        PairHMM secondAlignment = new PairHMM(firstProfile, third, tau, epsilon, delta, subMatrix);
        HashProfile secondProfile = secondAlignment.getMEAAlignment();
        System.out.println(secondProfile);


        System.out.println("\nThird alignment:");

        PairHMM thirdAlignment = new PairHMM(secondProfile, fourth, tau, epsilon, delta, subMatrix);
        HashProfile thirdProfile = thirdAlignment.getMEAAlignment();
        System.out.println(thirdProfile);



    }
}
