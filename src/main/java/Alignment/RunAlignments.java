package Alignment;

import SubstitutionModels.Blosum62Probs;
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
 * Created by gabe on 23/08/2016.
 */
public class RunAlignments {

    public static void main(String[] args) throws Exception {


        // Setup pairwise inputs
//        String pairwiseQuery = "AACT";
//        String pairwiseTarget = "CT";


//        String pairwiseQuery = "PPGKRNDTG";
//        String pairwiseTarget = "PPGNDTT";

        String pairwiseQuery = "AAGPWRDDFP";
        String pairwiseTarget = "AARDFT";


        int gapOpen = 2;
        int gapExtend = 1;

        //Setup HMM inputs
        double tau = 0.1;
        double epsilon = 0.1;
        double delta = 0.2;

        //Setup MSA inputs
        int multiGapOpen = -10;
        int multiGapExtend = -1;
//
//        String multiSeq1 = "AATG";
//        String multiSeq2 = "TTG";
//        String multiSeq3 = "PTG";
//        String multiSeq4 = "CCG";

        // Elongation issue with MEA / Viterbi
//        String multiSeq1 = "CAACCAGGCATG";
//        String multiSeq2 = "CAATTG";
//        String multiSeq3 = "ACTTCACTG";
//        String multiSeq4 = "CTTCG";


        // Elongation issue with MEA / Viterbi
//        String multiSeq1 = "CCAGCTAG";
//        String multiSeq2 = "AGCC";
//        String multiSeq3 = "CTGGGA";
//        String multiSeq4 = "ATGA";

        String multiSeq1 = "GMKKWPR";
        String multiSeq2 = "WPGMSVTNDC";
        String multiSeq3 = "NAREERED";
        String multiSeq4 = "DDPGA";


//        String multiSeq1 = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGEKKKKKKKKK";
//        String multiSeq2 = "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK";
//        String multiSeq3 = "MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK";
//        String multiSeq4 = "MHIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK";




        // Run pairwise alignments
//        runBioJavaPairwise(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend);
//        runPairwise(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend);
//        runViterbiPairwise(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);
//        runViterbiPairwiseOnProfile(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);
//        runMEAPairwise(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);
//        runMEAPairwiseOnProfile(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);


        // Run multiple sequence alignments
        runBioJavaMSA(multiSeq1, multiSeq2, multiSeq3, multiSeq4, multiGapOpen, multiGapExtend);
        runMSA(multiSeq1, multiSeq2, multiSeq3, multiSeq4, multiGapOpen, multiGapExtend);
        runViterbiMSA(multiSeq1, multiSeq2, multiSeq3, multiSeq4, tau, epsilon, delta);
        runMEAMSA(multiSeq1, multiSeq2, multiSeq3, multiSeq4, tau, epsilon, delta);

        // Run Baum Welch test

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

        double [] start = {.3, .7};

        double [][] transition = {
                {.5, .5},
                {.4, .6},
        };

        double [][] emission = {
                {.2, .8 },
                {.9, .1},
        };
//        BaumWelch bw = new BaumWelch(seqArray, start, transition, emission);
//        System.out.println(bw);


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

    private static void runPairwise(String pairwiseQuery, String pairwiseTarget, int gapOpen, int gapExtend){
        System.out.println("--------------------------");
        System.out.println("NW Alignment: ");


        //TODO: Fix up selecting substitution matrix
        Alignment alignment = new Alignment(pairwiseQuery, pairwiseTarget,
                -1 * gapOpen, -1 * gapExtend, Blosum62Probs.getMatrix(), false);

        System.out.println(alignment.getUpdatedProfile());

//        System.out.println(alignment.getUpdatedProfile());


    }

    private static void runViterbiPairwise(String pairwiseQuery, String pairwiseTarget, double tau, double epsilon,
                                           double delta){
        System.out.println("--------------------------");
        System.out.println("\nViterbi Pairwise: ");

        PairHMM pairHMM = new PairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);
        pairHMM.getViterbiAlignmnet();


    }

    private static void runMEAPairwise(String pairwiseQuery, String pairwiseTarget, double tau, double epsilon,
                                       double delta ){
        System.out.println("--------------------------");
        System.out.println("\nMaximum Expected Accuracy Alignment: ");

        PairHMM pairHMM = new PairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);
        Alignment meaAlignment = pairHMM.getMEAAlignment();

    }

    private static void runViterbiPairwiseOnProfile(String pairwiseQuery, String pairwiseTarget, double tau, double epsilon,
                                           double delta){
        System.out.println("--------------------------");
        System.out.println("\nViterbi Pairwise using ProfilePairHMM: ");

        //TODO: Handle this in the constructor
        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);

        ProfilePairHMM pairHMM = new ProfilePairHMM(profile1, profile2, tau, epsilon, delta);
        HashProfile alignment = pairHMM.getViterbiAlignmnet();
        System.out.println(alignment);

    }

    private static void runMEAPairwiseOnProfile(String pairwiseQuery, String pairwiseTarget, double tau, double epsilon,
                                       double delta ){
        System.out.println("--------------------------");
        System.out.println("\nMaximum Expected Accuracy Alignment using ProfilePairHMM: ");

        //TODO: Handle this in the constructor
        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);

        ProfilePairHMM pairHMM = new ProfilePairHMM(profile1, profile2, tau, epsilon, delta);
        HashProfile alignment = pairHMM.getMEAAlignment();

        System.out.println(alignment);

//        Alignment meaAlignment = pairHMM.getMEAAlignment();

//        System.out.println("Not being called correctly!");


//        ProfilePairHMM pairHMM = new ProfilePairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);

    }

    private static void runBioJavaMSA(String multiSeq1, String multiSeq2, String multiSeq3, String multiSeq4,
                                      int multiGapOpen, int multiGapExtend) throws CompoundNotFoundException {
        System.out.println("--------------------------");
        System.out.println("\nBioJava Multiple Sequence Alignment: ");

        ProteinSequence firstProt = new ProteinSequence(multiSeq1, AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence secondProt = new ProteinSequence(multiSeq2, AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence thirdProt = new ProteinSequence(multiSeq3, AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence fourthProt = new ProteinSequence(multiSeq4, AminoAcidCompoundSet.getAminoAcidCompoundSet());

        //TODO: Get this penalty value calling correctly
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
                               int multiGapOpen, int multiGapExtend){
        System.out.println("--------------------------");
        System.out.println("\nMultiple Sequence Alignment: ");
        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
        HashProfile fourth = new HashProfile(multiSeq4);




        System.out.println("\nFirst alignment:");
        Alignment firstAlignment = new Alignment(first, second, multiGapOpen, multiGapExtend, Blosum62Probs.getMatrix(), false);
        HashProfile firstProfile = firstAlignment.getUpdatedProfile();

        System.out.println(firstProfile);

        System.out.println("\nSecond alignment:");


        Alignment secondAlignment = new Alignment(firstProfile, third, multiGapOpen, multiGapExtend, Blosum62Probs.getMatrix(), false);
        HashProfile secondProfile = secondAlignment.getUpdatedProfile();
        System.out.println(secondProfile);


        System.out.println("\nThird alignment:");

        Alignment thirdAlignment = new Alignment(secondProfile, fourth, multiGapOpen, multiGapExtend, Blosum62Probs.getMatrix(), false);
        HashProfile thirdProfile = secondAlignment.getUpdatedProfile();

        System.out.println(thirdProfile);

    }


    private static void runViterbiMSA(String multiSeq1, String multiSeq2, String multiSeq3, String multiSeq4, double tau, double epsilon,
                                      double delta){
        System.out.println("--------------------------");
        System.out.println("\nViterbi Multiple Sequence Alignment: ");

        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
        HashProfile fourth = new HashProfile(multiSeq4);




        System.out.println("\nFirst alignment:");
        ProfilePairHMM firstAlignment = new ProfilePairHMM(first, second, tau, epsilon, delta);
        HashProfile firstProfile = firstAlignment.getViterbiAlignmnet();
        System.out.println(firstProfile);

        System.out.println("\nSecond alignment:");


        ProfilePairHMM secondAlignment = new ProfilePairHMM(firstProfile, third, tau, epsilon, delta);
        HashProfile secondProfile = secondAlignment.getViterbiAlignmnet();
        System.out.println(secondProfile);


        System.out.println("\nThird alignment:");

        ProfilePairHMM thirdAlignment = new ProfilePairHMM(secondProfile, fourth, tau, epsilon, delta);
        HashProfile thirdProfile = thirdAlignment.getViterbiAlignmnet();
        System.out.println(thirdProfile);

    }


    private static void runMEAMSA(String multiSeq1, String multiSeq2, String multiSeq3, String multiSeq4, double tau, double epsilon,
                                  double delta){
        System.out.println("--------------------------");
        System.out.println("\nMaximum Expected Accuracy Multiple Sequence Alignment: ");

        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
        HashProfile fourth = new HashProfile(multiSeq4);

        System.out.println("\nFirst alignment:");
        ProfilePairHMM firstAlignment = new ProfilePairHMM(first, second, tau, epsilon, delta);
        HashProfile firstProfile = firstAlignment.getMEAAlignment();
        System.out.println(firstProfile);

        System.out.println("\nSecond alignment:");


        ProfilePairHMM secondAlignment = new ProfilePairHMM(firstProfile, third, tau, epsilon, delta);
        HashProfile secondProfile = secondAlignment.getMEAAlignment();
        System.out.println(secondProfile);


        System.out.println("\nThird alignment:");

        ProfilePairHMM thirdAlignment = new ProfilePairHMM(secondProfile, fourth, tau, epsilon, delta);
        HashProfile thirdProfile = thirdAlignment.getMEAAlignment();
        System.out.println(thirdProfile);



    }
}
