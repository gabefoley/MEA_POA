package Alignment;

import SubstitutionModels.SubstitutionMatrix;
import jdk.internal.util.xml.impl.Pair;
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

import java.io.*;
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

        // Setup pairwise inputs
//        String pairwiseQuery = "TTACGACC";
//        String pairwiseTarget = "CCTTAACC";

//        String pairwiseQuery = "PPNGRNRR";
//        String pairwiseTarget = "RRPPNNRR";

        double tau = 0.01;
        double epsilon = 0.1;
        double delta = 0.2;
        double eta = 0.1111111111;
        double openPenalty = calculateOpenPenalty(tau, epsilon, delta, eta);
        double extensionPenalty = calculateExtensionPenalty(epsilon, eta);
        double calculatedEpsilon = calculateEpsilon(extensionPenalty, eta);

//        System.out.println("Open penalty is " + openPenalty);
//        System.out.println("Extension penalty is " + extensionPenalty);
//        System.out.println("Calculated epsilon is " + calculatedEpsilon);


        //Profile too large problem
        String pairwiseQuery = "MYSFPNSFRFGWSQAGFQSEMGTPGSEDPNTDWYKWVHDPENMAAGLVSGDLPENGPGYWGNYKTFHDNAQKMGLKIARLNSEW" +
                "SRQFPNPLPRPQNFDESKQDVTEVEINENELKRLDEYANKDALNHYREIFKDLKSRGLYFIQNMYHWPLPLWLHDPIRVRRGDFTGPSGWLSTRTVYEF" +
                "ARFSAYTAWKFDDLVDEYSTMNEPNVVGGLGYVGVKSGFPPGYLSFELSRRAMYNIIQAHARAYDGIKSVSKKPVGIIYANSSFQPLTDKDMEAVEMAE" +
                "NDNRWWFFDAIIRGEITRGNEKIVRDDLKGRLDWIGVNYYTRTVVKRTEKGYVSLGGYGHGCERNSVSLAGLPTSDFGWEFFPEGLYDVLTKYWNRYHL" +
                "YMYVTENGIADDADYQRPYYLVSHVYQVHRAINSGADVRGYLHWSLADNYEWASGFSMRFGLLKVDYNTKRLYWRPSALVYREIATNGAITDEIEHLNS" +
                "VPPVKPLRH";
        String pairwiseTarget = "IPRWRGFNLLEAFSIKSTGNFKEEDFLWXAQWDFNFVRIPXCHLLWSDRGNPFIIREDFFEKIDRVIFWGEKYGIHICISLHR" +
                "APGYSVNKEVEEKTNLWKDETAQEAFIHHWSFIARRYKGISSTHLSFNLINEPPFPDPQIXSVEDHNSLIKRTITEIRKIDPERLIIIDGLGYGNIPVD" +
                "DLTIENTVQSCRGYIPFSVTHYKAEWVDSKDFPVPEWPNGWHFGEYWNREKLLEHYLTWIKLRQKGIEVFCGEXGAYNKTPHDVVLKWLEDLLEIFKTL" +
                "NIGFALWNFRGPFGILDSERKDVEYEEWYGHKLDRKXLELLRKY";

//        String pairwiseQuery = "PPGRNAAKKKKKKKKPRETTANSGAGAGKKKKRDDWWWLHN";
//        String pairwiseTarget = "WWWWWWPPGRNAWGAGAGAAAAWWWRTEDNLKNFFLNWWWW";

//        String pairwiseQuery = "PPGRR";
//        String pairwiseTarget = "NNPP";

        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);

//        double gapOpen = openPenalty;
//        double gapExtend = extensionPenalty;

        int gapOpen = 4;
        int gapExtend = 1;

        //Setup HMM inputs
//        double tau = 0.1;
//        double epsilon = 0.1;
//        double delta = 0.2;



        double emissionX = 0.05;
        double emissionY = 0.05;

        //Setup MSA inputs
        int multiGapOpen = -4;
        int multiGapExtend = -1;


//        String multiSeq1 = "GMKKWPRDCC";
//        String multiSeq2 = "WPGMSVTNDC";
//        String multiSeq3 = "NAREEREDCA";
//        String multiSeq4 = "DDPGAEERED";

        String multiSeq1 = "GGPAWE";
        String multiSeq2 = "AWENNP";
        String multiSeq3 = "GRPWES";

//        String multiSeq4 = "DDPGAEERED";


//        String multiSeq1 = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGEKKKKKKKKK";
//        String multiSeq2 = "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK";
//        String multiSeq3 = "MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK";
//        String multiSeq4 = "MHIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK";



        SubstitutionMatrix blosum62 = new SubstitutionMatrix("blosum62");
        SubstitutionMatrix blosum62Probs = new SubstitutionMatrix("blosum62Probs");
        SubstitutionMatrix exampleModel = new SubstitutionMatrix("exampleModel");
        SubstitutionMatrix blosum62estimatedlambda = new SubstitutionMatrix("blosum62estimatedlambda");
        SubstitutionMatrix blosum62LatestProbs = new SubstitutionMatrix("blosum62LatestProbs");
        SubstitutionMatrix blosum62EstimatedWithX = new SubstitutionMatrix("blosum62EstimatedWithX");

        // Run pairwise alignments
//        runBioJavaPairwise(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend);
//        runPairwise(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend, blosum62);
//        runViterbiPairwiseOnProfile(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, emissionX, emissionY, blosum62EstimatedWithX);
//        runMEAPairwiseOnProfile(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, emissionX, emissionY, blosum62EstimatedWithX);

        // Run Baum Welch pairwise
//        runBaumWelch();
//        runBWDefault(pairwiseQuery, pairwiseTarget, blosum62estimatedlambda);
        runBWMulti(multiSeq1, multiSeq2, multiSeq3, blosum62EstimatedWithX);


        // Run multiple sequence alignments
//        runBioJavaMSA(multiSeq1, multiSeq2, multiSeq3, multiGapOpen, multiGapExtend);
//        runMSA(multiSeq1, multiSeq2, multiSeq3, multiGapOpen, multiGapExtend, blosum62);
//        runViterbiMSA(multiSeq1, multiSeq2, multiSeq3, tau, epsilon, delta, emissionX, emissionY, blosum62Probs);
//        runMEAMSA(multiSeq1, multiSeq2, multiSeq3, tau, epsilon, delta, emissionX, emissionY, blosum62Probs);








    }

//    private static void topsParameters(){
//        double[] start = {0.6814756989, 0.00008615339902, 0.1591759622 };
//
//        double[][] transition = { }
//    }
//}
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

    private static void runBWDefault(String pairwiseQuery, String pairwiseTarget, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("\nViterbi alignment with Baum Welch parameter estimation: ");
        String[] seqArray = new String[2];

//        String[] balibaseSeqs = getBalibaseSeqs();

        seqArray[0] = pairwiseQuery;
        seqArray[1] = pairwiseTarget;



//        BaumWelch bw = new BaumWelch(seqArray, "protein");
////
//        double[] start = bw.getStart();
//        double[][] transition = bw.getTransition();
//        double[][] emission = bw.getEmission();
//
//
//
//
//
//        HashProfile first = new HashProfile(multiSeq1);
//        HashProfile second = new HashProfile(multiSeq2);

        PairHMM pairHMM = new PairHMM(seqArray, subMatrix, "protein");
        HashProfile alignment = pairHMM.getViterbiAlignment();

        System.out.println(alignment);

    }

    public static void runBWMulti(String multiSeq1, String multiSeq2, String multiSeq3, SubstitutionMatrix subMatrix){
        String[] seqArray = new String[3];

        seqArray[0] = multiSeq1;
        seqArray[1] = multiSeq2;
        seqArray[2] = multiSeq3;

        PairHMM pairHMM = new PairHMM(seqArray, subMatrix, "protein");


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

    private static void runPairwise(String pairwiseQuery, String pairwiseTarget, int gapOpen, int gapExtend, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("NW Alignment: ");

        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);


        Alignment alignment = new Alignment(pairwiseQuery, pairwiseTarget,
                -1 * gapOpen, -1 * gapExtend, subMatrix, false);

        System.out.println(alignment.getUpdatedProfile().printSeqs());

//        System.out.println(alignment.getUpdatedProfile());


    }



    private static void runViterbiPairwiseOnProfile(String pairwiseQuery, String pairwiseTarget, double tau, double epsilon,
                                           double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("\nViterbi Pairwise using PairHMM: ");

        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);



        PairHMM pairHMM = new PairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, emissionX, emissionY, subMatrix);
        HashProfile alignment = pairHMM.getViterbiAlignment();
        System.out.println(alignment.printSeqs());

    }

    private static void runMEAPairwiseOnProfile(String pairwiseQuery, String pairwiseTarget, double tau, double epsilon,
                                       double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix ){
        System.out.println("--------------------------");
        System.out.println("\nMaximum Expected Accuracy Alignment using PairHMM: ");

        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);


        PairHMM pairHMM = new PairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, emissionX, emissionY, subMatrix);
        HashProfile alignment = pairHMM.getMEAAlignment();

        System.out.println(alignment.printSeqs());

//        Alignment meaAlignment = pairHMM.getMEAAlignment();



//        PairHMM pairHMM = new PairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);

    }

    private static void runBioJavaMSA(String multiSeq1, String multiSeq2, String multiSeq3,
                                      int multiGapOpen, int multiGapExtend) throws CompoundNotFoundException {
        System.out.println("--------------------------");
        System.out.println("\nBioJava Multiple Sequences Alignment: ");

        ProteinSequence firstProt = new ProteinSequence(multiSeq1, AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence secondProt = new ProteinSequence(multiSeq2, AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence thirdProt = new ProteinSequence(multiSeq3, AminoAcidCompoundSet.getAminoAcidCompoundSet());
//        ProteinSequence fourthProt = new ProteinSequence(multiSeq4, AminoAcidCompoundSet.getAminoAcidCompoundSet());

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
//        list.add(fourthProt);
        Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(list, penalty,
                SubstitutionMatrixHelper.getBlosum62());

        ConcurrencyTools.shutdown();

        System.out.println(profile);


    }

    private static void runMSA(String multiSeq1, String multiSeq2, String multiSeq3,
                               int multiGapOpen, int multiGapExtend, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("\nMultiple Sequences Alignment: ");
        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
//        HashProfile fourth = new HashProfile(multiSeq4);




        System.out.println("\nFirst alignment:");
        Alignment firstAlignment = new Alignment(first, second, multiGapOpen, multiGapExtend, subMatrix, false);
        HashProfile firstProfile = firstAlignment.getUpdatedProfile();

        System.out.println(firstProfile.printSeqs());

        System.out.println("\nSecond alignment:");


        Alignment secondAlignment = new Alignment(firstProfile, third, multiGapOpen, multiGapExtend,  subMatrix, false);
        HashProfile secondProfile = secondAlignment.getUpdatedProfile();
        System.out.println(secondProfile.printSeqs());


//        System.out.println("\nThird alignment:");
//
//        Alignment thirdAlignment = new Alignment(secondProfile, fourth, multiGapOpen, multiGapExtend,  subMatrix, false);
//        HashProfile thirdProfile = thirdAlignment.getUpdatedProfile();
//
//        System.out.println(thirdProfile.printSeqs());

    }


    private static void runViterbiMSA(String multiSeq1, String multiSeq2, String multiSeq3, double tau, double epsilon,
                                      double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix){

        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
//        HashProfile fourth = new HashProfile(multiSeq4);


        System.out.println("--------------------------");
        System.out.println("\nViterbi Multiple Sequences Alignment: ");



        System.out.println("\nFirst alignment:");
        PairHMM firstAlignment = new PairHMM(first, second, tau, epsilon, delta, emissionX, emissionY, subMatrix);
        HashProfile firstProfile = firstAlignment.getViterbiAlignment();
        System.out.println(firstProfile.printSeqs());

        System.out.println("\nSecond alignment:");


        PairHMM secondAlignment = new PairHMM(firstProfile, third, tau, epsilon, delta, emissionX, emissionY, subMatrix);
        HashProfile secondProfile = secondAlignment.getViterbiAlignment();
        System.out.println(secondProfile.printSeqs());


//        System.out.println("\nThird alignment:");
//
//        PairHMM thirdAlignment = new PairHMM(secondProfile, fourth, tau, epsilon, delta, emissionX, emissionY, subMatrix);
//        HashProfile thirdProfile = thirdAlignment.getViterbiAlignment();
//        System.out.println(thirdProfile.printSeqs());

    }


    private static void runMEAMSA(String multiSeq1, String multiSeq2, String multiSeq3, double tau, double epsilon,
                                  double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("\nMaximum Expected Accuracy Multiple Sequences Alignment: ");

        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
//        HashProfile fourth = new HashProfile(multiSeq4);

        System.out.println("\nFirst alignment:");
        PairHMM firstAlignment = new PairHMM(first, second, tau, epsilon, delta, emissionX, emissionY, subMatrix);
        HashProfile firstProfile = firstAlignment.getMEAAlignment();
        System.out.println(firstProfile.printSeqs());

        System.out.println("\nSecond alignment:");


        PairHMM secondAlignment = new PairHMM(firstProfile, third, tau, epsilon, delta, emissionX, emissionY, subMatrix);
        HashProfile secondProfile = secondAlignment.getMEAAlignment();
        System.out.println(secondProfile.printSeqs());

//
//        System.out.println("\nThird alignment:");
//
//        PairHMM thirdAlignment = new PairHMM(secondProfile, fourth, tau, epsilon, delta, emissionX, emissionY, subMatrix);
//        HashProfile thirdProfile = thirdAlignment.getMEAAlignment();
//        System.out.println(thirdProfile.printSeqs());



    }

    public static double calculateEta(String[] seqArray){
        double length = 0;
        for (String seq : seqArray){
            length += seq.length();
        }
        double meanLength = length / seqArray.length;

        double eta = 1 / (meanLength + 1);
        return eta;
    }

    public static double calculateOpenPenalty(double tau, double epsilon, double delta, double eta){



        double openPenalty = - Math.log((delta * (1 - epsilon - tau))/((1-eta) * (1-2*delta - tau)));

        return openPenalty;


    }

    public static double calculateExtensionPenalty(double epsilon, double eta ){

        double oneMinusEta = 1 - eta;
        double epsilondivided = epsilon / oneMinusEta;
        double negativelog = - Math.log(epsilondivided);

//        System.out.println("One minus eta:" + oneMinusEta);
//        System.out.println("epsilonDivided" + epsilondivided);
//        System.out.println("negativelog" + negativelog);

        double extensionPenalty = - Math.log((epsilon) / (1 - eta));

        return extensionPenalty;
    }

    public static double calculateEpsilon(double extensionPenalty, double eta){
        double epsilon = (1 - eta) / (Math.pow(Math.E, extensionPenalty));

        return epsilon;
    }

    public static String[] getBalibaseSeqs() throws  IOException{
        String dir = "/Users/gabe/Dropbox/PhD/Programs/bench1.0/bali3/in";
        File file = new File(dir);

        String[] inputs = file.list();

        for (String input: inputs){
            System.out.println("Reading from: " + input);

            FileReader fr = new FileReader(dir + "/" +  input);

            BufferedReader br = new BufferedReader(fr);
            String[] seqs;
            if (! br.readLine().startsWith(">")){
//                seqs(br.readLine();

            }
        }

        return inputs;

    }

}
