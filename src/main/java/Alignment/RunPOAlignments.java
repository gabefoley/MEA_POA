package Alignment;

import Alignment.Utilities.Sequence;
import SubstitutionModels.SubstitutionMatrix;
import bn.reconstruction.MSA;
import dat.POGraph;

import java.util.ArrayList;
import java.util.List;


/**
 * Created by gabe on 14/11/2016.
 */
public class RunPOAlignments {

    public static void main(String[] args) throws Exception {


        Sequence multiSeq1 = new Sequence("1", "AARS");
        Sequence multiSeq2 = new Sequence("2", "AAMGRS");
        Sequence multiSeq3 = new Sequence("3", "AAMRS");
//
//
//        Sequence[] multiArray = new Sequence[3];
//        multiArray[0] = multiSeq1;
//        multiArray[1] = multiSeq2;
//        multiArray[2] = multiSeq3;

        double tau = 0.1;
        double epsilon = 0.1;
        double delta = 0.2;

        SubstitutionMatrix blosum62EstimatedWithX = new SubstitutionMatrix("blosum62EstimatedWithX");

        double emissionX = 0.25;
        double emissionY = 0.25;

        int open = -4;
        int extend = -2;

        runPONW();
//        runNW(multiSeq1, multiSeq2, multiSeq3, open, extend, blosum62EstimatedWithX);
//        runPOViterbi(multiArray, tau, epsilon, delta, emissionX, emissionY, blosum62EstimatedWithX, "protein");
//        runViterbiMulti(multiArray, tau, epsilon, delta, emissionX, emissionY, blosum62EstimatedWithX, "protein");
//        runPOAlignment();
    }

    public static void runPONW() {

        MSA alignment = new MSA("/Users/gabe/Dropbox/Code/!Files/MEAPOA/singleDist.fasta");
        alignment.saveMSA("/Users/gabe/Dropbox/Code/!Files/MEAPOA/latestoutput");
    }

    private static void runNW(Sequence multiSeq1, Sequence multiSeq2, Sequence multiSeq3,
                               int multiGapOpen, int multiGapExtend, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("\nMultiple Sequences Alignment: ");
        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);


        System.out.println("\nFirst alignment:");
        Alignment firstAlignment = new Alignment(first, second, multiGapOpen, multiGapExtend, subMatrix, false);
        HashProfile firstProfile = firstAlignment.getUpdatedProfile();

        System.out.println(firstProfile.printSeqs());

        System.out.println("\nSecond alignment:");


        Alignment secondAlignment = new Alignment(firstProfile, third, multiGapOpen, multiGapExtend,  subMatrix, false);
        HashProfile secondProfile = secondAlignment.getUpdatedProfile();
        System.out.println(secondProfile.printSeqs());


    }

    public static void runPOViterbi(Sequence[] seqs, double tau, double epsilon,
                                    double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type) {

        System.out.println("-------------------------");
        System.out.println("\n Partial Order Viterbi Alignment: ");
        POGraph initialGraph = new POGraph();

        List<Integer> nodeIds = new ArrayList<>();
        for (int id = 0; id < seqs[0].toString().toCharArray().length; id++) {
            nodeIds.add(id);
        }
        initialGraph.addSequence(0, seqs[0].getID(), seqs[1].getSeq(), nodeIds);

//        POPairHMMUnderflow alignment = new POPairHMMUnderflow(seqArray, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
//        POGraph graph = alignment.getViterbiAlignment();
        System.out.println("Alignment");
//        System.out.println(graph.printSequences();


    }

    private static void runViterbiMulti(Sequence[] seqArray, double tau, double epsilon,
                                        double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type) {
        System.out.println("-------------------------");
        System.out.println("\nViterbi Alignment: ");

        PairHMMUnderflow alignment = new PairHMMUnderflow(seqArray, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile profile = alignment.getViterbiAlignment();

        System.out.println("Alignment");
        System.out.println(profile.printSeqs());


    }


    public static void runPOAlignment() {
        MSA msa = new MSA("/Users/gabe/Dropbox/Code/!Files/MEAPOA/test.fasta");
//        POAlignment alignment = new POAlignment("/Users/gabe/Dropbox/Code/!Files/MEAPOA/test.fasta");
    }
}
