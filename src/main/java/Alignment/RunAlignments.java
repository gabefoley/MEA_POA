package Alignment;

import SubstitutionModels.Blosum62;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;

/**
 * Created by gabe on 23/08/2016.
 */
public class RunAlignments {

    public static void main(String[] args) throws Exception {


        // Setup pairwise inputs
        String pairwiseQuery = "TAGGCC";
        String pairwiseTarget = "TACC";
        int gapOpen = 2;
        int gapExtend = 1;


        // Run pairwise alignments
        runBioJavaPairwise(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend);
        runPairwise(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend);
        runViterbiPairwise();
        runMEAPairwise();


        // Run multiple sequence alignments
        runBioJavaMSA();
        runMSA();
        runViterbiMSA();
        runMEAMSA();





    }

    private static void runBioJavaPairwise(String pairwiseQuery, String pairwiseTarget, int gapOpen, int gapExtend) throws Exception{

        GapPenalty penalty = new SimpleGapPenalty(gapOpen, gapExtend);

        // Setup BioJava Pairwise alignment
        PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments.getPairwiseAligner(
                new ProteinSequence(pairwiseQuery, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
                new ProteinSequence(pairwiseTarget, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
                Alignments.PairwiseSequenceAlignerType.GLOBAL,
                penalty, SubstitutionMatrixHelper.getBlosum50());

        SequencePair<ProteinSequence, AminoAcidCompound>
                alignment = aligner.getPair();

        System.out.println("\nBioJava NW Alignment: " + "\n" + alignment);

    }

    private static void runPairwise(String pairwiseQuery, String pairwiseTarget, int gapOpen, int gapExtend){

        System.out.println("\nNW Alignment: ");


        //TODO: Fix up selecting substitution matrix
        SequenceAligner sequenceAligner = new SequenceAligner(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend, Blosum62.getMatrix(), false);


    }

    private static void runBioJavaMSA(){
        System.out.println("\nBioJava Multiple Sequence Alignment: ");
        System.out.println("Not being called correctly!");


    }

    private static void runMSA(){
        System.out.println("\nMultiple Sequence Alignment: ");
        System.out.println("Not being called correctly!");

    }
    private static void runViterbiPairwise(){
        System.out.println("\nViterbi Pairwise: ");
        System.out.println("Not being called correctly!");

    }

    private static void runViterbiMSA(){
        System.out.println("\nViterbi Multiple Sequence Alignment: ");
        System.out.println("Not being called correctly!");
    }

    private static void runMEAPairwise(){
        System.out.println("\nMaximum Expected Accuracy Alignment: ");
        System.out.println("Not being called correctly!");

    }

    private static void runMEAMSA(){
        System.out.println("\nMaximum Expected Accuracy Multiple Sequence Alignment: ");
        System.out.println("Not being called correctly!");

    }
}
