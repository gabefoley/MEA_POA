/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package Alignment;

import SubstitutionModels.Blosum62;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;

public class AlignNW {

    public static void main(String[] args) throws Exception {

//        String query = "HEAGAWGHEE";
////
//        String target = "PAWHEAE";

        String query = "TAG";
//
        String target = "TA";

//        String query = "TTAAGTCC";

//        String target = "TATAAGCATTTA";



        GapPenalty penalty = new SimpleGapPenalty(2, 1);


        PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments.getPairwiseAligner(
                new ProteinSequence(query, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
                new ProteinSequence(target, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
                PairwiseSequenceAlignerType.GLOBAL,
                penalty, SubstitutionMatrixHelper.getBlosum50());

        SequencePair<ProteinSequence, AminoAcidCompound>
                alignment = aligner.getPair();

        System.out.println("BioJava NW Alignment: "+ "\n" + alignment);

        System.out.println("Gabe's NW Alignment: ");

        POGraphAlignment poGraphAlignment = new POGraphAlignment(query, target, 0, 0, Blosum62.getMatrix(), false);




//        int identical = alignment.getNumIdenticals();
//        System.out.println("Number of identical residues: "+ identical);
//        System.out.println("% identical query: "+ identical / (float) query.length());
//        System.out.println("% identical target: "+ identical / (float) target.length());
    }
}