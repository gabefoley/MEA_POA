package Alignment;

import java.io.File;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;


/**
 * Created by gabe on 18/07/2016.
 */
public class NeedlemanWunsch {

    public static void main(String[] args) throws CompoundNotFoundException{


        try {
            File subMatrix = new File(args[0]);
            alignPairGlobal(subMatrix);
        } catch (Exception e){
            e.printStackTrace();
        }

    }

    private static void alignPairGlobal(File subMatrix) throws Exception {
        ProteinSequence s1 = new ProteinSequence("AGKPPRRWARS"), s2 = new ProteinSequence("AGNNRWARS");
        SubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>
                (AminoAcidCompoundSet.getAminoAcidCompoundSet(), subMatrix);
        SequencePair<ProteinSequence, AminoAcidCompound> pair = Alignments.getPairwiseAlignment(s1, s2,
                Alignments.PairwiseSequenceAlignerType.GLOBAL, new SimpleGapPenalty(), matrix);

        System.out.println("Result is ");
        System.out.println(pair);

    }
}
