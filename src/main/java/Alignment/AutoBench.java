package Alignment;

/**
 * Automatic alignment, benchmarking and result visualisation for sequence alignment benchmarks
 */

import java.io.*;

import SubstitutionModels.SubstitutionMatrix;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;

import org.biojava.nbio.core.util.ConcurrencyTools;
import org.rosuda.JRI.Rengine;


import java.io.File;
import java.util.ArrayList;

public class AutoBench {
    public static void main(String args[]) throws Exception {

        String dir = args[0];
        String inputDir = dir + "/in";
        String refDir = dir + "/ref";
        String[] inputs;

        ArrayList<String> seqsList = new ArrayList<String>();
        ArrayList<Double> scoresList = new ArrayList<Double>();
        String seqString = "";
        String scoreString = "";

        try {
            File file = new File(refDir);

            // array of files and directory
            inputs = file.list();

            // For each file in the reference directory
            for (String input : inputs) {

                System.out.println("Reading from file : " + input);

                ArrayList<String[]> refAlignment = getSeqs(refDir, input);
                ArrayList<String[]> seqs = getSeqs(inputDir, input);

                // Get the character information from the conserved columns in the reference alignment
                ArrayList<String[]> goodCols = getGoodCols(refAlignment);

                // Align the input sequences
                int multiGapOpen = -10;
                int multiGapExtend = -1;
                SubstitutionMatrix subMatrix = new SubstitutionMatrix("blosum62");

                // Get BioJavaMSA
//                Profile<ProteinSequence, AminoAcidCompound> bioJavaMSA = getBioJavaMSA(seqs, multiGapOpen, multiGapExtend);

                // Get MSA
                HashProfile MSA = getMSA(seqs, multiGapOpen, multiGapExtend, subMatrix);


                // Get the good columns from calculated alignment
                int msaCount = 0;
//                int bioJavaMSACount = 0;

                for (String[] goodCol: goodCols) {

                    if (goodCol[1].equals(MSA.getColumn(Integer.valueOf(goodCol[0])))) {
                        msaCount += 1;
                    }


//                    String bioJavaMSAColOrdered = "";
//
//                    if (bioJavaMSA.getCompoundsAt(Integer.valueOf(goodCol[0]) + 1) != null){
//                        List<AminoAcidCompound> bioJavaMSACol = bioJavaMSA.getCompoundsAt(Integer.valueOf(goodCol[0]) + 1);
//
////                        bioJavaMSAColOrdered = "";
//                        for (AminoAcidCompound aa : bioJavaMSACol){
//                            bioJavaMSAColOrdered += aa;
//                        }
//
//                        char[] charArray = bioJavaMSAColOrdered.toCharArray();
//
//                        char temp = charArray[1];
//                        charArray[1] = charArray[2];
//                        charArray[2] = temp;
//
//                        bioJavaMSAColOrdered = new String (charArray);
//                    }
//
//                    if (goodCol[1].equals(bioJavaMSAColOrdered)) {
//                        bioJavaMSACount += 1;
//                    }


//                    System.out.println(goodCol[0]);
//                    System.out.println(goodCol[1]);
//
//                    System.out.println(bioJavaMSA.getCompoundsAt(Integer.valueOf(goodCol[0]) + 1));
//                    System.out.println(bioJavaMSAColOrdered);
//                    System.out.println(MSA.getColumn(Integer.valueOf(goodCol[0])));

                }

                double msaCS = (double) msaCount / goodCols.size() * 100.0;
//                double bioJavaMSACS = (double) bioJavaMSACount / goodCols.size() * 100.0;


                System.out.printf("MSA Column Score for file " + input + " is: %.2f%%%n ", msaCS);
//                System.out.printf("BioJava MSA Column Score for file " + input + " is: %.2f%%%n ",  bioJavaMSACS);
                System.out.println();

                seqsList.add(input);
                scoresList.add(msaCS);





                seqString = "c(";

                for (String seq: seqsList){
                    seqString += "'" + seq + "',";
                }

                seqString = seqString.substring(0, seqString.lastIndexOf(","));

                seqString += ")";

                scoreString = "c(";

                for (Double score: scoresList){
                    scoreString += score + ",";
                }

                scoreString = scoreString.substring(0, scoreString.lastIndexOf(","));

                scoreString += ")";

                System.out.println(seqString);
                System.out.println(scoreString);














                // Compare calculated alignment to reference alignment

            }

            Rengine engine = new Rengine(new String[] { "--no-save"}, false, null);

            engine.eval("require(ggplot2)");
            engine.eval("Seqs <-" + seqString);
            engine.eval("Scores <-" + scoreString);
            engine.eval("joined <- data.frame(Seqs, Scores)");
            engine.eval("plot <- ggplot(data=joined, aes(x=Seqs, y=Scores, group=1)) + geom_line()");
            engine.eval("ggsave('/Users/gabe/Dropbox/testScores.png', plot)");
            System.out.println("Done");

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Get the conserved column numbers and their content from the benchmark database
     * @param seqList List of sequences to get column numbers and content for
     * @return ArrayList of String arrays representing the column number and their content
     * @throws IOException
     */
    public static ArrayList<String[]> getGoodCols(ArrayList<String[]> seqList) throws IOException{

        String seq = seqList.get(0)[1];
        ArrayList<String[]> goodCols = new ArrayList<String[]>();


        for (int i = 0; i < seq.length(); i ++) {
            String colChar = "";
            if (Character.isUpperCase(seq.charAt(i))) {
                for (String[] aSeqList : seqList) {
                    colChar += aSeqList[1].charAt(i);
                }
                String[] colPos = new String[2];
                colPos[0] = (String.valueOf(i));
                colPos[1] = colChar;

                goodCols.add(colPos);



            }
        }

        return goodCols;

    }

    /**
     * Get the sequences from a file
     * @param dir Directory where the sequences are
     * @param input Filename of the file within the directory where the sequences are
     * @return ArrayList of String array representing the sequences from the file
     * @throws IOException
     */
    public static ArrayList<String[]> getSeqs(String dir, String input) throws IOException{

        FileReader fr = new FileReader(dir + "/" +  input);

        BufferedReader br = new BufferedReader(fr);
        String[] seqs;
        ArrayList<String[]> seqList = new ArrayList<String[]>();
        int count = 0;
        while (( seqs = getNextSequence(br)) != null){
            seqList.add(count, seqs);
            count++;
        }

        return seqList;
    }

    /**
     * Returns the next sequence in an opened file
     * @param br The BufferedReader object representing the open file stream
     * @return A string array with the sequence name and the full sequence
     * @throws IOException
     */
    public static String[] getNextSequence(BufferedReader br) throws IOException {
        String[] fastaSeq = new String[2];

        int BUFFER_SIZE = 8192;


        String line = br.readLine();

        if (line == null){
            return null;
        } else if (line.startsWith(">")){
            fastaSeq[0] = line.replaceFirst(">", "");

        } else {
            System.err.println(line);
            System.err.println("Invalid sequence name or format");
            System.exit(1);
        }

        StringBuilder sb = new StringBuilder();
        while(true){
            br.mark(BUFFER_SIZE);

            line = br.readLine();
            if(line == null || line.startsWith(">")){
                break;
            }
            else {
                sb.append(line);
            }
        }

        String sequence = sb.toString();
        fastaSeq[1] = sequence;
        br.reset();
        return fastaSeq;
    }

    /**
     * Calculate a multiple sequence alignment using BioJava
     * @param seqs Sequences to align
     * @param multiGapOpen Gap opening penalty
     * @param multiGapExtend Gap extension penalty
     * @return Profile representing the multiple sequence alignment
     * @throws CompoundNotFoundException
     */
    public static Profile<ProteinSequence, AminoAcidCompound> getBioJavaMSA(ArrayList<String[]> seqs, int multiGapOpen, int multiGapExtend) throws CompoundNotFoundException {
        ArrayList<ProteinSequence> seqList = new ArrayList<ProteinSequence>();

        GapPenalty penalty = new SimpleGapPenalty(multiGapOpen,multiGapExtend);

        for (String[] seq: seqs){
            ProteinSequence proteinSequence = new ProteinSequence(seq[1], AminoAcidCompoundSet.getAminoAcidCompoundSet());
            seqList.add(proteinSequence);

        }
        Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(seqList, penalty,
                SubstitutionMatrixHelper.getBlosum62());

        ConcurrencyTools.shutdown();

        return profile;
    }

    /**
     * Calculate a multiple sequence alignment using Needleman Wunsch dynamic programming
     * @param seqs Sequences to align
     * @param multiGapOpen Gap opening penalty
     * @param multiGapExtend Gap extension penalty
     * @param subMatrix Substitution matrix
     * @return Profile representing the multiple sequence alignment
     */
    public static HashProfile getMSA(ArrayList<String[]> seqs, int multiGapOpen, int multiGapExtend, SubstitutionMatrix subMatrix){

        ArrayList<HashProfile> seqList = new ArrayList<HashProfile>();

        for (String[] seq: seqs){
            HashProfile profile = new HashProfile(seq[1]);
            seqList.add(profile);
        }

        return alignProfilesRecursive(seqList.get(0), seqList.get(1), seqList, multiGapOpen, multiGapExtend, subMatrix, 1);

    }



//    public static HashProfile alignProfiles(HashProfile first, HashProfile second, int multiGapOpen, int multiGapExtend){
//        Alignment alignment = new Alignment(first, second, multiGapOpen, multiGapExtend, Blosum62Probs.getMatrix(), false);
//        return alignment.getUpdatedProfile();
//
//    }

    /**
     * Method for aligning multiple profiles together
     * @param first First profile to align
     * @param second Second profile to align
     * @param seqList List of all profiles to align
     * @param multiGapOpen Gap opening penalty
     * @param multiGapExtend Gap extension penalty
     * @param subMatrix Substituion matrix
     * @param counter Counter to keep track of which profile in the list to align
     * @return Profile representing the multiple sequence alignment
     */
    public static HashProfile alignProfilesRecursive (HashProfile first, HashProfile second, ArrayList<HashProfile> seqList, int multiGapOpen, int multiGapExtend, SubstitutionMatrix subMatrix, int counter){
        Alignment alignment = new Alignment(first, second, multiGapOpen, multiGapExtend, subMatrix, false);
        HashProfile profile = alignment.getUpdatedProfile();
        counter ++;
        if (counter < seqList.size()){
            return alignProfilesRecursive(profile, seqList.get(counter), seqList, multiGapOpen, multiGapExtend, subMatrix, counter);
        }

        else {

            return profile;


        }

    }

}



//                // Calculate the optimal column scores from the pairs of sequences
//
//                for (int i = 0; i < seqList.size(); i++) {
//                    for (int j = i + 1; j < seqList.size(); j++) {
//                        String seq1 = seqList.get(i)[1];
//                        String seq2 = seqList.get(j)[1];
//
//                        System.out.println(i + " is " + seqList.get(i)[1]);
//                        System.out.println(j + " is " + seqList.get(j)[1]);
//                        for (int k = 0; k < seq1.length(); k++){
//                            if (Character.isUpperCase(seq1.charAt(k))){
//                                System.out.println(k);
//                                System.out.println(seq1.charAt(k));
//                                System.out.println(seq2.charAt(k));
//                            }
//
//
//                        }
//
//
//
//                    }
//                }





