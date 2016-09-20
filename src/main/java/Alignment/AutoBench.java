package Alignment;

/**
 * Created by gabe on 6/09/2016.
 */

import java.io.*;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import SubstitutionModels.Blosum62Probs;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;

import Alignment.Utilities.Sequence;
import org.biojava.nbio.core.util.ConcurrencyTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class AutoBench {
    public static void main(String args[]) throws Exception {

        String dir = args[0];
        String inputDir = dir + "/in";
        String refDir = dir + "/ref";

        File file = null;
        String[] inputs;



        try {
            // create new file object
            file = new File(refDir);

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

                // Get BioJavaMSA
                Profile<ProteinSequence, AminoAcidCompound> bioJavaMSA = getBioJavaMSA(seqs, multiGapOpen, multiGapExtend);

                // Get MSA
                HashProfile MSA = getMSA(seqs, multiGapOpen, multiGapExtend);

//                System.out.println(bioJavaMSA);
//
//                System.out.println("--------------------------");
//
//
//                for (String[] seq: refAlignment){
//                    System.out.println(seq[1]);
//                }
//
//                System.out.println("--------------------------");
//
//
//                System.out.println(MSA);

                // Get the good columns from calculated alignment
                int msaCount = 0;
                int bioJavaMSACount = 0;

                for (String[] goodCol: goodCols) {

                    if (goodCol[1].equals(MSA.getColumn(Integer.valueOf(goodCol[0])))) {
                        msaCount += 1;
                    }


                    String bioJavaMSAColOrdered = "";

                    if (bioJavaMSA.getCompoundsAt(Integer.valueOf(goodCol[0]) + 1) != null){
                        List<AminoAcidCompound> bioJavaMSACol = bioJavaMSA.getCompoundsAt(Integer.valueOf(goodCol[0]) + 1);

//                        bioJavaMSAColOrdered = "";
                        for (AminoAcidCompound aa : bioJavaMSACol){
                            bioJavaMSAColOrdered += aa;
                        }

                        char[] charArray = bioJavaMSAColOrdered.toCharArray();

                        char temp = charArray[1];
                        charArray[1] = charArray[2];
                        charArray[2] = temp;

                        bioJavaMSAColOrdered = new String (charArray);
                    }

                    if (goodCol[1].equals(bioJavaMSAColOrdered)) {
                        bioJavaMSACount += 1;
                    }










//                    System.out.println(goodCol[0]);
//                    System.out.println(goodCol[1]);
//
//                    System.out.println(bioJavaMSA.getCompoundsAt(Integer.valueOf(goodCol[0]) + 1));
//                    System.out.println(bioJavaMSAColOrdered);
//                    System.out.println(MSA.getColumn(Integer.valueOf(goodCol[0])));

                }

                double msaCS = (double) msaCount / goodCols.size() * 100.0;
                double bioJavaMSACS = (double) bioJavaMSACount / goodCols.size() * 100.0;


                System.out.printf("MSA Column Score for file " + input + " is: %.2f%%%n ", msaCS);
                System.out.printf("BioJava MSA Column Score for file " + input + " is: %.2f%%%n ",  bioJavaMSACS);
                System.out.println();






                // Compare calculated alignment to reference alignment

            }

        } catch (Exception e) {
            // if any error occurs
            e.printStackTrace();
        }
    }

    public static ArrayList<String[]> getGoodCols(ArrayList<String[]> seqList) throws IOException{

        String seq = seqList.get(0)[1];
        ArrayList<String[]> goodCols = new ArrayList<String[]>();


        for (int i = 0; i < seq.length(); i ++) {
            String colChar = "";
            if (Character.isUpperCase(seq.charAt(i))) {
                for (int j=0; j < seqList.size(); j++){
                    colChar += seqList.get(j)[1].charAt(i);
                }
                String[] colPos = new String[2];
                colPos[0] = (String.valueOf(i));
                colPos[1] = colChar;

                goodCols.add(colPos);



            }
        }

        return goodCols;

    }

    public static ArrayList<String[]> getSeqs(String dir, String input) throws IOException{

        FileReader fr = new FileReader(dir + "/" +  input);

        BufferedReader br = new BufferedReader(fr);
        String[] seqs = new String[2];
        ArrayList<String[]> seqList = new ArrayList<String[]>();
        int count = 0;
        while (( seqs = getNextSequence(br)) != null){
            seqList.add(count, seqs);
            count++;
        }

        return seqList;
    }

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

    public static HashProfile getMSA(ArrayList<String[]> seqs, int multiGapOpen, int multiGapExtend){

        ArrayList<HashProfile> seqList = new ArrayList<HashProfile>();

        for (String[] seq: seqs){
            HashProfile profile = new HashProfile(seq[1]);
            seqList.add(profile);
        }

        HashProfile profile = alignProfilesRecursive(seqList.get(0), seqList.get(1), seqList, multiGapOpen, multiGapExtend, 1);


//        for (int i = 0; i + 1 < seqList.size(); i++){
//
//
//
//
//            HashProfile next = alignProfiles(seqList.get(i), seqList.get(i+1), multiGapOpen, multiGapExtend);
//            HashProfile again = alignProfiles(next, i + 1)
//
//        }


        return profile;
    }



//    public static HashProfile alignProfiles(HashProfile first, HashProfile second, int multiGapOpen, int multiGapExtend){
//        Alignment alignment = new Alignment(first, second, multiGapOpen, multiGapExtend, Blosum62Probs.getMatrix(), false);
//        return alignment.getUpdatedProfile();
//
//    }

    public static HashProfile alignProfilesRecursive (HashProfile first, HashProfile second, ArrayList<HashProfile> seqList, int multiGapOpen, int multiGapExtend, int counter){
        Alignment alignment = new Alignment(first, second, multiGapOpen, multiGapExtend, Blosum62Probs.getMatrix(), false);
        HashProfile profile = alignment.getUpdatedProfile();
        counter ++;
        if (counter < seqList.size()){
            return alignProfilesRecursive(profile, seqList.get(counter), seqList, multiGapOpen, multiGapExtend, counter);
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



//        // Create an R vector in the form of a string.
//        String javaVector = "c(1,2,3,4,5)";
//
//        // Start Rengine.
//        Rengine engine = new Rengine(new String[] { "--no-save" }, false, null);
//
//        // The vector that was created in JAVA context is stored in 'rVector' which is a variable in R context.
//        engine.eval("rVector=" + javaVector);
//
//        //Calculate MEAN of vector using R syntax.
//        engine.eval("meanVal=mean(rVector)");
//
//        //Retrieve MEAN value
//        double mean = engine.eval("meanVal").asDouble();
//
//        //Print output values
//        System.out.println("Mean of given vector is=" + mean);

