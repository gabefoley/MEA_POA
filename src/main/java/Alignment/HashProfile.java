package Alignment;

import sun.java2d.pipe.SpanShapeRenderer;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * HashMap implementation of a sequence profile
 *
 *
 */



public class HashProfile {


    private List<Map<Character,MutableInt>> profileArray;
    private List<String> sequences;

    /**
     * Constructor if using a single sequence
     * @param seq1 Sequence to build the profile from
     */
    public HashProfile(String seq1){

        this.profileArray = new ArrayList<Map<Character,MutableInt>>();
        this.sequences = new ArrayList<String>();
        this.sequences.add(seq1);


        for (int i = 0; i < seq1.length(); i++) {
            profileArray.add(i, new HashMap<Character, MutableInt>());

        }

        fillProfileArray(seq1);
    }

    /**
     * Constructor if adding a sequence to a profile
     * @param profile The existing profile
     * @param seq The sequence to add to the profile
     */
    public HashProfile(HashProfile profile, String seq){

        this.profileArray = profile.getProfileArray();
        this.sequences = profile.getSequences();
        this.sequences.add(seq);
        fillProfileArray(seq);

    }

    /**
     *  Constructor if joining two profiles together
     * @param profile1 The existing profile
     * @param profile2 The profile to add
     */

    public HashProfile(HashProfile profile1, HashProfile profile2){
        List<Map<Character,MutableInt>> newprofileArray = profile1.getProfileArray();
        this.profileArray = profile1.getProfileArray();
        this.sequences = profile1.getSequences();
        for (String seq: profile2.getSequences()){
            this.sequences.add(seq);
        }
        fillProfileArray(profile2);

    }


    public void fillProfileArray(String seq){


        for (int i = 0; i < seq.length(); i++) {
//            if (!(seq.charAt(i) == ('-'))) {
                MutableInt count = profileArray.get(i).get(seq.charAt(i));

                if (count == null) {
                    profileArray.get(i).put(seq.charAt(i), new MutableInt());

                } else {
                    count.increment();
                }
            }
//        }

    }

    public void fillProfileArray(HashProfile profile){
        for (int i = 0; i < profile.profileArray.size(); i++) {
            for (Character residue : profile.getProfileArray().get(i).keySet()) {
//                if (!residue.equals("-")) {
                    MutableInt count = profileArray.get(i).get(residue);
                    if (count == null) {
                        profileArray.get(i).put(residue, new MutableInt());

                    } else {
                        count.increment();
                    }
                }
            }
//        }

    }

    public void addGaps(List<Integer> gapPos){

        //TODO: Don't remove and recreate profileArray everytime we add gaps
        this.profileArray = new ArrayList<Map<Character,MutableInt>>();

        for (int i = 0; i < this.getSequences().size(); i++){
            String seq = this.getSequences().get(i);
            for (int pos : gapPos){
                seq = seq.substring(0, pos) + "-" + seq.substring(pos);
            }

            this.sequences.set(i, seq);


        }

        // Create a new profileArray the size of the first sequence
        for (int j = 0; j < this.getSequences().get(0).length(); j++) {
            profileArray.add(j, new HashMap<Character, MutableInt>());
        }

        // Fill in the profileArray with all the information from the seqs
        for (int i = 0; i < this.getSequences().size(); i++) {
            String seq = this.getSequences().get(i);
            this.fillProfileArray(seq);
        }

    }



    /**
     * Get a specific column from the profile
     * @param pos Position of the column to return
     * @return String representing the characters at the chosen column
     */
    public String getColumn(int pos){
        String columnList = "";
        for (String seq: this.getSequences()){

            if (pos <= seq.length()) {
                columnList += seq.charAt(pos);
            }
            else {
                columnList += "X";
            }

        }

        return columnList;
    }


    public List<Map<Character,MutableInt>> getProfileArray(){
        return profileArray;
    }

    public int getLength(){
        return profileArray.size();
    }

    public List<String> getSequences(){
        return this.sequences;
    }



    public String toString(){

        String seqOutput = "";
        for (String seq: this.getSequences()){
            seqOutput += seq + "\n";

        }
        // Remove the final new line
        return seqOutput.replaceAll("\\n+$", "");

    }

    public void printProfileArray(){

        for (int i = 0; i < profileArray.size(); i++){
            for (Character name: profileArray.get(i).keySet()){
                String key = name.toString();
                String value = profileArray.get(i).get(name).toString();
                String map = profileArray.get(i).toString();
                System.out.println("Column " + i + " Residue: " + key + " Count: " + value );

            }
        }
    }





}
