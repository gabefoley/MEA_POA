package Alignment;

import sun.java2d.pipe.SpanShapeRenderer;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by gabe on 15/08/2016.
 *
 *  # 0 1 2 3 4 5
 *  A
 *  G
 *  C
 *  T
 *
 */



public class HashProfile {


    private List<Map<Character,MutableInt>> profileArray;
    private List<String> sequences;

    public HashProfile(String seq1){
        this.profileArray = new ArrayList<Map<Character,MutableInt>>();
        this.sequences = new ArrayList<String>();
        this.sequences.add(seq1);
        for (int i = 0; i < seq1.length(); i++) {

            profileArray.add(i, new HashMap<Character, MutableInt>());

        }
        fillProfileArray(seq1);



    }


    public HashProfile(HashProfile profile, String seq){

        this.profileArray = profile.getProfileArray();
        this.sequences = profile.getSequences();
        this.sequences.add(seq);
        fillProfileArray(seq);

    }

    public HashProfile(HashProfile profile1, HashProfile profile2){
        this.profileArray = profile1.getProfileArray();
        this.sequences = profile1.getSequences();
        for (String seq: profile2.getSequences()){
            this.sequences.add(seq);
        }
        fillProfileArray(profile2);

    }

    public HashProfile(HashProfile profile1, HashProfile profile2, List<List<Integer>> matchesIndex){
//        for (Integer index: matchesIndex.get(0)){
//            if (index == null){
//                for (String seq: profile1.getSequences()){
//                    seq = seq.substring(0, index) + "-" + seq.substring(index);
//                }
//            }
//            System.out.println("Integer in 0 is " + index);
//
//        }
//
//        for (Integer index: matchesIndex.get(1)){
//            if (index == null){
//                for (String seq: profile2.getSequences()){
//                    seq = seq.substring(0, index) + "-" + seq.substring(index);
//                }
//            }
//            System.out.println("Integer in 1 is " + index);
//
//        }

    }

    public void fillProfileArray(String seq){


        for (int i = 0; i < seq.length(); i++) {
            MutableInt count = profileArray.get(i).get(seq.charAt(i));

            if (count == null) {
                profileArray.get(i).put(seq.charAt(i), new MutableInt());

            } else {
                count.increment();
            }
        }

    }

    public void fillProfileArray(HashProfile profile){
        for (int i = 0; i < profile.profileArray.size(); i++) {
            for (Character residue : profile.getProfileArray().get(i).keySet()) {
                if (!residue.equals("-")) {
                    MutableInt count = profileArray.get(i).get(residue);
                    if (count == null) {
                        profileArray.get(i).put(residue, new MutableInt());

                    } else {
                        count.increment();
                    }
                }
            }
        }

    }

    public void addGaps(List<Integer> gapPos){

        for (int pos: gapPos){
            System.out.println(pos);
        }

        for (int i = 0; i < this.getSequences().size(); i++){
            String seq = this.getSequences().get(i);
            for (int pos : gapPos){
                seq = seq.substring(0, pos) + "-" + seq.substring(pos);
            }

            this.sequences.set(i, seq);
        }



//        for (String seq : this.getSequences()) {
//            int index = this.getSequences().indexOf(seq);
//            for (int pos : gapPos){
//                seq = seq.substring(0, pos) + "-" + seq.substring(pos);
//            }
//            this.sequences.remove(index);
//            this.sequences.add(seq);
//
//        }

    }





    public List<Map<Character,MutableInt>> getProfileArray(){
        return this.profileArray;
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
//                System.out.println("Column " + i + " " + map );


            }
        }
    }





}
