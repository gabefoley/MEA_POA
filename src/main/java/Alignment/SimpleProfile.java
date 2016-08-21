package Alignment;

import sun.java2d.pipe.SpanShapeRenderer;

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


public class SimpleProfile {


    private int[][] profileMatrix;

    public SimpleProfile(String seq1){
        this.profileMatrix = new int[5][seq1.length()];
//        PairHMM.printMatrix(profileMatrix);

        fillProfileMatrix(seq1);

//        PairHMM.printMatrix(profileMatrix);



    }

    public SimpleProfile(String seq1, String seq2){
        this.profileMatrix = new int[5][seq1.length()];
        fillProfileMatrix(seq1);
        fillProfileMatrix(seq2);
//        PairHMM.printMatrix(profileMatrix);

    }

    public SimpleProfile(SimpleProfile profile, String seq){

        this.profileMatrix = profile.profileMatrix;
        fillProfileMatrix(seq);

    }

    public SimpleProfile(SimpleProfile profile1, SimpleProfile profile2){
        this.profileMatrix = profile1.profileMatrix;
        fillProfileMatrix(profile2);

    }

    public void fillProfileMatrix(String seq){
        for (int i = 0; i < seq.length(); i++){
            profileMatrix[getPositionFromResidue(seq.substring(i,i + 1))][i] +=1;
        }
    }

    public void fillProfileMatrix(SimpleProfile profile){

    }

    public int[][] getProfileMatrix(){
        return this.profileMatrix;
    }

    public int getPositionFromResidue(String residue){

        if (residue.equals("A")) {
            return 0;
        }

        else if (residue.equals("G")) {
            return 1;
        }
        else if (residue.equals("C")) {
            return 2;
        }

        else if (residue.equals("T")) {
            return 3;
        }

        else if(residue.equals("-")){
            return 4;
        }

        return -1 ;


    }

    public char getResidueFromPosition(int position){

        if (position == 0){
            return 'A';
        }

        else if (position == 1){
            return 'G';
        }
        else if (position == 2){
            return 'C';
        }
        else if (position == 3){
            return 'T';
        }

        return 'X';
    }

}
