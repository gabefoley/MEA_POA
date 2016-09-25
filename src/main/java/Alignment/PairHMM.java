package Alignment;

import SubstitutionModels.SubstitutionMatrix;



import java.util.ArrayList;
import java.util.List;

/**
 * PairHMM for aligning sequences
 */
public class PairHMM {

    private double transitionMM, transitionMX, transitionMY, transitionXX, transitionYY, transitionXM, transitionYM;
    private HashProfile profile1, profile2;
    private double[][] vM;
    private double[][] vX;
    private double[][] vY;
    private String[][] tracebackM, tracebackX, tracebackY;
//    private double forwardProb, backwardProb;

    private SubstitutionMatrix subMatrix;
    private double tau;

    public PairHMM(String seq1, String seq2, double tau, double epsilon, double delta, SubstitutionMatrix subMatrix) {

        this(new HashProfile(seq1), new HashProfile(seq2), tau, epsilon, delta, subMatrix);
    }

    public PairHMM(HashProfile profile1, HashProfile profile2, double[] start, double[][] transition, double[][] emission, SubstitutionMatrix subMatrix){
        this.profile1 = profile1;
        this.profile2 = profile2;

        this.subMatrix = subMatrix;

        //TODO: Current BW allows for transition between X and Y
        this.transitionMM = transition[0][0];
        this.transitionMX = transition[0][1];
        this.transitionMY = transition[0][2];
        this.transitionXM = transition[1][0];
        this.transitionXX = transition[1][1];
        this.transitionYM = transition[2][0];
        this.transitionYY = transition[2][2];

    }


    public PairHMM(HashProfile profile1, HashProfile profile2, double tau, double epsilon, double delta, SubstitutionMatrix subMatrix) {
        this.tau = tau;

        this.profile1 = profile1;
        this.profile2 = profile2;

        this.subMatrix = subMatrix;

        this.transitionMM = 1 - (2 * delta) - tau;
        this.transitionMX = this.transitionMY = delta;
        this.transitionXX = this.transitionYY = epsilon;
        this.transitionXM = this.transitionYM = 1 - epsilon - tau;


        this.vM = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.vX = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.vY = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];

        this.tracebackM = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.tracebackX = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.tracebackY = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];

        vM[0][0] = 1;
    }

    // Fill out the match matrix

    public void fillVM(int i, int j, double[][] vM, double[][] vX, double[][] vY, String[][] tracebackM) {


        if ( i == 0 || j == 0){
            return;
        }


        double totalCount = getTotalCount(profile1, profile2, i, j);
        double totalScore = getTotalScore(profile1, profile2, i, j, subMatrix);


        double emissionM = totalScore / (totalCount);



        // Get the actual costs for transitioning for each of the states
        double currentTransitionMM = transitionMM * vM[i - 1][j - 1];
        double currentTransitionXM = transitionXM * vX[i - 1][j - 1];
        double currentTransitionYM = transitionYM * vY[i - 1][j - 1];

        // Work out the optimal cost and set the cell
        if (currentTransitionMM > currentTransitionXM && currentTransitionMM > currentTransitionYM) {
            vM[i][j] = currentTransitionMM * emissionM;
            tracebackM[i][j] = "M";
        } else if (currentTransitionXM > currentTransitionMM && currentTransitionXM > currentTransitionYM) {
            vM[i][j] = currentTransitionXM * emissionM;
            tracebackM[i][j] = "X";
        } else {
            vM[i][j] = currentTransitionYM * emissionM;
            tracebackM[i][j] = "Y";

        }


    }


    // Fill out the gap in X matrix
    public void fillVX(int i, int j, double[][] vM, double[][] vX, String[][] tracebackX) {

        if (j == 0){
            return;
        }

        //TODO: Remove magic number here and in similar places
        double emissionX = 0.25;



        double currentTransitionXX = transitionXX * vX[i][j - 1];
        double currentTransitionMX = transitionMX * vM[i][j - 1];


        if (currentTransitionXX > currentTransitionMX) {
            vX[i][j] = currentTransitionXX * emissionX;
            tracebackX[i][j] = "X";

        } else {
            vX[i][j] = currentTransitionMX * emissionX;
            tracebackX[i][j] = "M";

        }

    }


    // Fill out the gap in Y matrix
    public  void fillVY(int i, int j, double[][] vM, double[][] vY, String[][] tracebackY) {

        if (i == 0 ){
            return;
        }

        double emissionY = 0.25;

        double currentTransitionYY = transitionYY * vY[i - 1][j];
        double currentTransitionMY = transitionMY * vM[i - 1][j];

        if (currentTransitionYY > currentTransitionMY) {
            vY[i][j] = currentTransitionYY * emissionY;
            tracebackY[i][j] = "Y";

        } else {
            vY[i][j] = currentTransitionMY * emissionY;
            tracebackY[i][j] = "M";

        }

    }

    public HashProfile traceback(){


        String lastState;
        int i = profile1.getProfileArray().size();
        int j = profile2.getProfileArray().size();
        List<Integer> profile1Matches = new ArrayList<Integer>();
        List<Integer > profile2Matches = new ArrayList<Integer>();
        int curstrIdx = 0;
        int curnodeIdx = 0;



        while ((i > 0) && (j > 0)) {
            if ((vM[i][j] > vX[i][j]) && (vM[i][j] > vY[i][j])) {
                profile1Matches.add(0, i - 1);
                profile2Matches.add(0, j - 1);

                lastState = tracebackM[i][j];
                curstrIdx = j - 1;
                curnodeIdx = i - 1;
                i--;
                j--;

            }

            else if((vX[i][j]) > vM[i][j] && (vX[i][j]) > vY[i][j]){
                profile1Matches.add(0, -1);
                profile2Matches.add(0, j - 1);

                lastState = tracebackX[i][j];
                j--;
            }

            else {

                profile1Matches.add(0, i - 1);
                profile2Matches.add(0, -1);
                lastState = tracebackY[i][j];
                i--;
            }


            while ((i > 0) && (j > 0)) {
                if (lastState.equals("M")){
                    profile1Matches.add(0, i - 1);
                    profile2Matches.add(0, j - 1);
                    lastState = tracebackM[i][j];
                    curnodeIdx = i - 1;
                    curstrIdx = j - 1;

                    i--;
                    j--;
                }
                else if (lastState.equals("Y")){
//                    seq1Output = profile1.charAt(i-1) + seq1Output;
//                    seq2Output = "-" + seq2Output;
                    profile1Matches.add(0, i - 1);
                    profile2Matches.add(0, -1);
                    lastState = tracebackY[i][j];
                    i--;
                }

                else {

                    profile1Matches.add(0, -1);
                    profile2Matches.add(0, j - 1);
                    lastState = tracebackX[i][j];
                    j--;
                }
            }

            // Fill out the remaining indexes of each profile
            while (profile1Matches.get(0) > 0 && curnodeIdx > 0){
                profile1Matches.add(0, curnodeIdx - 1);
                profile2Matches.add(0, -1);
                curnodeIdx -= 1;

            }

            while (profile2Matches.get(0) > 0 && curstrIdx > 0) {
                profile2Matches.add(0, curstrIdx - 1);
                profile1Matches.add(0, -1);
                curstrIdx -= 1;

            }

        }


        //TODO: Make calling gaps a function
        List<Integer> gapPos = new ArrayList<Integer>();
        List<Integer> gapPos2 = new ArrayList<Integer>();


        for (int k = 0; k < profile1Matches.size(); k++){
            if (profile1Matches.get(k) == -1){
                gapPos.add(k);

            }
        }

        for (int k = 0; k < profile2Matches.size(); k++){
            if (profile2Matches.get(k) == -1){
                gapPos2.add(k);
            }

        }

        if (gapPos.size() > 0){
            profile1.addGaps(gapPos);

        }

        if (gapPos2.size() > 0) {
            profile2.addGaps(gapPos2);
        }

        return new HashProfile(profile1, profile2);


    }

    public double[][] forwardAlgorithm() {

        double[][] fM = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        double[][] fX = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        double[][] fY = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];

        fM[0][0] = 1;

        for (int i = 0; i <= profile1.getProfileArray().size(); i++) {
            for (int j = 0; j <= profile2.getProfileArray().size(); j++) {
                if (i!= 0 || j!=0) {
                    sumfM(i, j, fM, fX, fY);
                    sumfX(i, j, fM, fX);
                    sumfY(i, j, fM, fY);
                }
            }

        }

//        double forwardProb = tau * (fM[profile1.getProfileArray().size()][profile2.getProfileArray().size()] + fX[profile1.getProfileArray().size()][profile2.getProfileArray().size()] +
//                fY[profile1.getProfileArray().size()][profile2.getProfileArray().size()]);

        return fM;
    }

    public void sumfM(int i, int j, double[][] fM, double[][]fX, double[][] fY){

        if((i - 1 < 0) || (j - 1 < 0)){
            fM[i][j] = 0;
        }

        else {


            double totalCount = getTotalCount(profile1, profile2, i, j);
            double totalScore = getTotalScore(profile1, profile2, i, j, subMatrix);

            double emissionM = totalScore / (totalCount);

            double forwardMM = transitionMM * fM[i - 1][j - 1];
            double forwardXM = transitionXM * fX[i - 1][j - 1];
            double forwardYM = transitionYM * fY[i - 1][j - 1];

            fM[i][j] = emissionM * (forwardMM + forwardXM + forwardYM);
        }



    }

    public void sumfX(int i, int j, double[][] fM, double[][]fX){
        double emissionX = 0.25;


        // If we're in the first column
        if (j - 1 < 0) {
            fX[i][j] = 0;
        } else {


            double forwardXX = transitionXX * fX[i][j - 1];
            double forwardMX = transitionMX * fM[i][j - 1];

            fX[i][j] = emissionX * (forwardMX + forwardXX);

        }

    }

    public void sumfY(int i, int j, double[][] fM, double[][] fY){

        double emissionY = 0.25;



        if (i - 1 < 0) {
            fY[i][j] = 0;
        } else {


            double forwardYY = transitionYY * fY[i - 1][j];
            double forwardMY = transitionMY * fM[i - 1][j];

            fY[i][j] = emissionY * (forwardMY + forwardYY);

        }


    }
    //
    public  double[][] backwardAlgorithm() {

        double[][] bM = new double[profile1.getProfileArray().size() + 2][profile2.getProfileArray().size() + 2];
        double[][] bX = new double[profile1.getProfileArray().size() + 2][profile2.getProfileArray().size() + 2];
        double[][] bY = new double[profile1.getProfileArray().size() + 2][profile2.getProfileArray().size() + 2];

        bM[profile1.getProfileArray().size()][profile2.getProfileArray().size()] = tau;
        bX[profile1.getProfileArray().size()][profile2.getProfileArray().size()] = tau;
        bY[profile1.getProfileArray().size()][profile2.getProfileArray().size()] = tau;


        for (int i = profile1.getProfileArray().size(); i > 0; i--) {
            for (int j = profile2.getProfileArray().size(); j > 0 ; j--) {
                if (i!= profile1.getProfileArray().size() || j!= profile2.getProfileArray().size()) {
                    sumbM(i, j, bM, bX, bY);
                    sumbX(i, j, bM, bX);
                    sumbY(i, j, bM, bY);
                }
            }

        }


//        backwardProb = bM[1][1] + bX[1][1] + bY[1][1];


        return bM;


    }

    public void sumbM(int i, int j, double[][] bM, double[][] bX, double[][] bY){

        double emissionM = getEmission(profile1, profile2, i, j);

        double emissionX = 0.25;
        double emissionY = 0.25;




        double backwardMM = emissionM * transitionMM * bM[i + 1][j + 1];
        double backwardXM = emissionX * transitionMX * bX[i][j + 1];
        double backwardYM = emissionY * transitionMY * bY[i + 1][j];

        bM[i][j] = backwardMM + backwardXM + backwardYM;

    }

    public void sumbX(int i, int j, double[][] bM, double[][] bX){

        double emissionM = getEmission(profile1, profile2, i, j);



        double emissionX = 0.25;


        double backwardMX = emissionM * transitionXM * bM[i + 1][j + 1];
        double backwardXX = emissionX * transitionXX * bX[i][j + 1];

        bX[i][j] = backwardMX + backwardXX;



    }

    public void sumbY(int i, int j, double[][] bM, double[][] bY) {

        double emissionM = getEmission(profile1, profile2, i, j);

        double emissionY = 0.25;

        double backwardMY = emissionM * transitionYM * bM[i + 1][j + 1];
        double backwardYY = emissionY * transitionYY * bY[i + 1][j];

        bY[i][j] = backwardMY + backwardYY;

    }


    public HashProfile getViterbiAlignment(){

        for (int i = 0; i <= profile1.getProfileArray().size(); i++) {
            for (int j = 0; j <= profile2.getProfileArray().size(); j ++) {
                if (!(i == 0) || !(j == 0)) {
                    fillVM(i, j, vM, vX, vY, tracebackM);
                    fillVX(i, j, vM, vX, tracebackX);
                    fillVY(i, j, vM, vY, tracebackY);
                }
            }

        }

        return traceback();



    }

    public HashProfile getMEAAlignment(){
        double[][] fM = this.forwardAlgorithm();
        double[][] bM = this.backwardAlgorithm();
        double[][] pM = calcPosteriorMatrix(fM, bM);

        SubstitutionMatrix subMatrix = new SubstitutionMatrix(pM);

        Alignment alignment = new Alignment(profile1, profile2, 0, 0, subMatrix, true);

        return alignment.getUpdatedProfile();
    }

    public static double[][] calcPosteriorMatrix(double[][] fM, double[][]bM){

        double[][] pM = new double[fM.length][fM[0].length];

        for (int i = 0; i < fM.length; i++) {
            for (int j = 0; j < fM[0].length; j++) {
                pM[i][j] = fM[i][j] * bM[i][j];
            }


        }

        return pM;
    }


    public double getTotalCount(HashProfile profile1, HashProfile profile2, int i, int j){
        double profile1Count = 0;
        double profile2Count = 0;
        //TODO: Check calculation of profileCounts - maybe add it back into the loop
        //TODO: Check maths on summing profiles

        for (Character residue : profile1.getProfileArray().get(i-1).keySet()){
            profile1Count += profile1.getProfileArray().get(i-1).get(residue).getValue();
        }

        for (Character residue : profile2.getProfileArray().get(j-1).keySet()){
            profile2Count += profile2.getProfileArray().get(j-1).get(residue).getValue();
        }

        return profile1Count * profile2Count;
    }

    public double getTotalScore(HashProfile profile1, HashProfile profile2, int i, int j, SubstitutionMatrix subMatrix){
        double totalScore = 0;

        for (Character name : profile1.getProfileArray().get(i - 1).keySet()) {
            if (name != '-') {


                int profile1Value = profile1.getProfileArray().get(i - 1).get(name).getValue();

                for (Character name2 : profile2.getProfileArray().get(j - 1).keySet()) {
                    if (name2 != '-') {
                        int profile2Value = profile2.getProfileArray().get(j - 1).get(name2).getValue();


                        double matchScore = subMatrix.getDistance(name, name2);
                        totalScore += profile1Value * profile2Value * matchScore;

                    }
                }
            }

        }

        return totalScore;

    }

    public double getEmission(HashProfile profile1, HashProfile profile2, int i, int j){
        double emission;
        if (i > profile1.getProfileArray().size() - 1 || j > profile2.getProfileArray().size() - 1){
            emission = 0;
        }
        else {

            double totalCount = getTotalCount(profile1, profile2, i, j);

            double totalScore = getTotalScore(profile1, profile2, i, j, subMatrix);

            emission = totalScore / (totalCount);
        }

        return emission;
    }



}
