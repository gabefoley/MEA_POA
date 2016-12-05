package Alignment;

import Alignment.Utilities.MatrixUtils;
import Alignment.Utilities.Sequence;
import SubstitutionModels.SubstitutionMatrix;
import dat.POGraph;


import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * PairHMM for aligning sequences
 */
public class POPairHMMUnderflow {

    private double transitionMM, transitionMX, transitionMY, transitionXX, transitionYY, transitionXM, transitionYM,
            transitionSM, transitionSX, transitionSY;
    private HashProfile profile1, profile2;
    private double[][] vM;
    private double[][] vX;
    private double[][] vY;
    private String[][] tracebackM, tracebackX, tracebackY;
    private String[] seqArray;
    private String seq1;
    private String seq2;
    private HashProfile[] profileArray;
    private double[][] emissions;
    private double emissionX;
    private double emissionY;
    private double fValue;
    private double bValue;
//    private double forwardProb, backwardProb;

    private SubstitutionMatrix subMatrix;
    private String type;
    //TODO: Make Tau dynamic
    private double tau = 0.1;

    public POPairHMMUnderflow(String seq1, String seq2, double tau, double epsilon, double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type) {

        this(new HashProfile(seq1), new HashProfile(seq2), tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
    }

    public POPairHMMUnderflow(HashProfile profile1, HashProfile profile2, double[] start, double[][] transition, double[][] emission, SubstitutionMatrix subMatrix, boolean runBaumWelch, String type){


        this.profile1 = profile1;
        this.profile2 = profile2;

        this.profileArray = new HashProfile[2];
        profileArray[0] = profile1;
        profileArray[1] = profile2;

        this.subMatrix = subMatrix;
        this.type = type;

        //TODO: Current BW allows for transition between X and Y
        this.transitionMM = transition[0][0];
        this.transitionMX = transition[0][1];
        this.transitionMY = transition[0][2];
        this.transitionXM = transition[1][0];
        this.transitionXX = transition[1][1];
        this.transitionYM = transition[2][0];
        this.transitionYY = transition[2][2];
        this.transitionSM = start[0];
        this.transitionSX = start[1];
        this.transitionSY = start[2];

    }

    public POPairHMMUnderflow(String [] seqArray, double[] start, double[][]transition, double[][]emission, SubstitutionMatrix subMatrix, boolean runBaumWelch, String type){

        this.profile1 = new HashProfile(seqArray[0]);
        profile2 = new HashProfile(seqArray[1]);

        if (runBaumWelch) {
            BaumWelchMulti bw = new BaumWelchMulti(seqArray, start, transition, emission, subMatrix, type);
//            return new PairHMM(profile1, profile2, bw.getStart(), bw.getTransition(), bw.getEmission(), subMatrix, false);

        }


        profileArray = new HashProfile[seqArray.length];

        // Turn all the sequences into profiles
        for (int i = 0; i < seqArray.length; i++){
            profileArray[i]=(new HashProfile(seqArray[i]));

        }

        this.subMatrix = subMatrix;
        this.type = type;

//        double[] start = bw.getStart();
//        double[][] transition = bw.getTransition();
        this.emissions = emission;

        this.transitionMM = transition[0][0];
        this.transitionMX = transition[0][1];
        this.transitionMY = transition[0][2];
        this.transitionXM = transition[1][0];
        this.transitionXX = transition[1][1];
        this.transitionYM = transition[2][0];
        this.transitionYY = transition[2][2];
        this.transitionSM = start[0];
        this.transitionSX = start[1];
        this.transitionSY = start[2];

        this.vM = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.vX = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.vY = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];

        this.tracebackM = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.tracebackX = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.tracebackY = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];


        for (int i = 0; i < profile1.getProfileArray().size() + 1; i++) {
            vX[i][0] = -1;
            vM[i][0] = -1;
        }

        for (int i = 0; i < profile2.getProfileArray().size() + 1; i++) {
            vY[0][i] = -1;
            vM[0][i] = -1;

        }

        vM[0][0] = Math.log(1);




    }

    public POPairHMMUnderflow(Sequence[] seqArray, double[] start, double[][]transition, double[][]emission, SubstitutionMatrix subMatrix, boolean runBaumWelch, String type){

//        this.profile1 = new HashProfile(seqArray[0]);
//        this.profile2 = new HashProfile(seqArray[1]);

        if (runBaumWelch) {
//            BaumWelchMulti bw = new BaumWelchMulti(seqArray, start, transition, emission, "protein");
//            return new PairHMM(profile1, profile2, bw.getStart(), bw.getTransition(), bw.getEmission(), subMatrix, false);

        }

        profile1 = new HashProfile(seqArray[0]);
        profile2 = new HashProfile(seqArray[1]);

        profileArray = new HashProfile[seqArray.length];

        // Turn all the sequences into profiles
        for (int i = 0; i < seqArray.length; i++){
            profileArray[i]=(new HashProfile(seqArray[i]));

        }


        this.subMatrix = subMatrix;
        this.type = type;

//        double[] start = bw.getStart();
//        double[][] transition = bw.getTransition();
        this.emissions = emission;

        this.transitionMM = transition[0][0];
        this.transitionMX = transition[0][1];
        this.transitionMY = transition[0][2];
        this.transitionXM = transition[1][0];
        this.transitionXX = transition[1][1];
        this.transitionYM = transition[2][0];
        this.transitionYY = transition[2][2];
        this.transitionSM = start[0];
        this.transitionSX = start[1];
        this.transitionSY = start[2];

        this.vM = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.vX = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.vY = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];

        this.tracebackM = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.tracebackX = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.tracebackY = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];

        for (int i = 0; i < profile1.getProfileArray().size() + 1; i++) {
            vX[i][0] = -1;
            vM[i][0] = -1;
        }

        for (int i = 0; i < profile2.getProfileArray().size() + 1; i++) {
            vY[0][i] = -1;
            vM[0][i] = -1;

        }

        vM[0][0] = Math.log(1);
    }

    public POPairHMMUnderflow(String[] seqArray, SubstitutionMatrix subMatrix, String type){
//        BaumWelch bw = new BaumWelch(seqArray, type);
        BaumWelchMulti bw = new BaumWelchMulti(seqArray, type);

//        for (double value: bw.getStart()){
//            System.out.println(value);
//
//        }
//        MatrixUtils.printMatrix(bw.getEmission());
//        MatrixUtils.printMatrix(bw.getTransition());
//
//        for (double value: bw.getStart()){
//            System.out.println(value);
//
//        }
//        MatrixUtils.printMatrix(bw.getEmission());
//        MatrixUtils.printMatrix(bw.getTransition());


        profile1 = new HashProfile(seqArray[0]);
        profile2 = new HashProfile(seqArray[1]);

        profileArray = new HashProfile[seqArray.length];

        // Turn all the sequences into profiles
        for (int i = 0; i < seqArray.length; i++){
            profileArray[i]=(new HashProfile(seqArray[i]));

        }

        this.subMatrix = subMatrix;
        this.type = type;

        double[] start = bw.getStart();
        double[][] transition = bw.getTransition();
        emissions = bw.getEmission();

        this.transitionMM = transition[0][0];
        this.transitionMX = transition[0][1];
        this.transitionMY = transition[0][2];
        this.transitionXM = transition[1][0];
        this.transitionXX = transition[1][1];
        this.transitionYM = transition[2][0];
        this.transitionYY = transition[2][2];
        this.transitionSM = start[0];
        this.transitionSX = start[1];
        this.transitionSY = start[2];

        this.vM = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.vX = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.vY = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];

        this.tracebackM = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.tracebackX = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        this.tracebackY = new String[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];

        vM[0][0] = Math.log(1);


//        PairHMM bwPairHMM = new PairHMM(profile1, profile2, bw.getStart(), bw.getTransition(), bw.getEmission(), subMatrix, false);





    }


    public POPairHMMUnderflow(String[] seqArray, double tau, double epsilon, double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type) {

//
        this.seqArray = seqArray;
        this.seq1 = seqArray[0];
        this.seq2 = seqArray[1];
        this.tau = tau;

        this.profile1 = profile1;
        this.profile2 = profile2;


        this.subMatrix = subMatrix;
        this.type = type;

        this.transitionMM = 1 - (2 * delta) - tau;
        this.transitionMX = this.transitionMY = delta;
        this.transitionXX = this.transitionYY = epsilon;
        this.transitionXM = this.transitionYM = 1 - epsilon - tau;


        this.emissionX = emissionX;
        this.emissionY = emissionY;

        createMatrices();



    }

    public POPairHMMUnderflow(Sequence[] seqArray, double tau, double epsilon, double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type) {

        profile1 = new HashProfile(seqArray[0]);
        profile2 = new HashProfile(seqArray[1]);

        profileArray = new HashProfile[seqArray.length];

        // Turn all the sequences into profiles
        for (int i = 0; i < seqArray.length; i++){
            profileArray[i]=(new HashProfile(seqArray[i]));

        }
        this.tau = tau;

        this.profile1 = profile1;
        this.profile2 = profile2;


        this.subMatrix = subMatrix;
        this.type = type;

        this.transitionMM = 1 - (2 * delta) - tau;
        this.transitionMX = this.transitionMY = delta;
        this.transitionXX = this.transitionYY = epsilon;
        this.transitionXM = this.transitionYM = 1 - epsilon - tau;

        this.emissionX = emissionX;
        this.emissionY = emissionY;

        createMatrices();



    }

    public POPairHMMUnderflow(HashProfile profile1, HashProfile profile2, double tau, double epsilon, double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type) {
        this.tau = tau;

        this.profile1 = profile1;
        this.profile2 = profile2;


        this.profileArray = new HashProfile[2];
        profileArray[0] = profile1;
        profileArray[1] = profile2;

        this.subMatrix = subMatrix;
        this.type = type;

        this.transitionMM = 1 - (2 * delta) - tau;
        this.transitionMX = this.transitionMY = delta;
        this.transitionXX = this.transitionYY = epsilon;
        this.transitionXM = this.transitionYM = 1 - epsilon - tau;

        this.emissionX = emissionX;
        this.emissionY = emissionY;

        createMatrices();



    }

    public void createMatrices(){
        this.vM = new double[seq1.length() + 1][seq2.length() + 1];
        this.vX = new double[seq1.length() + 1][seq2.length() + 1];
        this.vY = new double[seq1.length() + 1][seq2.length() + 1];

        this.tracebackM = new String[seq1.length() + 1][seq2.length() + 1];
        this.tracebackX = new String[seq1.length() + 1][seq2.length() + 1];
        this.tracebackY = new String[seq1.length() + 1][seq2.length() + 1];

        for (int i = 0; i < seq1.length() + 1; i++) {
            vX[i][0] = Double.NEGATIVE_INFINITY;
            vM[i][0] = Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < seq2.length() + 1; i++) {
            vY[0][i] = Double.NEGATIVE_INFINITY;
            vM[0][i] = Double.NEGATIVE_INFINITY;

        }

        this.vM[0][0] = Math.log(1);
    }

    /**
     *
     * @param i
     * @param j
     * @param vM
     * @param vX
     * @param vY
     * @param tracebackM
     */
    public void fillVM(int i, int j, double[][] vM, double[][] vX, double[][] vY, String[][] tracebackM) {

        if (i == 0){

        }


        if ( i == 0 || j == 0){
            return;
        }


//        double totalCount = getTotalCount(profile1, profile2, i, j);
//        double totalScore = getTotalScore(profile1, profile2, i, j, subMatrix);
//
//
//        double emissionM = totalScore / (totalCount);

        double emissionM = Math.log(getEmission(profile1, profile2, i, j));



        // Get the actual costs for transitioning for each of the states

        //TODO: Change these values -- CHANGED
        double currentTransitionMM = (vM[i-1][j-1] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionMM) + vM[i - 1][j - 1];
        double currentTransitionXM = (vX[i-1][j-1] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionXM) + vX[i - 1][j - 1];
        double currentTransitionYM = (vY[i-1][j-1] == -1) ? Double.NEGATIVE_INFINITY  : Math.log(transitionYM) + vY[i - 1][j - 1];

        // Work out the optimal cost and set the cell
        if (currentTransitionMM >= currentTransitionXM && currentTransitionMM >= currentTransitionYM) {
            vM[i][j] = currentTransitionMM + emissionM;
            tracebackM[i][j] = "M";
        } else if (currentTransitionXM >= currentTransitionMM && currentTransitionXM >= currentTransitionYM) {
            vM[i][j] = currentTransitionXM + emissionM;
            tracebackM[i][j] = "X";
        } else {
            vM[i][j] = currentTransitionYM + emissionM;
            tracebackM[i][j] = "Y";

        }


    }

    /**
     *
     * @param i
     * @param j
     * @param vM
     * @param vX
     * @param tracebackX
     */
    // Fill out the gap in X matrix
    public void fillVX(int i, int j, double[][] vM, double[][] vX, String[][] tracebackX) {

        if (j == 0){
            return;
        }

        if (emissions != null) {

            emissionX = getIndelEmission(profile2, j, 1);

//            for (Character character : profile2.getProfileArray().get(j - 1).keySet()) {
//                emissionX = emissions[1][MatrixUtils.returnIndex(character)];
//            }
        }

//        double totalCount = getTotalEmissionCount(profile1, profile2, i, j);
//        double totalScore = getTotalEmissionScore(profile1, profile2, i, j, subMatrix);
//
//        double emissionX = totalScore / (totalCount);

//        double emissionX = getIndelEmission(profile1, j, 1);



        //TODO: Change these values -- CHANGED -- SIGNED

        double currentTransitionXX = (vX[i][j-1] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionXX) + vX[i][j - 1];
        double currentTransitionMX = (vM[i][j-1] == -1) ? Double.NEGATIVE_INFINITY : Math.log(transitionMX) + vM[i][j - 1];


        if (currentTransitionXX >= currentTransitionMX) {
            vX[i][j] = currentTransitionXX + Math.log(emissionX);
            tracebackX[i][j] = "X";

        } else {
            vX[i][j] = currentTransitionMX + Math.log(emissionX);
            tracebackX[i][j] = "M";

        }

    }

    /**
     *
     * @param i
     * @param j
     * @param vM
     * @param vY
     * @param tracebackY
     */
    // Fill out the gap in Y matrix
    public  void fillVY(int i, int j, double[][] vM, double[][] vY, String[][] tracebackY) {

        if (i == 0){
            return;
        }

        if (emissions != null) {

            emissionY = getIndelEmission(profile1, i, 2);

//            for (Character character : profile1.getProfileArray().get(i- 1).keySet()) {
//                if (i == 6 && j == 6){
//                    System.out.println("sixes");
//                    System.out.println(character);
//                }
//                emissionY = emissions[2][MatrixUtils.returnIndex(character)];
//            }
        }


//        double emissionY = 0.25;
        //TODO: Change these values --CHANGED

        double currentTransitionYY = (vY[i-1][j] == -1) ? Double.NEGATIVE_INFINITY:  Math.log(transitionYY) + vY[i - 1][j];
        double currentTransitionMY = (vM[i-1][j] == -1) ? Double.NEGATIVE_INFINITY: Math.log(transitionMY) + vM[i - 1][j];

        if (currentTransitionYY >= currentTransitionMY) {
            vY[i][j] = currentTransitionYY + Math.log(emissionY);
            tracebackY[i][j] = "Y";

        } else {
            vY[i][j] = currentTransitionMY + Math.log(emissionY);
            tracebackY[i][j] = "M";

        }

    }

    /**
     *
     * @return
     */
    public POGraph traceback(){

        long tracestartTime = System.nanoTime();

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

            } else if ((vX[i][j]) > vM[i][j] && (vX[i][j]) > vY[i][j]) {
                profile1Matches.add(0, -1);
                profile2Matches.add(0, j - 1);

                lastState = tracebackX[i][j];
                j--;
            } else {

                profile1Matches.add(0, i - 1);
                profile2Matches.add(0, -1);
                lastState = tracebackY[i][j];
                i--;
            }


            while ((i > 0) && (j > 0)) {
                if (lastState.equals("M")) {
                    profile1Matches.add(0, i - 1);
                    profile2Matches.add(0, j - 1);
                    lastState = tracebackM[i][j];
                    curnodeIdx = i - 1;
                    curstrIdx = j - 1;

                    i--;
                    j--;
                } else if (lastState.equals("Y")) {
//                    seq1Output = profile1.charAt(i-1) + seq1Output;
//                    seq2Output = "-" + seq2Output;
                    profile1Matches.add(0, i - 1);
                    profile2Matches.add(0, -1);
                    lastState = tracebackY[i][j];
                    i--;
                } else {

                    profile1Matches.add(0, -1);
                    profile2Matches.add(0, j - 1);
                    lastState = tracebackX[i][j];
                    j--;
                }
            }
//            curnodeIdx -=1;
//            curstrIdx -=1;

        }
        // Fill out the remaining indexes of each profile
        while (profile1Matches.get(0) > 0 || i > 0){
            curnodeIdx = curnodeIdx == 0 ? 1 : curnodeIdx;

            profile1Matches.add(0, i - 1);
            profile2Matches.add(0, -1);
            i -= 1;

        }

        while (profile2Matches.get(0) > 0 ||j > 0) {
            curstrIdx = curstrIdx == 0 ? 1 : curstrIdx;
            profile2Matches.add(0, j - 1);
            profile1Matches.add(0, -1);
            j -= 1;

        }

        long traceendTime = System.nanoTime();
        long traceduration = (traceendTime - tracestartTime);

//        System.out.println("Do traceback " + TimeUnit.SECONDS.convert(traceduration, TimeUnit.NANOSECONDS)
//        );

//        }

        long startTime = System.nanoTime();




        addGapsToProfiles(profile1Matches, profile2Matches);

        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
//        System.out.println("Add gaps to profiles  " + TimeUnit.SECONDS.convert(duration, TimeUnit.NANOSECONDS)
//        );

        return new POGraph(seq1, seq2);



    }

    /**
     *
     * @return
     */
    public double[][] forwardAlgorithm() {

        double[][] fM = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        double[][] fX = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];
        double[][] fY = new double[profile1.getProfileArray().size() + 1][profile2.getProfileArray().size() + 1];

//        createMatrices();

        for (int i = 0; i < profile1.getProfileArray().size() + 1; i++) {
            fX[i][0] = Double.NEGATIVE_INFINITY;
            fM[i][0] = Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < profile2.getProfileArray().size() + 1; i++) {
            fY[0][i] = Double.NEGATIVE_INFINITY;
            fM[0][i] = Double.NEGATIVE_INFINITY;

        }
        fM[0][0] = Math.log(1);

        for (int i = 0; i <= profile1.getProfileArray().size(); i++) {
            for (int j = 0; j <= profile2.getProfileArray().size(); j++) {
                if (i!= 0 || j!=0) {
                    sumfM(i, j, fM, fX, fY);
                    sumfX(i, j, fM, fX);
                    sumfY(i, j, fM, fY);
                }
            }

        }

        fValue = log_add(fM[profile1.getLength()][profile2.getLength()], fX[profile1.getLength()][profile2.getLength()]);
        fValue = log_add(fValue, fY[profile1.getLength()][profile2.getLength()]);
        fValue += Math.log(tau);
//        System.out.println("fValue " + fValue);

//        double forwardProb = tau * (fM[profile1.getProfileArray().size()][profile2.getProfileArray().size()] + fX[profile1.getProfileArray().size()][profile2.getProfileArray().size()] +
//                fY[profile1.getProfileArray().size()][profile2.getProfileArray().size()]);

        MatrixUtils.printMatrix(fM);
        MatrixUtils.printMatrix(fX);
        MatrixUtils.printMatrix(fY);
        return fM;
    }

    /**
     *
     * @param i
     * @param j
     * @param fM
     * @param fX
     * @param fY
     */
    public void sumfM(int i, int j, double[][] fM, double[][]fX, double[][] fY){


        if((i - 1 < 0) || (j - 1 < 0)){
//            fM[i][j] = 0;
        }

        else {

            double emissionM = getEmission(profile1, profile2, i, j);
            emissionM = emissionM ;





//            double totalCount = getTotalCount(profile1, profile2, i, j);
//            double totalScore = getTotalScore(profile1, profile2, i, j, subMatrix);
//
////            double emissionM = totalScore / (totalCount);
            //TODO: Change these values -- CHANGED


            double forwardMM = (fM[i-1][j-1] == Double.NEGATIVE_INFINITY) ? Double.NEGATIVE_INFINITY : Math.log(transitionMM) + fM[i - 1][j - 1];
            double forwardXM = (fX[i-1][j-1] == Double.NEGATIVE_INFINITY) ? Double.NEGATIVE_INFINITY : Math.log(transitionXM) + fX[i - 1][j - 1];
            double forwardYM = (fY[i-1][j-1] == Double.NEGATIVE_INFINITY) ? Double.NEGATIVE_INFINITY : Math.log(transitionYM) + fY[i - 1][j - 1];

            double transitionProbs = log_add(forwardXM, forwardYM);

            double fi = Math.log(emissionM);
            double f2 = log_add(forwardMM, transitionProbs);



            fM[i][j] = Math.log(emissionM) + log_add(forwardMM, transitionProbs);

            // If we're transitioning from the Start state to Match state.
//            if (i == 1 && j == 1){
//                double startProb = Math.log(transitionSM) + Math.log(emissionM) + fM[i-1][j-1];
////                fM[i][j] = log_add(fM[i][j], startProb);
//                fM[i][j] = startProb;
//
//            }
        }



    }

    /**
     *
     * @param i
     * @param j
     * @param fM
     * @param fX
     */
    public void sumfX(int i, int j, double[][] fM, double[][]fX){
//        double tempEmissionX = emissionX;


        if (j - 1 < 0){

//        }
//        else if (j - 1 < 0) {
////            fX[i][j] = 0;
        } else {

            if (emissions != null) {

                emissionX = getIndelEmission(profile2, j, 1);


//                for (Character character : profile1.getProfileArray().get(j - 1).keySet()) {
//                    if (MatrixUtils.returnIndex(character) == -1){
//                        System.out.println("got here!!");
//                    }
//                    double emissionX = emissions[1][MatrixUtils.returnIndex(character)];
//                }

            }


            //TODO: Change these values -- CHANGED

            double forwardMX = Math.log(transitionMX) + fM[i][j-1];
            double forwardXX = Math.log(transitionXX) + fX[i][j-1];

            if (fM[i][j-1] == Double.NEGATIVE_INFINITY){
                forwardMX = Double.NEGATIVE_INFINITY;
            }


            if (fX[i][j-1] == Double.NEGATIVE_INFINITY){
                forwardXX = Double.NEGATIVE_INFINITY;
            }

            if (fM[i][j-1] == Double.NEGATIVE_INFINITY && fX[i][j-1] == Double.NEGATIVE_INFINITY){
                fX[i][j] = Double.NEGATIVE_INFINITY;
//                tempEmissionX = 0;
            } else {

                fX[i][j] = Math.log(emissionX) + log_add(forwardMX, forwardXX);

            }

//            if (i == 1 && j == 0){
//                double startProb = Math.log(transitionSX) + Math.log(emissionX) + fM[i][j-1];
////                fX[i][j] = log_add(fX[i][j], startProb);
////                fX[i][j] = fX[i][j] + startProb;
//                fX[i][j] = startProb;
//                System.out.println("Start Prob: " + startProb);
//
//            }




        }

    }

    /**
     *
     * @param i
     * @param j
     * @param fM
     * @param fY
     */
    public void sumfY(int i, int j, double[][] fM, double[][] fY){

//        double tempEmission = emissionY;

        if (i - 1 < 0){

//        }
//
//
//        else if (i - 1 < 0) {
//            fY[i][j] = 0;
        } else {

            if (emissions != null) {

                emissionY = getIndelEmission(profile1, i, 2);


//                for (Character character : profile2.getProfileArray().get(i - 1).keySet()) {
//                    emissionY = emissions[2][MatrixUtils.returnIndex(character)];
//                }

            }

            //TODO: Change these values -- CHANGED


            double forwardYY = Math.log(transitionYY) + fY[i-1][j];
            double forwardMY = Math.log(transitionMY) + fM[i-1][j];


            if (fY[i-1][j] == Double.NEGATIVE_INFINITY){
                forwardYY = Double.NEGATIVE_INFINITY;
            }
            if (fM[i-1][j] == Double.NEGATIVE_INFINITY){
                forwardMY = Double.NEGATIVE_INFINITY;
            }

            if (fY[i-1][j] == Double.NEGATIVE_INFINITY  && fM[i-1][j] == Double.NEGATIVE_INFINITY){
                fY[i][j] = Double.NEGATIVE_INFINITY;
//                tempEmission = 0;
            } else {
                fY[i][j] = Math.log(emissionY) + log_add(forwardMY, forwardYY);

            }

//            if (i == 0 && j == 1){
//                double startProb = Math.log(transitionSY) + Math.log(emissionY) + fM[i-1][j];
//                fY[i][j] = startProb;
////                fY[i][j] = log_add(fY[i][j], startProb);
////                System.out.println("Start Prob: " + startProb);
//
//            }


        }


    }

    /**
     *
     * @return
     */
    public  double[][] backwardAlgorithm() {

        double[][] bM = new double[profile1.getProfileArray().size() + 2][profile2.getProfileArray().size() + 2];
        double[][] bX = new double[profile1.getProfileArray().size() + 2][profile2.getProfileArray().size() + 2];
        double[][] bY = new double[profile1.getProfileArray().size() + 2][profile2.getProfileArray().size() + 2];

        for (int i = 0; i < profile1.getProfileArray().size() + 1; i++) {
            bX[i][profile2.getProfileArray().size() + 1] = Double.NEGATIVE_INFINITY;
            bM[i][profile2.getProfileArray().size() + 1] = Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < profile2.getProfileArray().size() + 2; i++) {
            bY[profile1.getProfileArray().size() + 1][i] = Double.NEGATIVE_INFINITY;
            bM[profile1.getProfileArray().size() + 1][i] = Double.NEGATIVE_INFINITY;

        }
        bM[profile1.getProfileArray().size() ][profile2.getProfileArray().size() ] = Math.log(tau);
        bX[profile1.getProfileArray().size() ][profile2.getProfileArray().size() ] = Math.log(tau);
        bY[profile1.getProfileArray().size() ][profile2.getProfileArray().size() ] = Math.log(tau);

//        bM[profile1.getLength()][profile2.getLength()] = Math.log(0.1);
//        bX[profile1.getLength()][profile2.getLength()] = Math.log(0.1);
//        bY[profile1.getLength()][profile2.getLength()] = Math.log(0.1);


        for (int i = profile1.getProfileArray().size()  ; i > 0; i--) {
            for (int j = profile2.getProfileArray().size() ; j > 0; j--) {
                if (!(i == profile1.getLength() && j == profile2.getLength())) {

                    sumbM(i, j, bM, bX, bY);
                    sumbX(i, j, bM, bX);
                    sumbY(i, j, bM, bY);
                }

            }
        }


//        backwardProb = bM[1][1] + bX[1][1] + bY[1][1];
//            System.out.println("Backwards Matrix");
        System.out.println("backwards:");
        MatrixUtils.printMatrix(bM);
        MatrixUtils.printMatrix(bX);
        MatrixUtils.printMatrix(bY);


        return bM;

    }


    /**
     *
     * @param i
     * @param j
     * @param bM
     * @param bX
     * @param bY
     */
    public void sumbM(int i, int j, double[][] bM, double[][] bX, double[][] bY) {
        double emissionM = getEmission(profile1, profile2, i + 1, j + 1);

        if (i == profile1.getLength() && j == profile2.getLength()) {
            bM[i][j] = tau;


        } else {


//        double emissionX = 0.25;
//        double emissionY = 0.25;

            if (emissions != null) {

                emissionX = getIndelEmission(profile2, j, 1);
                emissionY = getIndelEmission(profile1, i, 2);


//            for (Character character : profile1.getProfileArray().get(i - 1).keySet()) {
//                emissionX = emissions[1][MatrixUtils.returnIndex(character)];
//            }
//
//            for (Character character : profile2.getProfileArray().get(i - 1).keySet()) {
//                emissionY = emissions[2][MatrixUtils.returnIndex(character)];
//            }

            }


            //TODO: Change these values --CHANGED

            double backwardMM = Math.log(emissionM) + Math.log(transitionMM) + bM[i + 1][j + 1];
            double backwardXM = Math.log(emissionX) + Math.log(transitionMX) + bX[i][j + 1];
            double backwardYM = Math.log(emissionY) + Math.log(transitionMY) + bY[i + 1][j];

            double transitionProbs = log_add(backwardXM, backwardYM);

            bM[i][j] = log_add(backwardMM, transitionProbs);
        }
    }

    /**
     *
     * @param i
     * @param j
     * @param bM
     * @param bX
     */
    public void sumbX(int i, int j, double[][] bM, double[][] bX){

        double emissionM = getEmission(profile1, profile2, i + 1, j + 1);



//        double emissionX = 0.25;

        if (emissions != null) {

            emissionX = getIndelEmission(profile2, j, 1);

//            for (Character character : profile1.getProfileArray().get(i - 1).keySet()) {
//                emissionX = emissions[1][MatrixUtils.returnIndex(character)];
//            }

        }

        //TODO: Change these values -- CHANGED

        double backwardMX = Math.log(emissionM) + Math.log(transitionXM) + bM[i + 1][j + 1];
        double backwardXX = Math.log(emissionX) + Math.log(transitionXX) + bX[i][j + 1];

        bX[i][j] = log_add(backwardMX, backwardXX);



    }

    /**
     *
     * @param i
     * @param j
     * @param bM
     * @param bY
     */
    public void sumbY(int i, int j, double[][] bM, double[][] bY) {

        double emissionM = getEmission(profile1, profile2, i + 1, j + 1);

//        double emissionY = 0.25;

        if (emissions != null) {

            emissionY = getIndelEmission(profile1, i, 2);


//            for (Character character : profile2.getProfileArray().get(i - 1).keySet()) {
//                emissionY = emissions[2][MatrixUtils.returnIndex(character)];
//            }

        }
        //TODO: Change these values --CHANGED

        double backwardMY = Math.log(emissionM) + Math.log(transitionYM) + bM[i + 1][j + 1];
        double backwardYY = Math.log(emissionY) + Math.log(transitionYY) + bY[i + 1][j];

        bY[i][j] = log_add(backwardMY, backwardYY);

    }

    /**
     *
     * @return
     */
    public POGraph getViterbiAlignment(){

        POGraph alignment = calculateViterbiAlignment(profileArray[0], profileArray[1], 1);

        return alignment;


    }

    public POGraph calculateViterbiAlignment(HashProfile profile1, HashProfile profile2, int counter){

        long startTime = System.nanoTime();


        for (int i = 0; i <= seq1.length(); i++) {
            for (int j = 0; j <= seq2.length(); j ++) {
                if (!(i == 0) || !(j == 0)) {

                    fillVM(i, j, vM, vX, vY, tracebackM);
                    fillVX(i, j, vM, vX, tracebackX);
                    fillVY(i, j, vM, vY, tracebackY);
                }
            }

        }

        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
//        System.out.println("Fill out matrix " + TimeUnit.SECONDS.convert(duration, TimeUnit.NANOSECONDS)
//        ); //divide by 1000000 to get milliseconds.


//        MatrixUtils.printMatrix(this.vM);
//        MatrixUtils.printMatrix(this.vX);
//        MatrixUtils.printMatrix(this.vY);

        POGraph newProfile = traceback();

//        System.out.println("Profile is " + newProfile);

//        counter ++;
//        if (counter < seqArray.length){
////            System.out.println("counter is " + counter);
//
//
//
//            // Update the references to the profiles to align
//            this.profile1 = newProfile;
////            System.out.println("Profile is " + newProfile);
//
//            this.profile2 = profileArray[counter];
//            createMatrices();
//            return calculateViterbiAlignment(newProfile, profileArray[counter], counter);
//        }
//
//        else {
////            System.out.println("got here");

            return newProfile;
//        }

//        return traceback();

    }


    public HashProfile getMEAAlignment(int counter){
        return getMEAAlignment(profileArray[0], profileArray[1], counter);
    }

    /**
     *
     * @return
     */
    public HashProfile getMEAAlignment(HashProfile profile1, HashProfile profile2, int counter){

        Alignment alignment = calculateMEAAlignment(profile1, profile2, counter);
//        System.out.println("PROFILE NOW IS: " + alignment.getUpdatedProfile());

        counter ++;
        if (counter < profileArray.length) {

//            System.out.println("Counter is " + counter);

            this.profile1 = alignment.getUpdatedProfile();
            this.profile2 = profileArray[counter];
//            createMatrices();


            return getMEAAlignment(this.profile1, this.profile2, counter);
//            return calculateMEAAlignment(alignment.getUpdatedProfile(), profileArray[counter], counter).getUpdatedProfile();


        }

        else {

            return alignment.getUpdatedProfile();

        }


    }

    public Alignment calculateMEAAlignment(HashProfile profile1, HashProfile profile2, int counter){
        double[][] fM = this.forwardAlgorithm();
        double[][] bM = this.backwardAlgorithm();
        double[][] pM = calcPosteriorMatrix(fM, bM);


//        MatrixUtils.printMatrix(fM);


        SubstitutionMatrix subMatrix = new SubstitutionMatrix(pM);

        Alignment alignment = new Alignment(profile1, profile2, 0, 0, subMatrix, true);

        return alignment;

//        counter ++;
//        if (counter < profileArray.length) {
//            System.out.println("Counter is " + counter);
//
//            this.profile1 = alignment.getUpdatedProfile();
//            this.profile2 = profileArray[counter];
//            createMatrices();
//            return calculateMEAAlignment(alignment.getUpdatedProfile(), profileArray[counter], counter);
//
//
//        }
//
//        else {
//
//            return alignment.getUpdatedProfile();
//
//        }
    }

    /**
     *
     * @param fM
     * @param bM
     * @return
     */
    public  double[][] calcPosteriorMatrix(double[][] fM, double[][]bM){

        double[][] pM = new double[fM.length][fM[0].length];


        //TODO: Change these values -- CHANGED BUT CAREFUL HERE


        for (int i = 0; i < fM.length; i++) {
            for (int j = 0; j < fM[0].length; j++) {
//                System.out.println(i + " " + j + " : " + fM[i][j] + " & " + bM[i][j]);
//                pM[i][j] = fM[i][j] + bM[i][j];
                pM[i][j] = Math.pow(Math.E, 0.05 * (fM[i][j] + bM[i][j]));

            }


        }
//        MatrixUtils.printMatrix(fM);
//        MatrixUtils.printMatrix(bM);
//        MatrixUtils.printMatrix(pM);

//        System.out.println(pM[1][1]);
//        System.out.println(Math.pow(Math.E, pM[1][1]));

        return pM;
    }

    /**
     *
     * @param profile1
     * @param profile2
     * @param i
     * @param j
     * @return
     */
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

    /**
     *
     * @param profile1
     * @param profile2
     * @param i
     * @param j
     * @param subMatrix
     * @return
     */
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

    /**
     *
     * @return
     */
    public double getIndelEmission(HashProfile profile, int pos, int state){

        double totalCount = getTotalIndelCount(profile, pos, state);

        double totalScore = getTotalIndelScore(profile, pos, state);

        double emission = totalScore / (totalCount);

        return emission;

    }



    public double getTotalIndelScore(HashProfile profile, int pos, int state) {
        double totalScore = 0;

        for (Character character : profile.getProfileArray().get(pos - 1).keySet()) {
            if (character!= '-') {
//                if (MatrixUtils.returnNucIndex(character) == -1){
//                    System.out.println("neg one");
//                }
                if (MatrixUtils.returnIndex(character,type ) == -1){
                    System.out.println(type);
                    System.out.println(character);

                }
                totalScore += emissions[state][MatrixUtils.returnIndex(character, type)];
            }
        }

        return totalScore;
    }

    public double getTotalIndelCount(HashProfile profile, int pos, int state) {

        double profileCount = 0;
        if (pos == 84 && pos > profile.getProfileArray().size()){
            System.out.println("here");
        }

        for (Character character : profile.getProfileArray().get(pos- 1).keySet()) {
            profileCount =+ profile.getProfileArray().get(pos - 1).get(character).getValue();
        }

        return profileCount;
    }




    /**
     *
     * @param profile1
     * @param profile2
     * @param i
     * @param j
     * @return
     */
    public double getEmission(HashProfile profile1, HashProfile profile2, int i, int j){
        double emission;
        if (i > profile1.getProfileArray().size()  || j > profile2.getProfileArray().size()){
            emission = 0;
        }
        else {

            double totalCount = getTotalCount(profile1, profile2, i, j);

            double totalScore = getTotalScore(profile1, profile2, i, j, subMatrix);

            emission = totalScore / (totalCount);
        }

        return emission;
    }

    /**
     *
     * @param profile1Matches
     * @param profile2Matches
     */
    public void addGapsToProfiles(List<Integer> profile1Matches, List<Integer> profile2Matches){
        List<Integer> gapPos = new ArrayList<Integer>();
        List<Integer> gapPos2 = new ArrayList<Integer>();





        for (int k = 0; k < profile1Matches.size(); k++){
            if (profile1Matches.get(k) == -1){
                gapPos.add(k);

            }
        }

        for (int k = 0; k < profile2Matches.size(); k++){
            if (k== 564){
//                System.out.println("here");
            }
            if (profile2Matches.get(k) == -1){
                gapPos2.add(k);
            }

        }


//        System.out.println("Profile 1 Size: " + profile1.getProfileArray().size());
//        System.out.println("Profile 2 Size: " + profile2.getProfileArray().size());
//        System.out.println("Profile 1 Matches Size: " + profile1Matches.size());
//        System.out.println("Profile 2 Matches Size: " + profile2Matches.size());
//        System.out.println("Gap Pos 1 Size: " + gapPos.size());
//        System.out.println("Gap Pos 2 Size: " + gapPos2.size());




        if (gapPos.size() > 0){
            profile1.addGaps(gapPos);

        }

        if (gapPos2.size() > 0) {
            profile2.addGaps(gapPos2);
        }
    }

    double getFValue(){
        return this.fValue;
    }

    double getBValue(){
        return this.bValue;
    }

    double log_add(double x, double y) {
        if(x == Double.NEGATIVE_INFINITY)
            return y;
        if(y == Double.NEGATIVE_INFINITY)
            return x;
        else if (x < y)
            return y + Math.log (1 + Math.exp(x-y));
        else
            return x + Math.log (1 + Math.exp(y-x));
    }



}
