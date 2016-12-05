package Alignment;

import Alignment.Utilities.MatrixUtils;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import java.util.ArrayList;

/**
 * Class to perform parameter estimation on Hidden Markov Models using Baum Welch algorithm
 *
 */
public class BaumWelch {

//    int numSeqs = 4;
//    int numStates = 3;
//    int numAlpha = 4;
//    int seqLength = 5;


//    // Parameters for the Moss (2008) example
//    int numSeqs = 2;
//    int numStates = 2;
//    int numAlpha = 4;
//    int seqLength = 4;

    // Parameters for the dog example
//    int numSeqs = 9;
//    int numStates = 2;
//    int numAlpha = 4;
//    int seqLength = 2;

    int numSeqs;
    int numStates;
    int numAlpha;
    int seqLength;
    String[] seqArray;
    int[][] intArray;

    double[] start;
    double[][] emission;
    double[][] transition;

    double transitionDiff;
    double emissionDiff;

    private String type;



    double fValue = 0;
    double bValue = 0;
    double ZERO_CONSTANT = 0.00000000000000000000000000001;

    // Default constructor when start, transition, and emission stating values are not provided
    public BaumWelch(String[] seqArray, String type){


        double [] start = {.5, .25, .25};

        double [][] transition = {
                {.6, .2, .2},
                {.5, .4, .1},
                {.5, .1, .4}
        };
        this.numAlpha = 20;

        double[][] emission = makeEmissionArray(start.length, numAlpha);



        BaumWelch bw = new BaumWelch(seqArray, start, transition, emission, type);

        this.start = bw.getStart();
        this.transition = bw.getTransition();
        this.emission = bw.getEmission();
        this.type = type;


    }


    public BaumWelch(String[] seqArray, double[] start, double[][] transition, double[][] emission, String type) {

        //TODO: Work out which sequence length is needed

        this.numSeqs = seqArray.length;
        this.numStates = start.length;
        this.numAlpha = 3;
//        this.seqLength = seqArray[0].length();
        this.seqArray = seqArray;


        this.start = start;
        this.transition = transition;
        this.emission = emission;

        this.type = type;

        intArray = new int[numSeqs][seqLength];






        boolean converged = false;






        // Compute recurrance
        while (!converged) {
            MatrixUtils.printMatrix(transition);
            MatrixUtils.printMatrix(emission);
//            double[] startR = new double[3];
//            double[][] transitionR = new double[numStates][numStates];
//            double[][] emissionR = new double[numStates][numAlpha];

//            for (int i = 0; i < numStates; i++) {
//                startR[i] = ZERO_CONSTANT;
//            }
//
//            for (int i = 0; i < numStates; i++) {
//                for (int j = 0; j < numStates; j++) {
//                    transitionR[i][j] = ZERO_CONSTANT;
//                }
//            }
//
//            for (int i = 0; i < numStates; i++) {
//                for (int j = 0; j < numAlpha; j++) {
//                    emissionR[i][j] = ZERO_CONSTANT;
//                }
//            }

            for (int n = 0; n < numSeqs; n++) {
                String seq = seqArray[n];
                int seqLength = seq.length();
                double[][] f = new double[numStates][seqLength + 1];
                double[][] b = new double[numStates][seqLength + 1];

                for (int i = 0; i < numStates; i++) {
                    for (int j = 0; j < seqLength + 1; j++) {
                        f[i][j] = ZERO_CONSTANT;
                        b[i][j] = ZERO_CONSTANT;
                    }
                }

                // Calculate the forward and backwards values
                fValue = forward(n, f, start, transition, emission);
                bValue = backward(n, b, start, transition, emission);
                double[][] probMatrix = calculateProbMatrix(f, b);
                double[][][] transitionProbs = calculateTransitionProbs(f, b, seq, transition, emission);
                double[][] updatedTransition = calculateTransitionMatrix(transitionProbs, probMatrix, f, b, transition, emission);
                System.out.println("here?");

                double[][] updatedEmission = calculateEmissionMatrix(probMatrix, seq);
                System.out.println("here");
                double transitionDiff = getDiff(transition, updatedTransition);
                double emissionDiff = getDiff(emission, updatedEmission);

//                List<Object> updatedTransitionValues = getUpdatedValues(transition, updatedTransition, numStates, transitionDiff);
//                List<Object> updatedEmissionValues = getUpdatedValues(emission, updatedEmission, numAlpha, emissionDiff);

                double startDiff = 0;


                System.out.println(transitionDiff);

                // Check for convergence
                if (transitionDiff < 0.00000001 && emissionDiff < 0.00000001) {
                    converged = true;

                }
//                transitionDiff = (Double) updatedTransitionValues.get(1);
//
//                emissionDiff = (Double) updatedEmissionValues.get(1);
//
//                System.out.println(transitionDiff);


//                // Check for convergence
//                if (startDiff + transitionDiff + emissionDiff < 0.000000000001) {
//                    converged = true;
//                }

                start = start;
                transition = updatedTransition;
                emission = updatedEmission;

            }
        }

////                System.out.println("Original fValue " + fValue + " and bValue " + bValue);
//
//
//                //Update start matrix
//                for (int k = 0; k < numStates; k++) {
//                    startR[k] += (1 / fValue) * f[k][0] * b[k][0];
//
//
//                    //Update transition matrix
//                    for (int l = 0; l < numStates; l++) {
//                        for (int i = 0; i < seqLength - 2; i++) {
////                            System.out.println("Emission: " + emission[l][intArray[n][i+1]]);
//                            transitionR[k][l] += (1/fValue) * f[k][i] * transition[k][l] * emission[l][intArray[n][i+1]] * b[l][i+1];
////                            System.out.println(transitionR[k][l]);
//
//                        }
//                    }
//
//                    //Update emission matrix
//                    for (int i=0; i < seqLength; i++){
//                        emissionR[k][intArray[n][i]] += (1/fValue) * f[k][i] * b[k][i];
////                        System.out.println(emissionR[k][intArray[n][i]]);
//                    }
//
//
//                }
////                for (double value: startR){
////            System.out.println(value);
////
////        }
////        MatrixUtils.printMatrix(emissionR);
////        MatrixUtils.printMatrix(transitionR);
////
//            }
//
//            //Update main arrays and find max difference
//
//            double sum = 0;
//            double startDiff = 0;
//
//            for (int i = 0; i< numStates; i++){
//                sum += startR[i];
//            }
//
//            for (int i=0; i <numStates; i++){
//                if (Math.abs((startR[i] / sum) - start[i]) > startDiff){
//                    startDiff = Math.abs((startR[i] / sum) - start[i]);
//                }
//                start[i] = startR[i] / sum;
//            }
//
//            // Update transition
//            transitionDiff = 0.0;
//            emissionDiff = 0.0;
//
//
//
//            List<Object> updatedTransitionValues = getUpdatedValues(transition, transitionR, numStates, transitionDiff);
//            List<Object> updatedEmissionValues = getUpdatedValues(emission, emissionR, numAlpha, emissionDiff);
//
//
//
//            transition = (double[][]) updatedTransitionValues.get(0);
//            transitionDiff = (Double) updatedTransitionValues.get(1);
//
//            emission = (double[][]) updatedEmissionValues.get(0);
//            emissionDiff = (Double) updatedEmissionValues.get(1);
//
//
//
//            // Check for convergence
//            if (startDiff + transitionDiff + emissionDiff < 0.000000000001) {
//                converged = true;
//            }
//
//            this.start = start;
//            this.transition = transition;
//            this.emission = emission;
//
//        }



    }

    /**
     * Calculate the forward value
     * @param n Number of iterations
     * @param f Array of forward values
     * @param start Array of starting state probabilities
     * @param transition Array of transition probabilities
     * @param emission Array of emission probabilities
     * @return Double representing the updated forward probability
     */

    double forward (int n, double [][] f, double [] start, double[][] transition, double [][] emission){
        String seq = seqArray[n];
        int seqLength = seq.length();
        double startSum = 0;
        fValue = 0;

        for(int i=0; i< numStates; i++){
            f[i][0] = start[i];
            startSum += f[i][0];

        }

        // Normalise the starting cells
        for (int i=0; i < numStates; i++){
            f[i][0] = f[i][0] / startSum;
        }

//        // The starting cells are found by multiplying the starting probability by the emission and transition probabilities
//        for(int i=0; i < numStates; i++) {
//            for (int j = 0; j < numStates; j++) {
//
//
//                f[i][0] += start[j] * emission[i][intArray[n][0]] * transition[j][i];
//            }
//            startSum += f[i][0];
//        }



        for (int i=0; i< seqLength; i++){
            double cellSum = 0;
            for (int j=0; j< numStates; j++){
                for(int k=0; k< numStates; k++){
                    // Emission at new cell multipled by the cost of transitioning from previous cell
//                    System.out.println(emission[j][intArray[n][i]]);
//                    System.out.println(f[k][i]);
//                    System.out.println(transition[k][j]);
//                    System.out.println(emission[j][intArray[n][i]] * f[k][i] * transition[k][j]);
                    f[j][i + 1] += emission[j][MatrixUtils.returnIndex(seq.charAt(i), type)] * f[k][i] * transition[k][j];
                }
                cellSum += f[j][i+1];
            }

            // Normalise the cells
            for (int l=0; l < numStates; l++){
                f[l][i + 1] = f[l][i + 1] / cellSum;
            }

        }


        for (int i=0; i <numStates; i++){
            fValue += f[i][seqLength-1];
        }

        return fValue;

    }

    /**
     * Calculate the backwards value
     * @param n Number of iterations
     * @param b Array of backwards values
     * @param start Array of starting state probabilities
     * @param transition Array of transition probabilities
     * @param emission Array of emission probabilities
     * @return Double representing the updated backward probability
     */

    double backward (int n, double [][] b, double [] start, double[][] transition, double [][] emission){
        String seq = seqArray[n];
        seqLength = seq.length();

        bValue = 0;
        for(int i=0; i< numStates; i++){
            b[i][seqLength] = 1;
        }

        for (int i = seqLength - 1; i>= 0; i--){
            double cellSum = 0;
            for(int j=0; j<numStates; j++){
                for(int k=0; k<numStates; k++){

                    b[j][i] += emission[k][MatrixUtils.returnIndex(seq.charAt(i), type)] * b[k][i + 1] * transition[j][k];

//                    b[j][i] += emission[k][intArray[n][i]] * b[k][i + 1] * transition[j][k];
                }
                cellSum += b[j][i];
            }

            // Normalise the cells
            for (int l=0; l < numStates; l++){
                b[l][i] = b[l][i] / cellSum;
            }
        }

        for (int i = 0; i < numStates; i++){
            bValue += b[i][0] * emission[i][MatrixUtils.returnIndex(seq.charAt(0), type)] * start[i];
        }

        return bValue;
    }

    /**
     * Updates the current matrix and returns the updated matrix and the difference between the two matrices
     * @param currentMatrix The current matrix to update
     * @param updatedMatrix The newer matrix to update
     * @param length
     * @param diff
     * @return
     */
    public List<Object> getUpdatedValues(double[][] currentMatrix, double[][] updatedMatrix, int length, Double diff){
        for (int i = 0; i < numStates; i++){
            double sum = 0;
            for (int j = 0; j< length; j++){
                sum += updatedMatrix[i][j];
            }
            for (int j=0; j< length; j++){
                if (Math.abs((updatedMatrix[i][j] / sum) - currentMatrix[i][j]) > diff) {
                    diff = Math.abs((updatedMatrix[i][j] / sum) - currentMatrix[i][j]);
                }
                currentMatrix[i][j] = updatedMatrix[i][j] / sum;
            }
        }


        List<Object> updatedValues = new ArrayList<Object>();
        updatedValues.add(currentMatrix);
        updatedValues.add(diff);
        return updatedValues;
    }

    /**
     * Construct a generic emission array with uniform emission probabilities for each character in each state
     * @param states The number of states in the model
     * @param characters The number of characters in the model
     * @return Double array of emission probabilities
     */
    public double[][] makeEmissionArray(int states, int characters){

        double [][] emissionArray = new double[states][characters];

        for (double[] array : emissionArray){
            Arrays.fill(array, 1.0 / characters);

        }

        return emissionArray;

    }

    public double[][] calculateProbMatrix(double[][] f, double[][]b){
        double[][] p = new double[numStates][seqLength + 1];
        for (int j = 0; j < seqLength + 1; j++){
            double cellSum = 0;
            for (int i=0; i< numStates; i++){
                p[i][j] = f[i][j] * b[i][j];
                cellSum += p[i][j];

            }

            // Normalise the cells
            for (int k=0; k < numStates; k++){
                p[k][j] = p[k][j] / cellSum;
            }

        }

        return p;

    }

    public double[][][] calculateTransitionProbs(double[][] f, double[][] b, String seq, double[][] transition, double[][] emission){
        double[][][] updatedTransitions = new double[numStates][numStates][seqLength];

        for (int i=0; i < numStates; i++){
            for (int j=0; j < numStates; j++){
                for (int k = 0; k < seqLength; k++) {

                    updatedTransitions[i][j][k] = f[i][k] * b[j][k+1] * transition[i][j] * emission[j][MatrixUtils.returnIndex(seq.charAt(k), type)];

//                    updatedTransitions[i][j][k] = f[i][k] * b[j][k+1] * transition[i][j] * emission[j][intArray[0][k]];
                }
            }

        }

        return updatedTransitions;

    }

    public double[][] calculateTransitionMatrix(double[][][] transitionProbs, double[][] probMatrix, double[][] f, double[][] b, double[][]
                                               transition, double[][] emission){
        double[][] updatedTransition = new double[numStates][numStates];


        // Update to the sum of all possible transitions from state - state / probability of being in this state
        for (int i = 0; i < numStates; i++){
            for (int j = 0; j < numStates; j++){
                updatedTransition[i][j] = sumRow(transitionProbs[i][j], 0) / sumRow(probMatrix[i], 0);

            }
        }

            double[] sumRows = new double[numStates];
            for (int k = 0; k < numStates; k++){
                 sumRows[k] = sumRow(updatedTransition[k], 0);

            }
        for (int i = 0; i < numStates; i++) {

            double sumRow = sumRow(updatedTransition[i], 0);
            for (int k = 0; k < numStates; k++) {
//                updatedTransition[i][k] = updatedTransition[i][k] / sumRows[k];
                updatedTransition[i][k] = updatedTransition[i][k] / sumRow;

            }
        }
        return updatedTransition;
    }

    public double[][] calculateEmissionMatrix(double[][] probMatrix, String seq){
        double[][] updatedEmission = new double[numStates][numStates];

        System.out.println("and here?");

//        MatrixUtils.printMatrix(probMatrix);
        for (int i = 0; i < numStates; i++){

//                    Character[] alphabet = MatrixUtils.getAlphabet(type)
            HashSet<Character> alphabet = getAlphabet(seq);
                    for (Character letter : alphabet) {
                        double sum = 0;



                            for (int l = 0; l < seqLength; l++) {
                                if (letter.equals(seq.charAt(l))) {
//                                    System.out.println("MATCH");
//                                    System.out.println(letter);
//                                    System.out.println(probMatrix[i][l + 1]);
                                    sum += probMatrix[i][l + 1];
//                                    System.out.println(sum);
                                    updatedEmission[i][MatrixUtils.returnIndex(letter, type)] += probMatrix[i][l+1];
                                }


                        }
                        updatedEmission[i][MatrixUtils.returnIndex(letter, type)] = updatedEmission[i][MatrixUtils.returnIndex(letter, type)] / sumRow(probMatrix[i], 1);

                    }
        }
        return updatedEmission;
    }

    public double getDiff(double[][] original, double[][] updated){

        double sum = 0;
        for (int j = 0; j < original.length; j++) {
            for (int i = 0; i < original[0].length; i++) {
                sum += Math.pow(original[i][j] - updated[i][j], 2);
            }
        }
        double diff = Math.sqrt(sum);
        return diff;
    }

    /**
     * @return Array of doubles representing start probabilities
     */
    public double[] getStart(){
        return this.start;
    }

    /**
     * @return Array of doubles representing transition probabilities
     */
    public double[][] getTransition(){
        return this.transition;
    }

    /**
     * @return Array of doubles representing emission probabilities
     */
    public double[][] getEmission(){
        return this.emission;
    }

    public HashSet<Character> getAlphabet(String seq){
        HashSet<Character> alphabet = new HashSet<Character>();
        while (alphabet.size() < numAlpha){
        for (int i = 0; i < seqLength; i++) {
            if (!alphabet.contains(seq.charAt(i))) {
                alphabet.add(seq.charAt(i));
            }
        }

        }
        return alphabet;
    }

    public double sumRow(double[] row, int startPos){
        double total = 0;
        for (int val = startPos; val < row.length; val++){
            total += row[val];
        }
        return total;
    }



}
