package Alignment;

import Alignment.Utilities.MatrixUtils;
import SubstitutionModels.SubstitutionMatrix;

import java.util.Arrays;
import java.util.List;

import java.util.ArrayList;

/**
 * Class to perform parameter estimation on Hidden Markov Models using Baum Welch algorithm
 *
 */
public class BaumWelchMulti {

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
    //    int seqLength;
    String[] seqArray;
//    int[][] intArray;

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
    public BaumWelchMulti(String[] seqArray, String type){


        double [] start = {.5, .25, .25};

        double [][] transition = {
                {.8, .1, .1},
                {.5, .4, .1},
                {.5, .1, .4}
        };
        this.numAlpha = 21;

        double[][] emission = makeEmissionArray(start.length, numAlpha);

        SubstitutionMatrix subMatrix = new SubstitutionMatrix("Blosum62WithExtimatedX");



        BaumWelchMulti bw = new BaumWelchMulti(seqArray, start, transition, emission, subMatrix, type);

        this.start = bw.getStart();
        this.transition = bw.getTransition();
        this.emission = bw.getEmission();
        this.type = type;


    }


    public BaumWelchMulti(String[] seqArray, double[] start, double[][] transition, double[][] emission, SubstitutionMatrix subMatrix, String type) {

        //TODO: Work out which sequence length is needed

        this.numSeqs = seqArray.length;
        this.numStates = start.length;
        this.numAlpha = 4;
//        this.seqLength = seqArray[0].length();
        this.seqArray = seqArray;


        this.start = start;
        this.transition = transition;
        this.emission = emission;

        this.type = type;

//        intArray = new int[numSeqs][seqLength];





        boolean converged = false;





//        double [] start = {0.5, 0.25, 0.25};
//        double [][] transition = {
//                {0.70, 0.15, 0.15},
//                {0.70, 0.30, 0},
//                {0.70, 0, 0.30}};
//
////        double [][] emission = {
////                {0.25, 0.25, 0.25, 0.25},
////                {0.25, 0.35, 0.10, 0.30},
////                {0.25, 0.20, 0.05, 0.50}};
//
//        double [][] emission = {
//                {0.25, 0.25, 0.25, 0.25},
//                {0.25, 0.25, 0.25, 0.25},
//                {0.25, 0.25, 0.25, 0.25}};

//        String [] seqArray = {
//                "AGAAAGGTCTAGTGTTTGGTGATGTATCTATAGAGGGACGGGGGGGGGGGGGGGG",
//                "GGTCCTTTCAATATCAGTTGAATATGATGTGAGTGAGTTG",
//                "GGGGGGTGGGGCCTTGATAAGAAGGGCTGTCTTTTGGTAG",
//                "GTACCGGTATAGAAAAGACCGGATTCGAATTAATAATAAG",
//                "TATTACTTGTTCAGCGTTATAAGATTCAGGAGGAGGTGTG"};

//        String [] seqArray = {
//                "AAAAA",
//                "AAGGA",
//                "AACCA",
//                "AACGG"
//        };

//        // Parameters for the Moss (2008) example
//        String [] seqArray = {
//                "ATTA",
//                "TATT"
//        };
//
//        double [] start = {.85, .15};
//
//        double [][] transition = {
//                {.3, .7},
//                {.9, .1},
//        };
//
//        double [][] emission = {
//                {.4, .6 },
//                {.5, .5},
//        };

        // Parameters for the Moss (2008) example







        // Convert sequences

//        if (type.equals("nucleotide")){
//
//            for (int i=0; i < numSeqs; i++){
//                for (int j=0; j< seqArray[i].length(); j++) {
//                    switch (seqArray[i].charAt(j)) {
//                        case 'A':
//                            intArray[i][j] = 0;
//                            break;
//                        case 'T':
//                            intArray[i][j] = 1;
//                            break;
//                        case 'C':
//                            intArray[i][j] = 2;
//                            break;
//                        case 'G':
//                            intArray[i][j] = 3;
//                            break;
//                    }
//
//                }
//            }
//
//        }
//
//        else if (type.equals("protein")) {
//
//
//            for (int i = 0; i < numSeqs; i++) {
//                for (int j = 0; j < seqLength; j++) {
//                    switch (seqArray[i].charAt(j)) {
//                        case 'A':
//                            intArray[i][j] = 0;
//                            break;
//                        case 'R':
//                            intArray[i][j] = 1;
//                            break;
//                        case 'N':
//                            intArray[i][j] = 2;
//                            break;
//                        case 'D':
//                            intArray[i][j] = 3;
//                            break;
//                        case 'C':
//                            intArray[i][j] = 4;
//                            break;
//                        case 'Q':
//                            intArray[i][j] = 5;
//                            break;
//                        case 'E':
//                            intArray[i][j] = 6;
//                            break;
//                        case 'G':
//                            intArray[i][j] = 7;
//                            break;
//                        case 'H':
//                            intArray[i][j] = 8;
//                            break;
//                        case 'I':
//                            intArray[i][j] = 9;
//                            break;
//                        case 'L':
//                            intArray[i][j] = 10;
//                            break;
//                        case 'K':
//                            intArray[i][j] = 11;
//                            break;
//                        case 'M':
//                            intArray[i][j] = 12;
//                            break;
//                        case 'F':
//                            intArray[i][j] = 13;
//                            break;
//                        case 'P':
//                            intArray[i][j] = 14;
//                            break;
//                        case 'S':
//                            intArray[i][j] = 15;
//                            break;
//                        case 'T':
//                            intArray[i][j] = 16;
//                            break;
//                        case 'W':
//                            intArray[i][j] = 17;
//                            break;
//                        case 'Y':
//                            intArray[i][j] = 18;
//                            break;
//                        case 'V':
//                            intArray[i][j] = 19;
//                            break;
//                    }
//                }
//            }
//        }




        // Compute recurrance
        while (!converged){
            double [] startR = new double [3];
            double [][] transitionR = new double [numStates][numStates];
            double [][] emissionR = new double [numStates][numAlpha];

            for (int i = 0; i < numStates; i++){
                startR[i] = ZERO_CONSTANT;
            }

            for (int i = 0; i< numStates; i++){
                for (int j=0; j< numStates; j++) {
                    transitionR[i][j] = ZERO_CONSTANT;
                }
            }

            for (int i = 0; i < numStates; i++){
                for (int j=0; j< numAlpha; j++){
                    emissionR[i][j] = ZERO_CONSTANT;
                }
            }

            for (int n = 0; n < numSeqs; n++){
                String seq = seqArray[n];

                double [][] f = new double[numStates][seq.length()];
                double [][] b = new double [numStates][seq.length()];

                for (int i = 0; i < numStates; i++){
                    for (int j = 0; j < seq.length(); j++){
                        f[i][j] = ZERO_CONSTANT;
                        b[i][j] = ZERO_CONSTANT;
                    }
                }

                // Calculate the forward and backwards values
//                fValue = forward (n, f, start, transition, emission);
                PairHMMUnderflow pairHMM = new PairHMMUnderflow(seqArray, start, transition, emission, subMatrix, false, type);
                pairHMM.forwardAlgorithm();
                double fValue2 = pairHMM.getFValue();
//                fValue2 = Math.pow(Math.E, fValue2);
//                pairHMM.backwardAlgorithm();
                bValue = backward(n, b, start, transition, emission);


                MatrixUtils.printMatrix(f);


                System.out.println("Multi fValue " + fValue + " and fValue2 " + fValue2);
                System.out.println("Done");


                //Update start matrix
                for (int k = 0; k < numStates; k++) {
                    startR[k] += (1 / fValue) * f[k][0] * b[k][0];


                    //Update transition matrix
                    for (int l = 0; l < numStates; l++) {
                        for (int i = 0; i < seq.length() - 2; i++) {
//                            System.out.println("Emission: " +  emission[l][MatrixUtils.returnIndex(seq.charAt(i+1))]);
                            transitionR[k][l] += (1/fValue) * f[k][i] * transition[k][l] * emission[l][MatrixUtils.returnIndex(seq.charAt(i+1), type)] * b[l][i+1];
//                            System.out.println(transitionR[k][l]);

                        }
                    }

                    //Update emission matrix
                    for (int i=0; i < seq.length(); i++){
                        emissionR[k][MatrixUtils.returnIndex(seq.charAt(i), type)] += (1/fValue) * f[k][i] * b[k][i];
//                        System.out.println(emissionR[k][MatrixUtils.returnIndex(seq.charAt(1))]);
                    }


                }

//                for (double value: startR){
//                    System.out.println(value);
//
//                }
//                MatrixUtils.printMatrix(emissionR);
//                MatrixUtils.printMatrix(transitionR);


            }

            //Update main arrays and find max difference

            double sum = 0;
            double startDiff = 0;

            for (int i = 0; i< numStates; i++){
                sum += startR[i];
            }

            for (int i=0; i <numStates; i++){
                if (Math.abs((startR[i] / sum) - start[i]) > startDiff){
                    startDiff = Math.abs((startR[i] / sum) - start[i]);
                }
                start[i] = startR[i] / sum;
            }

            // Update transition
            transitionDiff = 0.0;
            emissionDiff = 0.0;



            List<Object> updatedTransitionValues = getUpdatedValues(transition, transitionR, numStates, transitionDiff);
            List<Object> updatedEmissionValues = getUpdatedValues(emission, emissionR, numAlpha, emissionDiff);



            transition = (double[][]) updatedTransitionValues.get(0);
            transitionDiff = (Double) updatedTransitionValues.get(1);

            emission = (double[][]) updatedEmissionValues.get(0);
            emissionDiff = (Double) updatedEmissionValues.get(1);



            // Check for convergence
            if (startDiff + transitionDiff + emissionDiff < 0.000000000001) {
                converged = true;
            }

            this.start = start;
            this.transition = transition;
            this.emission = emission;

        }



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
        HashProfile profile1 = new HashProfile(seqArray[0]);
        HashProfile profile2 = new HashProfile(seqArray[1]);
        int seqLength = seq.length();
        fValue = 0;

        f[0][0] = 1;

//        for (int j = 0; j < profile1.getLength(); j++){
//            for (int i = 0; i < profile2.getLength(); i++){
//                if ((i == 1) && ( j==1)){
//                    System.out.println(i + " " + j + " Start - Match");
////                    double emissionM = PairHMM.getEmission(profile1, profile2, i, j);
////                    f[i][j] = start[0] * emission[0][MatrixUtils.returnIndex(profile2.getProfileArray().get() type)] * f[i][j-1];
//
//
//                }
//
//                if (i >= 1 && j >=1){
//                    System.out.println(i + " " + j + " Match, Insert, Delete");
//
//                }
//
//                if (i == 0 && j == 1){
//                    System.out.println(i + " " + j + " Start - Delete");
//                    f[i][j] = start[2] * emission[2][MatrixUtils.returnIndex(profile2.charAt(0), type)] * f[i][j-1];
//
//
//                }
//
//                if (j >= 1){
//                    System.out.println(i + " " + j + " Delete");
//
//                }
//
//                if (i==1 && j==0){
//                    System.out.println(i + " " + j + " Start - Insert");
//                    f[i][j] = start[1] * emission[1][MatrixUtils.returnIndex(profile1.charAt(0), type)] * f[i-1][j];
//
//
//                }
//
//                if ( i>= 1){
//                    System.out.println(i + " " + j + " Insert");
//
//                }
//            }
//        }


        for(int i=0; i < numStates; i++){


            f[i][0] = start[i] * emission[i][MatrixUtils.returnIndex(seq.charAt(0), type)];
        }

        for (int i=1; i< seqLength; i++){
            for (int j=0; j< numStates; j++){
                for(int k=0; k< numStates; k++){
                    // Emission at new cell multipled by the cost of transitioning from previous cell
//                    System.out.println(emission[j][intArray[n][i]]);
//                    System.out.println(f[k][i-1]);
//                    System.out.println(transition[k][j]);
//                    System.out.println(emission[j][intArray[n][i]] * f[k][i-1] * transition[k][j]);
//                    System.out.println("SEQ HERE:" + seq + " " + i);
                    f[j][i] += emission[j][MatrixUtils.returnIndex(seq.charAt(i), type)] * f[k][i-1] * transition[k][j];
                }
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
        int seqLength = seq.length();

        bValue = 0;
        for(int i=0; i< numStates; i++){
            b[i][seqLength-1] = 1;
        }

        for (int i = seqLength - 2; i>= 0; i--){
            for(int j=0; j<numStates; j++){
                for(int k=0; k<numStates; k++){
                    b[j][i] += emission[k][MatrixUtils.returnIndex(seq.charAt(i+1), type)] * b[k][i+1] * transition[j][k];
                }
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


}
