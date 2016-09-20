package Alignment;

import com.sun.corba.se.spi.extension.ZeroPortPolicy;

/**
 * Created by gabe on 13/09/2016.
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
    int[][] intArray;

    double fValue = 0;
    double bValue = 0;
    double ZERO_CONSTANT = 0.00000000000000000000000000001;


    public BaumWelch(String[] seqArray, double[] start, double[][] transition, double[][] emission) {

        //TODO: Work out which sequence length is needed

        this.numSeqs = seqArray.length;
        this.numStates = start.length;
        this.numAlpha = emission.length;
        this.seqLength = seqArray[0].length();

        int[][] intArray = new int[numSeqs][seqLength];




        boolean converged = false;
        int iteration = 0;





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

        for (int i=0; i < numSeqs; i++){
            for (int j=0; j< seqLength; j++){
                switch(seqArray[i].charAt(j)) {
                    case 'A':
                        intArray[i][j] = 0;
                        break;
                    case 'T':
                        intArray[i][j] = 1;
                        break;
                    case 'C':
                        intArray[i][j] = 2;
                        break;
                    case 'G':
                        intArray[i][j] = 3;
                        break;
                }
            }
        }

        // Compute recurrance
        while (converged == false){
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
                double [][] f = new double[numStates][seqLength];
                double [][] b = new double [numStates][seqLength];

                for (int i = 0; i < numStates; i++){
                    for (int j = 0; j < seqLength; j++){
                        f[i][j] = ZERO_CONSTANT;
                        b[i][j] = ZERO_CONSTANT;
                    }
                }

                fValue = forward (n, f, start, transition, emission, fValue);
                bValue = backward(n, b, start, transition, emission, bValue);

                // Update recurrance arrays
                //Update start matrix
                for (int k = 0; k < numStates; k++) {
                    startR[k] += (1 / fValue) * f[k][0] * b[k][0];


                    //Update transition matrix
                    for (int l = 0; l < numStates; l++) {
                        for (int i = 0; i < seqLength - 2; i++) {
                            transitionR[k][l] += (1/fValue) * f[k][i] * transition[k][l] * emission[l][intArray[n][i+1]] * b[l][i+1];
                        }
                    }

                    //Update emission matrix
                    for (int i=0; i < seqLength; i++){
                        emissionR[k][intArray[n][i]] += (1/fValue) * f[k][i] * b[k][i];
                    }


                }

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
            double transitionDiff = 0;
            for (int i=0; i< numStates; i++){
                sum = 0;
                for (int j=0; j< numStates; j++){
                    sum += transitionR[i][j];
                }
                for (int j = 0; j < numStates; j++){
                    if (Math.abs((transitionR[i][j] / sum) - transition[i][j]) > transitionDiff){
                        transitionDiff = Math.abs((transitionR[i][j] / sum) - transition[i][j]);
                    }
                    transition[i][j] = transitionR[i][j] / sum;
                }
            }

            // Update emission
            double emissionDiff = 0;
            for (int i = 0; i < numStates; i++){
                sum = 0;
                for (int j = 0; j< numAlpha; j++){
                    sum += emissionR[i][j];
                }
                for (int j=0; j< numAlpha; j++){
                    if (Math.abs((emissionR[i][j] / sum) - emission[i][j]) > emissionDiff) {
                        emissionDiff = Math.abs((emissionR[i][j] / sum) - emission[i][j]);
                    }
                    emission[i][j] = emissionR[i][j] / sum;
                }
            }

            // Check for convergence
            if (startDiff + transitionDiff + emissionDiff < 0.000000000001) {
                converged = true;
            }
            iteration++;

        }

        System.out.println("Done");



    }

    double forward (int n, double [][] f, double [] start, double[][] transition, double [][] emission, double fValue){
        fValue = 0;
        for(int i=0; i < numStates; i++){
            f[i][0] = start[i] * emission[i][intArray[n][0]];
        }

        for (int i=1; i< seqLength; i++){
            for (int j=0; j< numStates; j++){
                for(int k=0; k< numStates; k++){
                    f[j][i] += emission[j][intArray[n][i]] * f[k][i-1] * transition[k][j];
                }
            }
        }

        for (int i=0; i <numStates; i++){
            fValue += f[i][seqLength-1];
        }

        return fValue;

    }

    double backward (int n, double [][] b, double [] start, double[][] transition, double [][] emission, double bValue){

        bValue = 0;
        for(int i=0; i< numStates; i++){
            b[i][seqLength-1] = 1;
        }

        for (int i = seqLength - 2; i>= 0; i--){
            for(int j=0; j<numStates; j++){
                for(int k=0; k<numStates; k++){
                    b[j][i] += emission[k][intArray[n][i+1]] * b[k][i+1] * transition[j][k];
                }
            }
        }

        for (int i = 0; i < numStates; i++){
            bValue += b[i][0] * emission[i][intArray[n][0]] * start[i];
        }

        return bValue;
    }


}
