package Alignment.Utilities;

import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * Utilities for printing and adding matrices
 */
public class MatrixUtils {



    public static void printMatrix(double[][] matrix) {
        for (double[] aMatrix : matrix) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", aMatrix[j]) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }

    public static void printMatrix2(double[][] matrix) {
        NumberFormat formatter = new DecimalFormat("0.#####E0");
        for (double[] aMatrix : matrix) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", formatter.format(aMatrix[j])) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }


    public static void printMatrix(int[][] matrix) {
        NumberFormat formatter = new DecimalFormat("0.#####E0");
        for (int[] aMatrix : matrix) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", formatter.format(aMatrix[j])) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }

    public static void printMatrix(String[][] matrix) {


        for (String[] aMatrix : matrix) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", aMatrix[j]) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }

    /**
     * Method to add two arrays together by adding the contents of each respective index in each array together.
     * For example
     * [1,4,9,6] + [2,2,3] = [3,7,12,6]
     *
     * @param firstArray  first array to add
     * @param secondArray second array to add
     * @return int array with the add
     */

    public static double[] addArrays(double[] firstArray, double[] secondArray) {

        double[] shorterArray = (firstArray.length < secondArray.length ? firstArray : secondArray);
        double[] longerArray = (firstArray.length > secondArray.length ? firstArray : secondArray);
        double[] finalArray = new double[longerArray.length];


        for (int i = 0; i < shorterArray.length; i++) {
            finalArray[i] = (firstArray[i] + secondArray[i]);
        }
        System.arraycopy(longerArray, shorterArray.length, finalArray, shorterArray.length, longerArray.length - shorterArray.length);
        return finalArray;
    }
    
    public static int returnIndex(Character character){
        
        switch (character) {
            case 'A':
               return 0;
            case 'R':
               return 1;
            case 'N':
               return 2;
            case 'D':
               return 3;
            case 'C':
               return 4;
            case 'Q':
               return 5;
            case 'E':
               return 6;
            case 'G':
               return 7;
            case 'H':
               return 8;
            case 'I':
               return 9;
            case 'L':
               return 10;
            case 'K':
               return 11;
            case 'M':
               return 12;
            case 'F':
               return 13;
            case 'P':
               return 14;
            case 'S':
               return 15;
            case 'T':
               return 16;
            case 'W':
               return 17;
            case 'Y':
               return 18;
            case 'V':
               return 19;
            default:
                return -1;
        }
    }


}


