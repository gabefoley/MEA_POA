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
        NumberFormat formatter = new DecimalFormat();
        formatter = new DecimalFormat("0.#####E0");

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", formatter.format(matrix[i][j])) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }


    public static void printMatrix(int[][] matrix) {
        NumberFormat formatter = new DecimalFormat();
        formatter = new DecimalFormat("0.#####E0");

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

    //TODO: Move this helper method to the matrix helper methods
    public static double[] addArrays(double[] firstArray, double[] secondArray) {

        double[] shorterArray = (firstArray.length < secondArray.length ? firstArray : secondArray);
        double[] longerArray = (firstArray.length > secondArray.length ? firstArray : secondArray);
        double[] finalArray = new double[longerArray.length];


        for (int i = 0; i < shorterArray.length; i++) {
            finalArray[i] = (firstArray[i] + secondArray[i]);
        }
        for (int i = shorterArray.length; i < longerArray.length; i++) {
            finalArray[i] = longerArray[i];
        }
        return finalArray;
    }


}


