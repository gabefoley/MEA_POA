//package SubstitutionModels;
//
///**
// * Created by gabe on 31/07/2016.
// */
//public class ExampleModel {
//
//    // A,G,C,T
//    private static final double[][] matrix = new double[][]{
//            {0.50, 0.05, 0.15, 0.30},
//            {0.05, 0.50, 0.30, 0.15},
//            {0.15, 0.30, 0.50, 0.05},
//            {0.30, 0.15, 0.05, 0.50}};
//
//    private static int getIndex(char a) {
//
//        // check for upper and lowercase characters
//        switch ((String.valueOf(a)).toUpperCase().charAt(0)) {
//
//            case 'A': return 0;
//            case 'G': return 1;
//            case 'C': return 2;
//            case 'T': return 3;
//            default: return -1;
//        }
//    }
//
//    public static double getDistance(char a1, char a2) {
//        return matrix[getIndex(a1)][getIndex(a2)];
//    }
//
//}
