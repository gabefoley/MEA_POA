package Alignment;


import SubstitutionModels.Blosum62;
import SubstitutionModels.Blosum62Probs;

import java.util.*;

import static java.lang.Math.max;


public class Alignment {

    /**
     * Class for aligning a sequence to a PO Graph
     */

    private int openGapPenalty;
    private int extendGapPenalty;
    private String seq1;
    private String seq2;
    private HashProfile profile1;
    private HashProfile profile2;
    private List<Integer> stringIndexes;
    private List<Integer> nodeIndexes;
    private double[][] matrix;
    private boolean MEA;
    private int profileMatrixHeight;
    private HashProfile updatedProfile;

    public Alignment(String seq1, String seq2, int openGapPenalty, int extendGapPenalty, double[][] matrix, boolean MEA) {

        this(new HashProfile(seq1), new HashProfile(seq2), openGapPenalty, extendGapPenalty, matrix, MEA);
    }

    public Alignment(HashProfile profile1, HashProfile profile2, int openGapPenalty, int extendGapPenalty, double[][] matrix, boolean MEA) {
        this.profile1 = profile1;
        this.profile2 = profile2;
        this.openGapPenalty = openGapPenalty;
        this.extendGapPenalty = extendGapPenalty;
        this.stringIndexes = new ArrayList<Integer>();
        this.nodeIndexes = new ArrayList<Integer>();
        this.matrix = matrix;
        this.MEA = MEA;
//        this.profileMatrixHeight = profile1.getProfileArray()[0].length;
//        List<List<Integer>> matches = this.alignProfileToProfile();
        this.updatedProfile = this.alignProfileToProfile();

    }

    public HashProfile getUpdatedProfile(){
        return this.updatedProfile;
    }
    /**
     * Return the scores associated with matching a single character in a PO Graph to an
     * array of characters.
     * Used to get the scores for matching each row in the score matrix to each column
     * @param base the character in the PO Graph to be matched against
     * @param seqvec characters in the sequence to be matched against
     * @return an int array containing the match scores
     */
    public double[] getMatchScore(char base, String seqvec) {
        double[] matches = new double[seqvec.length()];
        for (int i = 0; i < seqvec.length(); i++) {
            double matchScore = Blosum62.getDistance(seqvec.charAt(i), base);
            matches[i] = matchScore;
        }

        return matches;
    }

    public double[] getMEAMatchScore(int i, String seqVec){
        double[] matches = new double[seqVec.length()];
        for (int j = 0; j < seqVec.length(); j++) {
            double matchScore = matrix[i+1][j+1];
            matches[j] = matchScore;

        }

        return matches;
    }

    public double[] getMEAMatchScore(int i, HashProfile profile){
        double[] matches = new double[profile.getProfileArray().size()];
        for (int j = 0; j < profile.getProfileArray().size(); j++) {
            double matchScore = matrix[i+1][j+1];
            matches[j] = matchScore;

        }

        return matches;
    }

    public double[] getProfileMatchScore(HashProfile profile1, HashProfile profile2, int index ){
        //TODO: Create a get length from profile method

        int profile2Length = profile2.getProfileArray().size();

        double[] matches = new double[profile2Length];


        //TODO: Normalise match score at individual index
        for (int i = 0; i < profile2Length; i++) {
            double totalScore = 0;
            double totalCount = 0;
            double profile1Count = 0;
            double profile2Count = 0;

            for (Character residue : profile1.getProfileArray().get(index).keySet()){
                profile1Count += profile1.getProfileArray().get(index).get(residue).getValue();
            }

            for (Character residue : profile2.getProfileArray().get(i).keySet()){
                profile2Count += profile2.getProfileArray().get(i).get(residue).getValue();
            }

            totalCount = profile1Count * profile2Count;

            for (Character name : profile1.getProfileArray().get(index).keySet()) {
                if (name != '-') {
                    Character profile1Name = name;
                    int profile1Value = profile1.getProfileArray().get(index).get(name).getValue();

                    for (Character name2 : profile2.getProfileArray().get(i).keySet()) {
                        if (name2 != '-') {
                            Character profile2Name = name2;
                            int profile2Value = profile2.getProfileArray().get(i).get(name2).getValue();
//                            System.out.println(profile1Value);
//                            System.out.println(profile2Value);

                            if (profile1Name == -1 || profile2Name == -1){
                                System.out.println("UMMM");
                            }

                            double matchScore = Blosum62.getDistance(profile1Name, profile2Name);
                            totalScore += profile1Value * profile2Value * matchScore;
                            matches[i] = totalScore / totalCount;

//                            System.out.println("Postion: " + i);
//                            System.out.println("Profile 1 residue: " + profile1Name + " and count: " + profile1Value);
//                            System.out.println("Profile 2 residue: " + profile2Name + " and count: " + profile2Value);
                        }
                    }
                }

            }
        }

//        System.out.println("------------");

//        System.out.println(matches);
        return matches;
    }


    /**
     * Get characters in the sequence
     * @return sequence string
     */
    public String getSequence() {
        return this.seq1;
    }

    /**
     * Get list of Node IDs representing current sequence after alignment, null for gaps in sequence
     * @return list of Node IDS
     */
    public List<Integer> getStringIndexes() {
        return this.stringIndexes;
    }


    /**
     * Get list of Node IDs that new sequence is being matched to, null for no match
     * @return list of Node IDS
     */
    public List<Integer> getNodeIndexes() {
        return this.nodeIndexes;
    }


//    /**
//     * Returns the indexes of all predecessor nodes
//     * @param node the node to get predecessors for
//     * @param nodeIDtoIndex a mapping of IDs to index positions
//     * @return list of indexes of all predecessor nodes
//     */
//
//    public List<Integer> getPrevIndices(Node node, Map<Integer, Integer> nodeIDtoIndex) {
//        List<Integer> prev = new ArrayList<Integer>();
//
//
//        List<Integer> edgeList = new ArrayList<Integer>();
//
//        for (Integer edge : node.getInEdges().keySet()) {
//            edgeList.add(edge);
//        }
//
//        for (Integer predID : edgeList) {
//            prev.add(nodeIDtoIndex.get(predID));
//        }
//
//        if (prev.size() == 0) {
//            prev.add(-1);
//        }
//
//        return prev;
//    }

    /**
     * Set up the initial sizes and values for the scores matrix, the matrices to
     * record the optimal moves through the score matrix, and the two mappings from
     * node ID to index and index to node ID
     * @return Object containing scores matrix, backSeq matrix, backGraph matrix, nodeID to
     * index, and index to node ID
     */
    public List<Object> initialiseDyncamicProgrammingData(int l1, int l2) {
//        l1 = this.seq1.length();
//        Integer l2 = this.seq2.length();
        Map<Integer, Integer> nodeIDToIndex = new HashMap<Integer, Integer>();
        Map<Integer, Integer> nodeIndexToID = new HashMap<Integer, Integer>();

//        for (Integer index : this.graph.getNodeDict().keySet()) {
//            nodeIDToIndex.put(this.graph.getNodeIDList().get(index), index);
//            nodeIndexToID.put(index, this.graph.getNodeIDList().get(index));
//        }

        // Create score matrix
        double[][] scores = new double[l1+1][l2+1];

        // Fill top row with gap penalties
        for (int i = 0; i < l2; i++) {
            scores[0][i+1] = openGapPenalty + (i * extendGapPenalty);
        }


        for (int i = 0; i < l1; ++i) {
            scores[i+1][0] = openGapPenalty + (i * extendGapPenalty);
        }

        int[][] backStrIdx = new int[l1 + 1][l2 + 1];

        for (int i = 0; i <l2; i++){
            backStrIdx[0][i+1] = i;
        }

        int[][] backGraphIdx = new int[l1 + 1][l2 + 1];

        for (int i = 0; i <l1; i++){
            backGraphIdx[i + 1][0] = i;
        }





        List<Object> initialisedData = new ArrayList<Object>();
        initialisedData.add(nodeIDToIndex);
        initialisedData.add(nodeIndexToID);
        initialisedData.add(scores);
        initialisedData.add(backStrIdx);
        initialisedData.add(backGraphIdx);

        return initialisedData;

    }

    /**
     * Trace back through the scores matrix, building the optimal alignment
     * @param scores
     * @param backStrIdx
     * @param backGrphIdx
     * @return
     */

//    public List<List<Integer>> backtrack(double[][] scores, int[][] backStrIdx, int[][] backGrphIdx) {
    public HashProfile backtrack(double[][] scores, int[][] backStrIdx, int[][] backGrphIdx) {


        //TODO: Scores only recording one optimal path - see simpleProfile vs simplerProfile 19/08/2016

        int besti = scores.length;
        int bestj = scores[0].length;
        String seq1output = "";
        String seq2output = "";
        besti -= 1;
        bestj -= 1;
        List<Integer> terminalIndices = new ArrayList<Integer>();


        double bestScore = scores[besti][bestj];

        for (int i = 0; i < terminalIndices.size(); i++) {
            double score = scores[terminalIndices.get(i)][bestj];
            if (score > bestScore) {
                bestScore = score;
                besti = terminalIndices.get(i);
            }
        }

        List<Integer> matches = new ArrayList<Integer>();
        List<Integer> strIndexes = new ArrayList<Integer>();

        int curstrIdx = 0;
        int curnodeIdx = 0;


        while (!(besti == 0 && bestj == 0)) {
            // Get the i and j coordinates of the cell to move to
            int nexti = backGrphIdx[besti][bestj];
            int nextj = backStrIdx[besti][bestj];

            curstrIdx = bestj - 1;

//            besti = (besti == 0) ? 1 : besti;
            curnodeIdx = besti - 1;

            if (nextj != bestj) {
                strIndexes.add(0, curstrIdx);
            } else {
                strIndexes.add(0, -1);

            }
            if (nexti != besti) {
                matches.add(0, curnodeIdx);

            } else {
                matches.add(0, -1);
            }

            besti = nexti;
            bestj = nextj;
        }

        // Fill out the remaining indexes of each profile
        while (matches.get(0) > 0 && curnodeIdx > 0){
            matches.add(0, curnodeIdx - 1);
            strIndexes.add(0, -1);
            curnodeIdx -= 1;

        }

        while (strIndexes.get(0) > 0 && curstrIdx > 0) {
            strIndexes.add(0, curstrIdx - 1);
            matches.add(0, -1);
            curstrIdx -= 1;

        }


        List<List<Integer>> matchesIndex = new ArrayList<List<Integer>>();

        matchesIndex.add(matches);
        matchesIndex.add(strIndexes);

        //TODO: Make calling gaps a function
        List<Integer> gapPos = new ArrayList<Integer>();
        List<Integer> gapPos2 = new ArrayList<Integer>();


        for (int i = 0; i < matches.size(); i++){
            if (matches.get(i) == -1) {
                gapPos.add(i);
            }
        }

        for (int i = 0; i < strIndexes.size(); i++){
            if (strIndexes.get(i) == -1){
                gapPos2.add(i);
            }

        }

        if (gapPos.size() > 0){
            profile1.addGaps(gapPos);

        }

        if (gapPos2.size() > 0) {
            profile2.addGaps(gapPos2);
        }


//        for (Integer index: matchesIndex.get(1)){
//            if (index == null){
//                profile1.addGaps(index);
//            }
//
//        }

        HashProfile updatedProfile = new HashProfile(profile1, profile2);
//        System.out.println(updatedProfile.toString());
//        System.out.println("strIndex " + strIndexes);
//        System.out.println("matchs " + matchesIndex);

//        for (Integer pos: matchesIndex.get(1)){
//            if (pos == null){
//                seq1output += "-";
//            }
//            else{
//                seq1output += seq1.charAt(pos);
//
//            }
//
//
//            }
//
//        for (Integer pos: matchesIndex.get(0)) {
//            if (pos == null) {
//                seq2output += "-";
//            } else {
//                seq2output += seq2.charAt(pos);
//
//            }
//        }

//        System.out.println(seq1output);
//        System.out.println(seq2output);


        return updatedProfile;
//        return matchesIndex;
    }

    /**
     * Align a new sequence to the PO Graph
     * @return
     */


    public HashProfile alignSeqToGraph() {
//        public List<List<Integer>> alignSeqToGraph() {



            int l1 = this.seq1.length();
        int l2 = this.seq2.length();


        List<Object> initialisedData = this.initialiseDyncamicProgrammingData(l1, l2);

        double[][] scores = (double[][]) initialisedData.get(2);
        int[][] backStrIdx = (int[][]) initialisedData.get(3);
        int[][] backGraphIdx = (int[][]) initialisedData.get(4);
        double[][] insertCost = new double[l1 + 1][l2 + 1];
        double[][] deleteCost = new double[l1 + 1][l2 + 1];

        for (double[] row : insertCost)
            Arrays.fill(row, this.openGapPenalty);

        for (double[] row : deleteCost)
            Arrays.fill(row, this.openGapPenalty);


        // Array to keep track of whether we should insert
        boolean[] inserted = new boolean[l2 + 1];
        Arrays.fill(inserted, false);


        double[] insscores = new double[l2 + 2];

        for (int i = 0; i < this.seq1.length(); i++) {

            // Get character of node
            char pbase = seq1.charAt(i);

            // Get predecessors of node
//            List<Integer> predecessors = this.getPrevIndices(this.graph.getNodeDict().get(this.graph.getNodeIDList().get(i)), nodeIDToIndex);


            // Get array of scores for matching current node with each position in sequence
            double[] matchPoints = new double[seq2.length()];

            if (MEA){
                matchPoints = this.getMEAMatchScore(i, seq2);
            }

            else {

                matchPoints = this.getMatchScore(pbase, seq2);

            }

            // Get array of scores equal to previous row + the cost of a deletion at each position
            double[] deleteScore = Arrays.copyOfRange(scores[i], 1, scores[0].length);
            deleteScore = addArrays(deleteScore, Arrays.copyOfRange(deleteCost[i], 1, deleteCost[0].length));


            // Array to store the best possible score for deleting at each position
            int[] bestDelete = new int[l2];

            // Fill array with first predecessor position
            Arrays.fill(bestDelete, i);



            double[] matchScore = addArrays(Arrays.copyOfRange(scores[i], 0, scores[0].length - 1), matchPoints);
            int[] bestMatch = new int[l2];
            Arrays.fill(bestMatch, i);


//            for (int j = 1; j < predecessors.size(); j++) {
//                int predecessor = predecessors.get(j);
//                double[] newDeleteScore = Arrays.copyOfRange(scores[predecessors.get(j) + 1], 1, scores[0].length);
//
//
//                for (int k = 0; k < bestDelete.length; k++) {
//                    int bestDeleteElement = newDeleteScore[k] > deleteScore[k] ? predecessor + 1 : bestDelete[k];
//                    bestDelete[k] = bestDeleteElement;
//
//                    int bestDeleteScoreElement = max(newDeleteScore[k], deleteScore[k]);
//                    deleteScore[k] = bestDeleteScoreElement;
//                }
//
//
//                double[] newMatchScore = addArrays(Arrays.copyOfRange(scores[predecessors.get(j) + 1], 0, scores[0].length - 1), matchPoints);
//
//                for (int n = 0; n < bestMatch.length; n++) {
//                    int bestMatchElement = newMatchScore[n] > matchScore[n] ? predecessor + 1 : bestMatch[n];
//                    bestMatch[n] = bestMatchElement;
//
//                    int bestMatchScoreElement = max(newMatchScore[n], matchScore[n]);
//                    matchScore[n] = bestMatchScoreElement;
//                }
//            }


            boolean[] deleted = new boolean[deleteScore.length];

            for (int k = 0; k < deleteScore.length; k++) {
                deleted[k] = (deleteScore[k]) >= matchScore[k];
                scores[i + 1][k + 1] = max(deleteScore[k], matchScore[k]);
            }


            for (int l = 1; l < scores[0].length; l++) {
                scores[i + 1][l] = max(deleteScore[l - 1], matchScore[l - 1]);

                if (deleted[l - 1]) {
                    backGraphIdx[i + 1][l] = bestDelete[l - 1];
                    backStrIdx[i + 1][l] = l;
                    deleteCost[i + 1][l] = extendGapPenalty;

                } else {
                    backGraphIdx[i + 1][l] = bestMatch[l - 1];
                    backStrIdx[i + 1][l] = l - 1;
                }
            }

            double[] scoreRow = scores[i + 1];
            double[] insertRow = insertCost[i + 1];
            int[] backStrRow = backStrIdx[i + 1];
            Arrays.fill(inserted, false);
            insscores = addArrays(scoreRow, insertRow);

            for (int n = 0; n < l2; n++) {
                if (insscores[n] >= scoreRow[n + 1]) {
                    scores[i + 1][n + 1] = insscores[n];
                    scoreRow[n + 1] = insscores[n];

                    insertCost[i + 1][n + 1] = this.extendGapPenalty;
                    insertRow[n + 1] = this.extendGapPenalty;
                    backStrIdx[i + 1][n + 1] = n;
                    backStrRow[n + 1] = n;
                    inserted[n + 1] = true;
                    insscores[n + 1] = scores[i + 1][n + 1] + insertCost[i + 1][n + 1];
                }
            }

            for (int o = 0; o < inserted.length; o++) {
                if (inserted[o]) {
                    insertCost[i + 1][o] = this.extendGapPenalty;
                    deleteCost[i + 1][o] = this.openGapPenalty;
                    backGraphIdx[i + 1][o] = i + 1;
                }
            }
        }
//        System.out.println("Scores matrix is ");
//        PairwisePairHMM.printMatrix(scores);

        return backtrack(scores, backStrIdx, backGraphIdx);

    }

    /**
     * Align a new sequence to the PO Graph
     * @return
     */


    public HashProfile alignProfileToProfile() {
//        public List<List<Integer>> alignProfileToProfile() {

        int profile1Length = this.profile1.getProfileArray().size();
        int profile2Length = this.profile2.getProfileArray().size();

        int l1 = profile1Length;
        int l2 = profile2Length;

        //TODO: Fix up how initialisedData is returned
        List<Object> initialisedData = this.initialiseDyncamicProgrammingData(l1, l2);

        double[][] scores = (double[][]) initialisedData.get(2);
        int[][] backStrIdx = (int[][]) initialisedData.get(3);
        int[][] backGraphIdx = (int[][]) initialisedData.get(4);
        double[][] insertCost = new double[l1 + 1][l2 + 1];
        double[][] deleteCost = new double[l1 + 1][l2 + 1];

        for (double[] row : insertCost)
            Arrays.fill(row, this.openGapPenalty);

        for (double[] row : deleteCost)
            Arrays.fill(row, this.openGapPenalty);


        // Array to keep track of whether we should insert
        boolean[] inserted = new boolean[l2 + 1];
        Arrays.fill(inserted, false);


        double[] insscores = new double[l2 + 2];

        for (int i = 0; i < profile1Length; i++) {

//            profile1.getProfileMatrix()[][i];



            double[] matchPoints;

            if (MEA){
                matchPoints = this.getMEAMatchScore(i, profile2);
            }

            else {

                matchPoints = this.getProfileMatchScore(profile1, profile2, i);

            }



            // Get character of node
//            char pbase = seq1.charAt(i);

            // Get predecessors of node
//            List<Integer> predecessors = this.getPrevIndices(this.graph.getNodeDict().get(this.graph.getNodeIDList().get(i)), nodeIDToIndex);


            // Get array of scores for matching current node with each position in sequence
//            double[] matchPoints = new double[seq2.length()];

//            if (MEA){
//                matchPoints = this.getMEAMatchScore(i, seq2);
//            }
//
//            else {
//
//                matchPoints = this.getMatchScore(pbase, seq2);
//
//            }

            // Get array of scores equal to previous row + the cost of a deletion at each position
            double[] deleteScore = Arrays.copyOfRange(scores[i], 1, scores[0].length);
            deleteScore = addArrays(deleteScore, Arrays.copyOfRange(deleteCost[i], 1, deleteCost[0].length));


            // Array to store the best possible score for deleting at each position
            int[] bestDelete = new int[l2];

            // Fill array with first predecessor position
            Arrays.fill(bestDelete, i);



            double[] matchScore = addArrays(Arrays.copyOfRange(scores[i], 0, scores[0].length - 1), matchPoints);
            int[] bestMatch = new int[l2];
            Arrays.fill(bestMatch, i);


//            for (int j = 1; j < predecessors.size(); j++) {
//                int predecessor = predecessors.get(j);
//                double[] newDeleteScore = Arrays.copyOfRange(scores[predecessors.get(j) + 1], 1, scores[0].length);
//
//
//                for (int k = 0; k < bestDelete.length; k++) {
//                    int bestDeleteElement = newDeleteScore[k] > deleteScore[k] ? predecessor + 1 : bestDelete[k];
//                    bestDelete[k] = bestDeleteElement;
//
//                    int bestDeleteScoreElement = max(newDeleteScore[k], deleteScore[k]);
//                    deleteScore[k] = bestDeleteScoreElement;
//                }
//
//
//                double[] newMatchScore = addArrays(Arrays.copyOfRange(scores[predecessors.get(j) + 1], 0, scores[0].length - 1), matchPoints);
//
//                for (int n = 0; n < bestMatch.length; n++) {
//                    int bestMatchElement = newMatchScore[n] > matchScore[n] ? predecessor + 1 : bestMatch[n];
//                    bestMatch[n] = bestMatchElement;
//
//                    int bestMatchScoreElement = max(newMatchScore[n], matchScore[n]);
//                    matchScore[n] = bestMatchScoreElement;
//                }
//            }


            boolean[] deleted = new boolean[deleteScore.length];

            for (int k = 0; k < deleteScore.length; k++) {
                deleted[k] = (deleteScore[k]) >= matchScore[k];
                scores[i + 1][k + 1] = max(deleteScore[k], matchScore[k]);
            }


            for (int l = 1; l < scores[0].length; l++) {
                scores[i + 1][l] = max(deleteScore[l - 1], matchScore[l - 1]);

                if (deleted[l - 1]) {
                    backGraphIdx[i + 1][l] = bestDelete[l - 1];
                    backStrIdx[i + 1][l] = l;
                    deleteCost[i + 1][l] = extendGapPenalty;

                } else {
                    backGraphIdx[i + 1][l] = bestMatch[l - 1];
                    backStrIdx[i + 1][l] = l - 1;
                }
            }

            double[] scoreRow = scores[i + 1];
            double[] insertRow = insertCost[i + 1];
            int[] backStrRow = backStrIdx[i + 1];
            Arrays.fill(inserted, false);
            insscores = addArrays(scoreRow, insertRow);

            for (int n = 0; n < l2; n++) {
                if (insscores[n] >= scoreRow[n + 1]) {
                    scores[i + 1][n + 1] = insscores[n];
                    scoreRow[n + 1] = insscores[n];

                    insertCost[i + 1][n + 1] = this.extendGapPenalty;
                    insertRow[n + 1] = this.extendGapPenalty;
                    backStrIdx[i + 1][n + 1] = n;
                    backStrRow[n + 1] = n;
                    inserted[n + 1] = true;
                    insscores[n + 1] = scores[i + 1][n + 1] + insertCost[i + 1][n + 1];
                }
            }

            for (int o = 0; o < inserted.length; o++) {
                if (inserted[o]) {
                    insertCost[i + 1][o] = this.extendGapPenalty;
                    deleteCost[i + 1][o] = this.openGapPenalty;
                    backGraphIdx[i + 1][o] = i + 1;
                }
            }
        }
//        System.out.println("Scores matrix is ");
//        PairwisePairHMM.printMatrix(scores);

        return backtrack(scores, backStrIdx, backGraphIdx);

    }




    /**
     * Method to add two arrays together by adding the contents of each respective index in each array together.
     * For example
     * [1,4,9,6] + [2,2,3] = [3,7,12,6]
     * @param firstArray first array to add
     * @param secondArray second array to add
     * @return int array with the add
     *
     *
     */

    //TODO: Move this helper method to the matrix helper methods
    public static double[] addArrays ( double[] firstArray, double[] secondArray){

        double[] shorterArray = (firstArray.length < secondArray.length ? firstArray: secondArray);
        double[] longerArray = (firstArray.length > secondArray.length ? firstArray: secondArray);
        double[] finalArray = new double[longerArray.length];


        for (int i = 0; i < shorterArray.length; i++) {
            finalArray[i] = (firstArray[i] + secondArray[i]);
        }
        for (int i = shorterArray.length; i < longerArray.length; i++){
            finalArray[i] = longerArray[i];
        }
        return finalArray;
    }


    //TODO: Remove printMatrix to helper Class
    public static void printMatrix(double[][] matrix) {


        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", matrix[i][j]) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }

    //TODO: Remove printMatrix to helper Class
    public static void printMatrix(int[][] matrix) {


        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(String.format("  %-10s", matrix[i][j]) + "|");
            }
            System.out.println();

        }

        System.out.println();

    }
}