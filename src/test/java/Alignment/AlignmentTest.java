package Alignment;

import SubstitutionModels.SubstitutionMatrix;
import org.junit.Test;

import java.util.ArrayList;

/**
 * Created by gabe on 22/08/2016.
 *
 *
 */

public class AlignmentTest {
//    private HashProfile first;
//    private HashProfile second;
//    private HashProfile third;
//    private HashProfile fourth;
//    private HashProfile firstJoined;
//    private HashProfile secondJoined;

    private SubstitutionMatrix blosum62 = new SubstitutionMatrix("blosum62");




    @Test
    public void initialiseData(){

        HashProfile first = new HashProfile("WPG");
        HashProfile second = new HashProfile("PPG");
        HashProfile third = new HashProfile("AAPG");
        HashProfile fourth = new HashProfile("-WWG");

        HashProfile firstJoined = new HashProfile(first, second);
        HashProfile secondJoined = new HashProfile(third, fourth);
    }

    @Test
    public void testFirstAlignment(){
        HashProfile first = new HashProfile("WPG");
        HashProfile second = new HashProfile("PPG");
        HashProfile third = new HashProfile("AAPG");
        HashProfile fourth = new HashProfile("-WWG");

        HashProfile firstJoined = new HashProfile(first, second);
        HashProfile secondJoined = new HashProfile(third, fourth);        Alignment firstAlignment = new Alignment(first, second, -2, -1, blosum62, false);


    }

    @Test
    public void testSecondAlignment(){
        HashProfile first = new HashProfile("WPG");
        HashProfile second = new HashProfile("PPG");
        HashProfile third = new HashProfile("AAPG");
        HashProfile fourth = new HashProfile("-WWG");

        HashProfile firstJoined = new HashProfile(first, second);
        HashProfile secondJoined = new HashProfile(third, fourth);        Alignment secondAlignment = new Alignment(third, fourth, -2, -1, blosum62, false);
    }



    @Test
    public void testProfileProfileAlignment(){
        HashProfile first = new HashProfile("WPG");
        HashProfile second = new HashProfile("PPG");
        HashProfile third = new HashProfile("AAPG");
        HashProfile fourth = new HashProfile("-WWG");

        HashProfile firstJoined = new HashProfile(first, second);
        HashProfile secondJoined = new HashProfile(third, fourth);
        Alignment joinedAlignment = new Alignment(secondJoined, firstJoined, -2, -1, blosum62, false);

        ArrayList<Integer> expectedStringMatches = new ArrayList<Integer>();
        expectedStringMatches.add(0);
        expectedStringMatches.add(1);
        expectedStringMatches.add(2);
        expectedStringMatches.add(3);


        ArrayList<Integer> expectedNodeMatches = new ArrayList<Integer>();
        expectedNodeMatches.add(-1);
        expectedNodeMatches.add(0);
        expectedNodeMatches.add(1);
        expectedNodeMatches.add(2);

        // Assert that the matches array is correct
//        assertEquals(joinedAlignment.getStringIndexes(), expectedStringMatches);
//        assertEquals(joinedAlignment.getNodeIndexes(), expectedNodeMatches);


    }

    @Test
    public void testNucleotideAlignment(){
        HashProfile first = new HashProfile("AATGTGGG");
        HashProfile second = new HashProfile("TTGGG");
        HashProfile third = new HashProfile("CCAATG");
        HashProfile fourth = new HashProfile("CAACCG");



        HashProfile firstJoined = new HashProfile(first, second);
        HashProfile secondJoined = new HashProfile(third, fourth);
        Alignment joinedAlignment = new Alignment(first, second, -2, -1, blosum62, false);



    }

}