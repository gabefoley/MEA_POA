package Alignment;

import SubstitutionModels.Blosum62;
import org.junit.Test;

import java.util.ArrayList;

import static org.junit.Assert.*;

/**
 * Created by gabe on 22/08/2016.
 *
 *
 */

public class SequenceAlignerTest {
//    private HashProfile first;
//    private HashProfile second;
//    private HashProfile third;
//    private HashProfile fourth;
//    private HashProfile firstJoined;
//    private HashProfile secondJoined;




    //TODO: Not working when aligning two profiles and the first profile is shorter
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
        HashProfile secondJoined = new HashProfile(third, fourth);        SequenceAligner firstAlignment = new SequenceAligner(first, second, -2, -1, Blosum62.getMatrix(), false);


    }

    @Test
    public void testSecondAlignment(){
        HashProfile first = new HashProfile("WPG");
        HashProfile second = new HashProfile("PPG");
        HashProfile third = new HashProfile("AAPG");
        HashProfile fourth = new HashProfile("-WWG");

        HashProfile firstJoined = new HashProfile(first, second);
        HashProfile secondJoined = new HashProfile(third, fourth);        SequenceAligner secondAlignment = new SequenceAligner(third, fourth, -2, -1, Blosum62.getMatrix(), false);
    }



    @Test
    public void testProfileProfileAlignment(){
        HashProfile first = new HashProfile("WPG");
        HashProfile second = new HashProfile("PPG");
        HashProfile third = new HashProfile("AAPG");
        HashProfile fourth = new HashProfile("-WWG");

        HashProfile firstJoined = new HashProfile(first, second);
        HashProfile secondJoined = new HashProfile(third, fourth);
        SequenceAligner joinedAlignment = new SequenceAligner(secondJoined, firstJoined, -2, -1, Blosum62.getMatrix(), false);

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
        assertEquals(joinedAlignment.getStringIndexes(), expectedStringMatches);
        assertEquals(joinedAlignment.getNodeIndexes(), expectedNodeMatches);


    }

    @Test
    public void testNucleotideAlignment(){
        HashProfile first = new HashProfile("AATGTGGG");
        HashProfile second = new HashProfile("TTGGG");
        HashProfile third = new HashProfile("CCAATG");
        HashProfile fourth = new HashProfile("CAACCG");



        HashProfile firstJoined = new HashProfile(first, second);
        HashProfile secondJoined = new HashProfile(third, fourth);
        SequenceAligner joinedAlignment = new SequenceAligner(first, second, -2, -1, Blosum62.getMatrix(), false);



    }

}