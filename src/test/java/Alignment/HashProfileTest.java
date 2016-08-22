package Alignment;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

/**
 * Created by gabe on 22/08/2016.
 */
public class HashProfileTest {

    private HashProfile setupSingle(){
        HashProfile hashProfile = new HashProfile("AGCAG");
        return hashProfile;
    }

    @Test
    public void testprintSingle() throws Exception {
        HashProfile hashProfile = setupSingle();
        assertEquals("AGCAG", hashProfile.toString());
    }



    @Test
    public void testFillProfileArray() throws Exception {

    }
    

    @Test
    public void testFillProfileArray1() throws Exception {

    }

    @Test
    public void testAddGaps() throws Exception {

    }

    @Test
    public void testGetProfileArray() throws Exception {

    }

    @Test
    public void testGetSequences() throws Exception {
        

    }

    @Test
    public void testToString() throws Exception {

    }

    @Test
    public void testPrintProfileArray() throws Exception {

    }
}