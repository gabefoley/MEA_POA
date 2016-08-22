package Alignment;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by gabe on 22/08/2016.
 */
public class ProfilePairHMMTest {


    public ProfilePairHMM createSimpleProfile(){

        HashProfile first = new HashProfile("TTG");
        HashProfile second = new HashProfile("CTG");
        HashProfile third = new HashProfile("AATG");
        HashProfile fourth = new HashProfile("-CCG");

        HashProfile firstJoined = new HashProfile(first, second);
        HashProfile secondJoined = new HashProfile(third, fourth);

        double tau = 0.1;
        double epsilon = 0.1;
        double delta = 0.2;



        ProfilePairHMM pairHMM = new ProfilePairHMM(firstJoined, secondJoined, tau, epsilon, delta);
        return pairHMM;

    }

    @Test
    public void testMatchArrayCorrectProfilePairHMM(){
        ProfilePairHMM pairHMM = createSimpleProfile();
        assertNotNull(pairHMM);
    }

    @Test
    public void testMEA(){
        ProfilePairHMM pairHMM = createSimpleProfile();
        

    }



}