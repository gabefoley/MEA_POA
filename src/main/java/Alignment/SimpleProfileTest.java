package Alignment;

import SubstitutionModels.Blosum62;

/**
 * Created by gabe on 15/08/2016.
 */
public class SimpleProfileTest {
    public static void main(String[] args) {
        SimpleProfile simpleProfile = new SimpleProfile("AGCAG");

        SimpleProfile first = new SimpleProfile("AGCAG", "AG-TG");

        SimpleProfile second = new SimpleProfile("GCAG", "-CAG");

//        POGraphAlignment alignment = new POGraphAlignment(first, second, -2, -1, Blosum62.getMatrix(), false);

    }

}
