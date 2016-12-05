package Alignment;

import Alignment.Utilities.MatrixUtils;
import Alignment.Utilities.Sequence;
import SubstitutionModels.SubstitutionMatrix;
import jdk.internal.util.xml.impl.Pair;
//import org.biojava.nbio.alignment.Alignments;
//import org.biojava.nbio.alignment.SimpleGapPenalty;
//import org.biojava.nbio.alignment.template.GapPenalty;
//import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
//import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
//import org.biojava.nbio.core.alignment.template.Profile;
//import org.biojava.nbio.core.alignment.template.SequencePair;
//import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
//import org.biojava.nbio.core.sequence.ProteinSequence;
//import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
//import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
//import org.biojava.nbio.core.util.ConcurrencyTools;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Helper method to test running alignments
 */
public class RunAlignments {

    public static void main(String[] args) throws Exception {


        // Setup pairwise inputs
//        String pairwiseQuery = "TTACG";
//        String pairwiseTarget = "TAG";

        // Setup pairwise inputs
//        String pairwiseQuery = "TTACGACC";
//        String pairwiseTarget = "CCTTAACC";

//        String pairwiseQuery = "V";
//        String pairwiseTarget = "N";

//        double tau = 0.01;
//        double epsilon = 0.1;
//        double delta = 0.2;
//        double eta = 0.1111111111;
//        double openPenalty = calculateOpenPenalty(tau, epsilon, delta, eta);
//        double extensionPenalty = calculateExtensionPenalty(epsilon, eta);
//        double calculatedEpsilon = calculateEpsilon(extensionPenalty, eta);

//        System.out.println("Open penalty is " + openPenalty);
//        System.out.println("Extension penalty is " + extensionPenalty);
//        System.out.println("Calculated epsilon is " + calculatedEpsilon);


        //Profile too large problem
//        String pairwiseQuery = "MYSFPNSFRFGWSQAGFQSEMGTPGSEDPNTDWYKWVHDPENMAAGLVSGDLPENGPGYWGNYKTFHDNAQKMGLKIARLNSEW" +
//                "SRQFPNPLPRPQNFDESKQDVTEVEINENELKRLDEYANKDALNHYREIFKDLKSRGLYFIQNMYHWPLPLWLHDPIRVRRGDFTGPSGWLSTRTVYEF" +
//                "ARFSAYTAWKFDDLVDEYSTMNEPNVVGGLGYVGVKSGFPPGYLSFELSRRAMYNIIQAHARAYDGIKSVSKKPVGIIYANSSFQPLTDKDMEAVEMAE" +
//                "NDNRWWFFDAIIRGEITRGNEKIVRDDLKGRLDWIGVNYYTRTVVKRTEKGYVSLGGYGHGCERNSVSLAGLPTSDFGWEFFPEGLYDVLTKYWNRYHL" +
//                "APGYSVNKEVEEKTNLWKDETAQEAFIHHWSFIARRYKGISSTHLSFNLINEPPFPDPQIXSVEDHNSLIKRTITEIRKIDPERLIIIDGLGYGNIPVD" +
//                "YMYVTENGIADDADYQRPYYLVSHVYQVHRAINSGADVRGYLHWSLADNYEWASGFSMRFGLLKVDYNTKRLYWRPSALVYREIATNGAITDEIEHLNS" +
//                "VPPVKPLRH";
//        String pairwiseTarget = "IPRWRGFNLLEAFSIKSTGNFKEEDFLWXAQWDFNFVRIPXCHLLWSDRGNPFIIREDFFEKIDRVIFWGEKYGIHICISLHR" +
//                "APGYSVNKEVEEKTNLWKDETAQEAFIHHWSFIARRYKGISSTHLSFNLINEPPFPDPQIXSVEDHNSLIKRTITEIRKIDPERLIIIDGLGYGNIPVD" +
//                "DLTIENTVQSCRGYIPFSVTHYKAEWVDSKDFPVPEWPNGWHFGEYWNREKLLEHYLTWIKLRQKGIEVFCGEXGAYNKTPHDVVLKWLEDLLEIFKTL" +
//                "NIGFALWNFRGPFGILDSERKDVEYEEWYGHKLDRKXLELLRKY";

//        String pairwiseQuery = "PRANGAFRANGAVANNN";
//        String pairwiseTarget = "PRANGAKKKKKKKFRANGAKKKKVANNNNN";

//        String pairwiseQuery = "PPGRR";
//        String pairwiseTarget = "NNPP";

//        HashProfile profile1 = new HashProfile(pairwiseQuery);
//        HashProfile profile2 = new HashProfile(pairwiseTarget);

//        double gapOpen = openPenalty;
//        double gapExtend = extensionPenalty;

        int gapOpen = 10;
        int gapExtend = 1;




//        double tau = 0.0001;
//        double epsilon = 0.6;
//        double delta = 0.05;
//
//        String pairwiseQuery = "TACATACACCC";
//        String pairwiseTarget = "CATGATCATCAT";

//        checkMEA();


//        String pairwiseQuery = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE";
//        String pairwiseTarget = "MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK";


//        String pairwiseQuery = "MLSRNATSSYFLGWQEYEKNPYHEVHNTNGIIQMGLAENQLCFDLLESWLAKNPEAAAFKKNGESIFAELALFQDYHGLPAFKK" +
//                "AMVDFMAEIRGNKVTFDPNHLVLTAGATSANETFIFCLADPGEAVLIPTPYYPGFDRDLKWRTGVEIVPIHCTSSNGFQITETALEEAYQEAEKRNLRV" +
//                "KGVLVTNPSNPLGTTMTRNELYLLLSFVEDKGIHLISDEIYSGTAFSSPSFISVMEVLKDRNCDENSEVWQRVHVVYSLSKDLGLPGFRVGAIYSNDDM" +
//                "VVAAATKMSSFGLVSSQTQHLLSAMLSDKKLTKNYIAENHKRLKQRQKKLVSGLQKSGISCLNGNAGLFCWVDMRHLLRSNTFEAEMELWKKIVYEVHL" +
//                "NISPGSSCHCTEPGWFRVCFANLPERTLDLAMQRLKAFVG";
//
//        String pairwiseTarget = "LFNTAHGGNIREPATVLGISPDQLLDFSANINPLGMPVSVKRALIDNLDCIERYPDADYFHLHQALARHHQVPASWILAGNGE" +
//                "TESIFTVASGLKPRRAMIVTPGFAEYGRALAQSGCEIRRWSLREADGWQLTDAILEALTPDLDCLFLCTPNNPTGLLPERPLLQAIADRCKSLNINLIL" +
//                "DEAFIDFIPHETGFIPALKDNPHIWVLRSLTKFYAIPGLRLGYLVNSDDAAMARMRRQQMPWSVNALAALAGEVALQDSAWQQATWHWLREEGARFYQA" +
//                "LCQLPLLTVYPGRANYLLLRCEREDIDLQRRLLTQRILIRSCANYPGLDSRYYRVAIRSAAQNERLLAALRNVL";


//        double emissionX = 0.05;
//        double emissionY = 0.05;

        //Setup MSA inputs
        int multiGapOpen = -4;
        int multiGapExtend = -2;


        String multiSeq1 = "AARS";
        String multiSeq2 = "AAMGRS";
        String multiSeq3 = "AARN";

//        String multiSeq1 = "WAA";
//        String multiSeq2 = "NWA";
//        String multiSeq3 = "NA";

//        String multiSeq1 = "CAA";
//        String multiSeq2 = "TCA";
//        String multiSeq3 = "TA";




//        String multiSeq4 = "FRANKKKWAN";

//        String multiSeq1 = "GGPAWKKKKKKKRE";
//        String multiSeq2 = "AWENNP";
//        String multiSeq3 = "GRPWESASS";
//
//        String multiSeq4 = "DDPGAEERED";

        String[] multiArray = new String[3];
        multiArray[0] = multiSeq1;
        multiArray[1] = multiSeq2;
        multiArray[2] = multiSeq3;


//        String multiSeq1 = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGEKKKKKKKKK";
//        String multiSeq2 = "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK";
//        String multiSeq3 = "MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK";
//        String multiSeq4 = "MHIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK";



        SubstitutionMatrix blosum62 = new SubstitutionMatrix("blosum62");
        SubstitutionMatrix blosum62Probs = new SubstitutionMatrix("blosum62Probs");
        SubstitutionMatrix exampleModel = new SubstitutionMatrix("exampleModel");
        SubstitutionMatrix blosum62estimatedlambda = new SubstitutionMatrix("blosum62estimatedlambda");
        SubstitutionMatrix blosum62LatestProbs = new SubstitutionMatrix("blosum62LatestProbs");
        SubstitutionMatrix blosum62EstimatedWithX = new SubstitutionMatrix("blosum62EstimatedWithX");
        SubstitutionMatrix hmmocModel = new SubstitutionMatrix("hmmocModel");
        //Setup HMM inputs
        double tau = 0.1;
        double epsilon = 0.1;
        double delta = 0.2;

        String pairwiseQuery = "TTACG";
        String pairwiseTarget = "TAG";
        double emissionX = 0.25;
        double emissionY = 0.25;


        // Run pairwise alignments
//        runBioJavaPairwise(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend);
//        runPairwise(pairwiseQuery, pairwiseTarget, gapOpen, gapExtend, blosum62);


//        runViterbiPairwiseOnProfile(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, emissionX, emissionY, exampleModel, "nucleotide");
//        runMEAPairwiseOnProfile(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, emissionX, emissionY, hmmocModel, "nucleotide");

//        runDurbinComparison(exampleModel);
        // Run Baum Welch pairwise
//        runBaumWelch();
//        exampleAlignment();



//        runBWDefault(pairwiseQuery, pairwiseTarget, blosum62estimatedlambda);
//        runBWMulti(multiSeq1, multiSeq2, multiSeq3, blosum62EstimatedWithX);
//
        runHMMOCComparison(exampleModel);
//        runHMMOCComparisonUnderflow(blosum62EstimatedWithX);

//        runBWTest();
//        runBWAlign();
//        runIceCream();


        // Run multiple sequence alignments
//        runBioJavaMSA(multiSeq1, multiSeq2, multiSeq3, multiGapOpen, multiGapExtend);
//        runMSA(multiSeq1, multiSeq2, multiSeq3, multiGapOpen, multiGapExtend, blosum62);




        // this one
//        runViterbiMSA(multiSeq1, multiSeq2, multiSeq3, tau, epsilon, delta, emissionX, emissionY, blosum62EstimatedWithX, "protein");
//        runViterbiMulti(multiArray, tau, epsilon, delta, emissionX, emissionY, blosum62EstimatedWithX);
//        runMEAMSA(multiSeq1, multiSeq2, multiSeq3, tau, epsilon, delta, emissionX, emissionY, blosum62Probs);
//        runMEAMSAMulti(multiArray, tau, epsilon, delta, emissionX, emissionY, blosum62EstimatedWithX);







    }

    private static void runBWAlign(){
        String seq1 = "TA";
        String seq2 = "CCTA";
        String[] seqs = new String[]{seq1, seq2};


        double [] start = {0.33333, 0.33333, 0.33333};

        double [][] transition = {
                {.1, .45, .45},
                {.1, .8, .1},
                {.1, .1, .8}};

        double [][] emission = {
                {0.4, 0.2, 0.2, 0.2},
                {0.25, 0.25, 0.25, 0.25},
                {0.25, 0.25, 0.25, 0.25}};

        BaumWelch bw = new BaumWelch(seqs, start, transition, emission, "nucleotide");
        System.out.println("done");
        start = bw.getStart();
        transition = bw.getTransition();
        emission = bw.getEmission();
        SubstitutionMatrix subMatrix = new SubstitutionMatrix("blosum62EstimatedWithX");
        PairHMM  pairHMM = new PairHMM(seqs, start, transition, emission, subMatrix, false, "nucleotide");
        HashProfile alignment = pairHMM.getViterbiAlignment();
//        HashProfile mea = pairHMM.getMEAAlignment(1);

        System.out.println(alignment);


    }

    private static void runIceCream(){
        String seq = "ACCACACAACGCCGGGAGGGCGAGGGACCACAA";
        String[] seqs = new String[]{seq};

        double [] start = {0.5, 0.5, 0};

        double [][] transition = {
                {.8, .1, .1},
                {.1, .8, .1},
                {0,0,1}};

        double [][] emission = {
                {0.7, 0.2, 0.1},
                {0.1, 0.2, 0.7},
                {0,0,0}};

        BaumWelch bw = new BaumWelch(seqs, start, transition, emission, "nucleotide");

    }


    private static void runBWTest(){
        String seq1 = "CGAAC";
        String[] seqs = new String[]{seq1};



        double [] start = {0.33333, 0.33333, 0.33333};

        double [][] transition = {
                {.1, .45, .45},
                {.1, .8, .1},
                {.1, .1, .8}};

        double [][] emission = {
                {0.4, 0.3, 0.3},
                {0.33333, 0.33333, 0.3333},
                {0.33333, 0.33333, 0.33333}};

        BaumWelch bw = new BaumWelch(seqs, start, transition, emission, "nucleotide");

    }

    private static void runBWTestRabiner(){
        String seq1 = "GAGGGGGAAGAAAGAGA";
        String[] seqs = new String[]{seq1};



        double [] start = {0.5, 0.5};

        double [][] transition = {
                {.85, .15},
                {.12, .88}};

        double [][] emission = {
                {0.8, 0.1, 0.1},
                {0, 0, 0.1}};

        BaumWelch bw = new BaumWelch(seqs, start, transition, emission, "nucleotide");

    }

    private static void runDurbinComparison(SubstitutionMatrix subMatrix){

        double [] start = {.5, .2, .2};

        double [][] transition = {
                {.5, .2, .2, .1},
                {.8, .1, 0, 0.1},
                {.8, 0, .1, 0.1}};


        double [][] emission = {
                {0.25, 0.25, 0.25, 0.25},
                {0.25, 0.25, 0.25, 0.25},
                {0.25, 0.25, 0.25, 0.25}};



        String pairwiseQuery = "TAG";
////
        String pairwiseTarget = "TTACG";

//        String pairwiseQuery = "T";
////
//        String pairwiseTarget = "TG";

        String[] seqArray = {pairwiseQuery, pairwiseTarget};

        String type = "nucleotide";

        PairHMMUnderflow  pairHMM = new PairHMMUnderflow(seqArray, start, transition, emission, subMatrix, false, type);
//        HashProfile alignment = pairHMM.getViterbiAlignment();

//        System.out.println(alignment);
        HashProfile mea = pairHMM.getMEAAlignment(1);
        System.out.println(mea);




    }

    private static void runHMMOCComparison(SubstitutionMatrix subMatrix){

        double [] start = {.6, 0.4, 0.4};

        double [][] transition = {
                {0.5, 0.2, 0.2, 0.1},
                {0.8, 0.1, 0, 0.1},
                {0.8, 0, 0.1, 0.1}};


        double [][] emission = {
                {0.25, 0.25, 0.25, 0.25},
                {0.25, 0.25, 0.25, 0.25},
                {0.25, 0.25, 0.25, 0.25}};



        String pairwiseQuery = "TA";
//
        String pairwiseTarget = "CCTA";

        String[] seqArray = {pairwiseQuery, pairwiseTarget};

        String type = "nucleotide";

        PairHMMUnderflow  pairHMM = new PairHMMUnderflow(seqArray, start, transition, emission, subMatrix, true, type);
//        HashProfile alignment = pairHMM.getViterbiAlignment();

//        System.out.println(alignment);
        HashProfile mea = pairHMM.getMEAAlignment(1);
        System.out.println(mea);




    }

    private static void checkMEA(){
//        String query = "TAG";
//        String target = "TTACG";

        String query = "NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNS";
        String target = "PLALLLDSSLEGEFDLVQRIIYEVDDPSLPNDEGITALHNAVCAGHTEIVKFLVQFGVNVNAADSDGWTPLHCAASCNNVQVCKFLVESGAAVFAMTYSDMQTAADKCEEMEEGYTQCSQFLYGVQEKMGIMNKGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGYVPRNLLGLYP";

        String[] seqs = new String[]{query, target};

        double tau = 0.1;
        double epsilon = 0.1;
        double delta = 0.2;

        Double emissionX = 0.25;
        Double emissionY = 0.25;

        double [] start = {.95, 0.025, 0.025};

//        double [][] transition = {
//                {.5, .2, .2, 0.1},
//                {.8, .1, 0, 0.1},
//                {.8, 0, .1, 0.1}};

        double [][] transition = {
                {.95, .025, .025, 0.00009},
                {.4, .6, 0, 0.00009},
                {.4, 0, .6, 0.00009}};

//        double [][] emission = {
//                {0.25, 0.25, 0.25, 0.25},
//                {0.25, 0.25, 0.25, 0.25},
//                {0.25, 0.25, 0.25, 0.25}};

//        double [][] emission = {
//                {.25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25},
//                {.25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25},
//                {.25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25}};



        double [][] emission = {
                {.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05},
                {.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05},
                {.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05}};

        String type = "protein";
        SubstitutionMatrix subMatrix = new SubstitutionMatrix("blosum62EstimatedWithX");
//        SubstitutionMatrix subMatrix = new SubstitutionMatrix("exampleModel");
        PairHMMUnderflow alignment = new PairHMMUnderflow(seqs, start, transition, emission, subMatrix, false, type);

//        PairHMMUnderflow alignment = new PairHMMUnderflow(seqs, tau, epsilon, delta, emissionX, emissionY, subMatrix);

//        PairHMM alignment = new PairHMM(seqs, start, transition, emission, subMatrix, false);

//        PairHMM alignment = new PairHMM(seqs, tau, epsilon, delta, emissionX, emissionY, subMatrix);
        HashProfile mea = alignment.getMEAAlignment(1);
        System.out.println(mea);

    }

    private static void runHMMOCComparisonUnderflow(SubstitutionMatrix subMatrix){

        double [] start = {.95, 0.025, 0.025};

        double [][] transition = {
                {.95, .025, .025, 0.00009},
                {.4, .6, 0, 0.00009},
                {.4, 0, .6, 0.00009}};


        double [][] emission = {
                {.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05},
                {.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05},
                {.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05}};



        String pairwiseQuery = "PPGNK";

        String pairwiseTarget = "GNK";

//        String pairwiseQuery = "SISDTVKRAREAFNSGKTRSLQFRIQQLEALQRMINENLKSISGALASDLGKNEWTSYYEEVAHVLEELDTTIKELPDWAEDEPVAKTRQTQQDDLYIHSEPLGVVLVIGAWNYPFNLTIQPMVGAVAAGNAVILKPSEVSGHMADLLATLIPQYMDQNLYLVVKGGVPETTELLKERFDHIMYTGSTAVGKIVMAAAAKHLTPVTLELGGKSPCYVDKDCDLDVACRRIAWGKFMNSGQTCVAPDYILCDPSIQNQIVEKLKKSLKDFYGEDAKQSRDYGRIINDRHFQRVKGLIDNQKVAHGGTWDQSSRYIAPTILVDVDPQSPVMQEEIFGPVMPIVCVRSLEEAIQFINQREKPLALYVFSNNEKVIKKMIAETSSGGVTANDVIVHITVPTLPFGGVGNSGMGAYHGKKSFETFSHRRSCLVKSLLNEEAHKARYPPSPA";
//////
//        String pairwiseTarget = "MTVEPFRNEPIETFQTEEARRAMREALRRVREEFGRHYPLYIGGEWVDTKERMVSLNPSAPSEVVGTTAKAGKAEAEAALEAAWKAFKTWKDWPQEDRSRLLLKAAALMRRRKRELEATLVYEVGKNWVEASADVAEAIDFIEYYARAALRYRYPAVEVVPYPGEDNESFYVPLGAGVVIAPWNFPVAIFTGMIVGPVAVGNTVIAKPAEDAVVVGAKVFEIFHEAGFPPGVVNFLPGVGEEVGAYLVEHPRIRFINFTGSLEVGLKIYEAAGRLAPGQTWFKRAYVETGGKNAIIVDETADFDLAAEGVVVSAYGFQGQKCSAASRLILTQGAYEPVLERVLKRAERLSVGPAEENPDLGPVVSAEQERKVLSYIEIGKNEGQLVLGGKRLEGEGYFIAPTVFTEVPPKARIAQEEIFGPVLSVIRVKDFAEALEVANDTPYGLTGGVYSRKREHLEWARREFHVGNLYFNRKITGALVGVQPFGGFKLSGTNAKTGALDYLRLFLEMKAVAERF";

//        String pairwiseQuery = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE";
//
//        String pairwiseTarget = "MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK";

        String sequence3 = "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK";

        Sequence seq1 = new Sequence("first", pairwiseQuery);
        Sequence seq2 = new Sequence("second", pairwiseTarget);
        Sequence seq3 = new Sequence("third", sequence3);

        Sequence[] seqs = new Sequence[3];
        seqs[0] = seq1;
        seqs[1] = seq2;
        seqs[2] = seq3;

        String type = "protein";




        String[] seqArray = {pairwiseQuery, pairwiseTarget};

//        PairHMMUnderflow  pairHMM = new PairHMMUnderflow(seqArray, start, transition, emission, subMatrix, false);

        PairHMMUnderflow  seqsPairHMM = new PairHMMUnderflow(seqs, start, transition, emission, subMatrix, false, type);


//        HashProfile alignment = pairHMM.getViterbiAlignment();

        HashProfile seqsAlignmnet = seqsPairHMM.getMEAAlignment(1);
//        HashProfile mea = pairHMM.getMEAAlignment(1);

//        System.out.println(alignment);
        System.out.println(seqsAlignmnet);

//        System.out.println(mea);



    }


    private static void runTopsComparison(SubstitutionMatrix subMatrix){

        double [] start = {.68, 0.16, 0.16};

        double [][] transition = {
                {.95, .025, .025, 0.00009},
                {.4, .6, 0, 0.00009},
                {.4, 0, .6, 0.00009}};


        double [][] emission = {
                {.25, .25, .25, .25},
                {.25, .25, .25, .25,},
                {.25, .25, .25, .25,}};




        String [] seqArray = {
                "TACATACACCCCACCGGGACGACGATAGTAGTAGTAGAAAA",
                "CATGATCATCATGACGATATATGAAAAA"
        };
        String type = "nucleotide";

        PairHMM  pairHMM = new PairHMM(seqArray, start, transition, emission, subMatrix, false, type);
        HashProfile alignment = pairHMM.getViterbiAlignment();

        System.out.println(alignment);



    }


//    private static void topsParameters(){
//        double[] start = {0.6814756989, 0.00008615339902, 0.1591759622 };
//
//        double[][] transition = { }
//    }
//}
    private static void runBaumWelch(){
        double [] start = {.95, 0.025, 0.025};

        double [][] transition = {
                {.95, .025, .025, 0.00009},
                {.4, .6, 0, 0.00009},
                {.4, 0, .6, 0.00009}};


        double [][] emission = {
                {.25, .25, .25, .25},
                {.25, .25, .25, .25,},
                {.25, .25, .25, .25,}};


        String [] seqArray = {
                "TA",
                "CCTA"
        };



//        BaumWelchMulti bw = new BaumWelchMulti(seqArray, start, transition, emission, "nucleotide");
        SubstitutionMatrix subMatrix = new SubstitutionMatrix("hmmocModel");
//
//        MatrixUtils.printMatrix(bw.getEmission());
//        MatrixUtils.printMatrix(bw.getTransition());

        PairHMMUnderflow initialPairhmm = new PairHMMUnderflow(seqArray, start, transition, emission, subMatrix, true, "nucleotide");

//        PairHMMUnderflow pairhmm = new PairHMMUnderflow(seqArray, start, bw.getTransition(), bw.getEmission(), subMatrix, true, "nucleotide");


        HashProfile alignment = initialPairhmm.getViterbiAlignment();
        System.out.println(alignment);


    }

    private static void runBWDefault(String pairwiseQuery, String pairwiseTarget, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("\nViterbi alignment with Baum Welch parameter estimation: ");
        String[] seqArray = new String[2];


        seqArray[0] = pairwiseQuery;
        seqArray[1] = pairwiseTarget;


        PairHMM pairHMM = new PairHMM(seqArray, subMatrix, "protein");

        HashProfile alignment = pairHMM.getViterbiAlignment();

        System.out.println(alignment);

    }



    public static void runBWMulti(String multiSeq1, String multiSeq2, String multiSeq3, SubstitutionMatrix subMatrix){
        String[] seqArray = new String[3];

        seqArray[0] = multiSeq1;
        seqArray[1] = multiSeq2;
        seqArray[2] = multiSeq3;

        PairHMM pairHMM = new PairHMM(seqArray, subMatrix, "protein");
        HashProfile alignment = pairHMM.getViterbiAlignment();

        System.out.println(alignment);


    }

//    private static void runBioJavaPairwise(String pairwiseQuery, String pairwiseTarget, int gapOpen, int gapExtend)
//            throws Exception{
//        System.out.println("--------------------------");
//        System.out.println("\nBioJava NW Alignment: ");
//
//        GapPenalty penalty = new SimpleGapPenalty(gapOpen, gapExtend);
//
//        // Setup BioJava Pairwise alignment
//        PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments.getPairwiseAligner(
//                new ProteinSequence(pairwiseQuery, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
//                new ProteinSequence(pairwiseTarget, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
//                Alignments.PairwiseSequenceAlignerType.GLOBAL,
//                penalty, SubstitutionMatrixHelper.getBlosum62());
//
//
//        SequencePair<ProteinSequence, AminoAcidCompound>
//                alignment = aligner.getPair();
//
//
//        System.out.println(alignment);
//
//    }

    private static void runPairwise(String pairwiseQuery, String pairwiseTarget, int gapOpen, int gapExtend, SubstitutionMatrix subMatrix){
        System.out.println("--------------------------");
        System.out.println("NW Alignment: ");

        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);


        Alignment alignment = new Alignment(pairwiseQuery, pairwiseTarget,
                -1 * gapOpen, -1 * gapExtend, subMatrix, false);

        System.out.println(alignment.getUpdatedProfile().printSeqs());

//        System.out.println(alignment.getUpdatedProfile());


    }



    private static void runViterbiPairwiseOnProfile(String pairwiseQuery, String pairwiseTarget, double tau, double epsilon,
                                           double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type){
        System.out.println("--------------------------");
        System.out.println("\nViterbi Pairwise using PairHMM: ");

        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);



        PairHMMUnderflow pairHMM = new PairHMMUnderflow(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
//        HashProfile alignment = pairHMM.getViterbiAlignment();
        HashProfile mea = pairHMM.getMEAAlignment(1);
        System.out.println(mea.printSeqs());

    }

    private static void runMEAPairwiseOnProfile(String pairwiseQuery, String pairwiseTarget, double tau, double epsilon,
                                       double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type ){
        System.out.println("--------------------------");
        System.out.println("\nMaximum Expected Accuracy Alignment using PairHMM: ");

        HashProfile profile1 = new HashProfile(pairwiseQuery);
        HashProfile profile2 = new HashProfile(pairwiseTarget);


        PairHMMUnderflow pairHMM = new PairHMMUnderflow(pairwiseQuery, pairwiseTarget, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile alignment = pairHMM.getMEAAlignment(1);

        System.out.println(alignment.printSeqs());

//        Alignment meaAlignment = pairHMM.getMEAAlignment();



//        PairHMM pairHMM = new PairHMM(pairwiseQuery, pairwiseTarget, tau, epsilon, delta);

    }

//    private static void runBioJavaMSA(String multiSeq1, String multiSeq2, String multiSeq3,
//                                      int multiGapOpen, int multiGapExtend) throws CompoundNotFoundException {
//        System.out.println("--------------------------");
//        System.out.println("\nBioJava Multiple Sequences Alignment: ");
//
//        ProteinSequence firstProt = new ProteinSequence(multiSeq1, AminoAcidCompoundSet.getAminoAcidCompoundSet());
//        ProteinSequence secondProt = new ProteinSequence(multiSeq2, AminoAcidCompoundSet.getAminoAcidCompoundSet());
//        ProteinSequence thirdProt = new ProteinSequence(multiSeq3, AminoAcidCompoundSet.getAminoAcidCompoundSet());
////        ProteinSequence fourthProt = new ProteinSequence(multiSeq4, AminoAcidCompoundSet.getAminoAcidCompoundSet());
//
//        GapPenalty penalty = new SimpleGapPenalty(multiGapOpen,multiGapExtend);
//
////        Alignments.getMultipleSequenceAlignment()
////
////        PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments.getPairwiseAligner(
////                new ProteinSequence(pairwiseQuery, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
////                new ProteinSequence(pairwiseTarget, AminoAcidCompoundSet.getAminoAcidCompoundSet()),
////                Alignments.PairwiseSequenceAlignerType.GLOBAL,
////                penalty, SubstitutionMatrixHelper.getBlosum62());
//
//        List<ProteinSequence> list = new ArrayList<ProteinSequence>();
//        list.add(firstProt);
//        list.add(secondProt);
//        list.add(thirdProt);
////        list.add(fourthProt);
//        Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(list, penalty,
//                SubstitutionMatrixHelper.getBlosum62());
//
//        ConcurrencyTools.shutdown();
//
//        System.out.println(profile);
//
//
//    }

//    private static void runMSA(String multiSeq1, String multiSeq2, String multiSeq3,
//                               int multiGapOpen, int multiGapExtend, SubstitutionMatrix subMatrix){
//        System.out.println("--------------------------");
//        System.out.println("\nMultiple Sequences Alignment: ");
//        HashProfile first = new HashProfile(multiSeq1);
//        HashProfile second = new HashProfile(multiSeq2);
//        HashProfile third = new HashProfile(multiSeq3);
////        HashProfile fourth = new HashProfile(multiSeq4);
//
//
//
//
//        System.out.println("\nFirst alignment:");
//        Alignment firstAlignment = new Alignment(first, second, multiGapOpen, multiGapExtend, subMatrix, false);
//        HashProfile firstProfile = firstAlignment.getUpdatedProfile();
//
//        System.out.println(firstProfile.printSeqs());
//
//        System.out.println("\nSecond alignment:");
//
//
//        Alignment secondAlignment = new Alignment(firstProfile, third, multiGapOpen, multiGapExtend,  subMatrix, false);
//        HashProfile secondProfile = secondAlignment.getUpdatedProfile();
//        System.out.println(secondProfile.printSeqs());
//
//
////        System.out.println("\nThird alignment:");
////
////        Alignment thirdAlignment = new Alignment(secondProfile, fourth, multiGapOpen, multiGapExtend,  subMatrix, false);
////        HashProfile thirdProfile = thirdAlignment.getUpdatedProfile();
////
////        System.out.println(thirdProfile.printSeqs());
//
//    }



    private static void runViterbiMSA(String multiSeq1, String multiSeq2, String multiSeq3, double tau, double epsilon,
                                      double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type){

        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
//        HashProfile fourth = new HashProfile(multiSeq4);


        System.out.println("--------------------------");
        System.out.println("\nViterbi Multiple Sequences Alignment: ");



        System.out.println("\nFirst alignment:");
        PairHMM firstAlignment = new PairHMM(first, second, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile firstProfile = firstAlignment.getViterbiAlignment();
        System.out.println(firstProfile.printSeqs());

        System.out.println("\nSecond alignment:");


        PairHMM secondAlignment = new PairHMM(firstProfile, third, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile secondProfile = secondAlignment.getViterbiAlignment();
        System.out.println(secondProfile.printSeqs());


//        System.out.println("\nThird alignment:");
//
//        PairHMM thirdAlignment = new PairHMM(secondProfile, fourth, tau, epsilon, delta, emissionX, emissionY, subMatrix);
//        HashProfile thirdProfile = thirdAlignment.getViterbiAlignment();
//        System.out.println(thirdProfile.printSeqs());

    }

//    private static void runViterbiMulti(String[] seqArray, double tau, double epsilon,
//                                        double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type){
//        System.out.println("--------------------------");
//        System.out.println("\nViterbi Alignment using seqArray method: ");
//
//        PairHMM alignment = new PairHMM(seqArray, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
//        HashProfile profile = alignment.getViterbiAlignment();
//
//        System.out.println("Alignment");
//        System.out.println(profile.printSeqs());
//
//
//
//
//    }



    private static void runMEAMSA(String multiSeq1, String multiSeq2, String multiSeq3, double tau, double epsilon,
                                  double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type){
        System.out.println("--------------------------");
        System.out.println("\nMaximum Expected Accuracy Multiple Sequences Alignment: ");

        HashProfile first = new HashProfile(multiSeq1);
        HashProfile second = new HashProfile(multiSeq2);
        HashProfile third = new HashProfile(multiSeq3);
//        HashProfile fourth = new HashProfile(multiSeq4);

        System.out.println("\nFirst alignment:");
        PairHMM firstAlignment = new PairHMM(first, second, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile firstProfile = firstAlignment.getMEAAlignment(1);
        System.out.println(firstProfile.printSeqs());

        System.out.println("\nSecond alignment:");


        PairHMM secondAlignment = new PairHMM(firstProfile, third, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile secondProfile = secondAlignment.getMEAAlignment(1);
        System.out.println(secondProfile.printSeqs());

//
//        System.out.println("\nThird alignment:");
//
//        PairHMM thirdAlignment = new PairHMM(secondProfile, fourth, tau, epsilon, delta, emissionX, emissionY, subMatrix);
//        HashProfile thirdProfile = thirdAlignment.getMEAAlignment();
//        System.out.println(thirdProfile.printSeqs());



    }

    private static void runMEAMSAMulti(String[] seqArray, double tau, double epsilon,
                                  double delta, double emissionX, double emissionY, SubstitutionMatrix subMatrix, String type) {

        System.out.println("--------------------------");
        System.out.println("\nMEA Alignment using seqArray method: ");

        PairHMM alignment = new PairHMM(seqArray, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile profile = alignment.getMEAAlignment(1);

        System.out.println("Alignment");
        System.out.println(profile.printSeqs());


    }

    public static double calculateEta(String[] seqArray){
        double length = 0;
        for (String seq : seqArray){
            length += seq.length();
        }
        double meanLength = length / seqArray.length;

        double eta = 1 / (meanLength + 1);
        return eta;
    }

    public static double calculateOpenPenalty(double tau, double epsilon, double delta, double eta){



        double openPenalty = - Math.log((delta * (1 - epsilon - tau))/((1-eta) * (1-2*delta - tau)));

        return openPenalty;


    }

    public static double calculateExtensionPenalty(double epsilon, double eta ){

        double oneMinusEta = 1 - eta;
        double epsilondivided = epsilon / oneMinusEta;
        double negativelog = - Math.log(epsilondivided);

//        System.out.println("One minus eta:" + oneMinusEta);
//        System.out.println("epsilonDivided" + epsilondivided);
//        System.out.println("negativelog" + negativelog);

        double extensionPenalty = - Math.log((epsilon) / (1 - eta));

        return extensionPenalty;
    }

    public static double calculateEpsilon(double extensionPenalty, double eta){
        double epsilon = (1 - eta) / (Math.pow(Math.E, extensionPenalty));

        return epsilon;
    }

    public static String[] getBalibaseSeqs() throws  IOException{
        String dir = "/Users/gabe/Dropbox/PhD/Programs/bench1.0/bali3/in";
        File file = new File(dir);

        String[] inputs = file.list();

        for (String input: inputs){
            System.out.println("Reading from: " + input);

            FileReader fr = new FileReader(dir + "/" +  input);

            BufferedReader br = new BufferedReader(fr);
            String[] seqs;
            if (! br.readLine().startsWith(">")){
//                seqs(br.readLine();

            }
        }

        return inputs;

    }


    public static void exampleAlignment(){

        double tau = 0.1;
        double epsilon = 0.1;
        double delta = 0.2;
        double emissionX = 0.25;
        double emissionY = 0.25;



        HashProfile a = new HashProfile("PARDW");
        HashProfile b = new HashProfile("PANDWST");
        HashProfile c = new HashProfile("SADFGPRETR");
        HashProfile d = new HashProfile("SKKDETR");
        HashProfile e = new HashProfile("RVWNMP");
        HashProfile f = new HashProfile("RVMPATR");


        String type = "protein";
        SubstitutionMatrix subMatrix = new SubstitutionMatrix("blosum62EstimatedWithX");


        PairHMMUnderflow pairHMM1 = new PairHMMUnderflow(a, b, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile alignment1 = pairHMM1.getViterbiAlignment();
        System.out.println("Alignment 1");
        System.out.println(alignment1);
        System.out.println();


        PairHMMUnderflow pairHMM2 = new PairHMMUnderflow(c, d, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile alignment2 = pairHMM2.getViterbiAlignment();
        System.out.println("Alignment 2");

        System.out.println(alignment2);
        System.out.println();


        PairHMMUnderflow pairHMM3 = new PairHMMUnderflow(e, f, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile alignment3 = pairHMM3.getViterbiAlignment();
        System.out.println("Alignment 3");

        System.out.println(alignment3);
        System.out.println();


        PairHMMUnderflow pairHMM4 = new PairHMMUnderflow(alignment2, alignment3, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        HashProfile alignment4 = pairHMM4.getViterbiAlignment();
        System.out.println("Alignment 4");

        System.out.println(alignment4);
        System.out.println();


        PairHMMUnderflow pairHMM5 = new PairHMMUnderflow(alignment1, alignment4, tau, epsilon, delta, emissionX, emissionY, subMatrix, type);
        System.out.println("Alignment 5");

        HashProfile alignment5 = pairHMM5.getViterbiAlignment();


        System.out.println(alignment5);
        System.out.println();

    }

}
