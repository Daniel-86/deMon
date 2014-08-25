package demon;

import java.io.IOException;

/**
 * User: daniel
 * Date: 4/06/13
 * Time: 01:32 AM
 */
public class DeMon {

    public DeMonParser parser;
    public DeMonPrinter printer;
    public Timer timer;
    public Atomo[] atoms;

    String trj_type = "NEW";
    String whichGeom = "NORMAL";
    
    String title;
    int charge = 0;
    int multiplicity;
    int n_elec;
    String matdia = "DSYEV";
    String startGuess = "HUECKEL";
    String scfTheory = "DEFAULT";
    double scfTol = 1e-5;
    boolean scfAdj = false;
    boolean cdxc = true;
    String fileToWrite = "MOLDEN";
    int scfMax = 100;
    boolean diis = true;
    boolean port = true;
    boolean l_shift = false;

    public DeMon() {
        try {
            parser = new DeMonParser();
            printer = new DeMonPrinter();
        } catch (IOException e) {
            System.out.println("Error deMon constructor:" + e.getMessage());
        }
    }




    public void initMPIExecution() {

    }




    public void printCopyright() {
        String copyright = "" +
                " ******************************************************************************\n" +
                " ******************************************************************************\n" +
                " ***                                                                        ***\n" +
                " ***                              PROGRAM deMon2k                           ***\n" +
                " ***                                                                        ***\n" +
                " ***                        (VERSION 4.1.4, July 2012)                      ***\n" +
                " ***                                                                        ***\n" +
                " ***     COPYRIGHT (C) BY THE INTERNATIONAL DEMON DEVELOPERS COMMUNITY      ***\n" +
                " ***                                                                        ***\n" +
                " ***    AUTHORS: A.M. Koster, G. Geudtner, P. Calaminici, M.E. Casida,      ***\n" +
                " ***             J. Carmona-Espindola, V.D. Dominguez, R. Flores-Moreno,    ***\n" +
                " ***             G.U. Gamboa, A. Goursot, T. Heine, A. Ipatov, F. Janetzko, ***\n" +
                " ***             J.M. del Campo, J.U. Reveles, J.M. Vasquez-Perez, A. Vela, ***\n" +
                " ***             B. Zuniga-Gutierrez and D.R. Salahub                       ***\n" +
                " ***                                                                        ***\n" +
                " ******************************************************************************\n" +
                " ******************************************************************************\n" +
                "\n\n";
        try {
            printer.printMessage(copyright);
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }





    public void startTimer() {
        timer = new Timer("FIRST TIMING");
        timer.start();
        timer.doMark(false);
    }






    public void doSymmetryAnalysis() {
        Symmetry symmetryDriver = new Symmetry();

        symmetryDriver.init();
        symmetryDriver.molsym(atoms);
//        symmetryDriver.findSymmetry(atoms);
        symmetryDriver.determineOrderPointGroup();
        symmetryDriver.generateSymmetryList(atoms);
        symmetryDriver.generetaSymmetryTransformationMatrix();
//        symmetryDriver.transformDummyAtoms();
        Globals.COORDINATE_SET = 1;
        Globals.nVibrations = 3 * atoms.length - 5;
    }





    public void printInputGeometry() throws IOException {
        String header;

        if(trj_type.equals("NEW") || trj_type.equals("RESTART")) {
            header = "\n\n *** GEOMETRY ***\n";
            switch(whichGeom) {
                case "ORIGINAL_O":
                    header += "\n FILE deMon.mol WILL BE WRITTEN";
                    break;
                case "ORIGINAL_PES":
                    header += "\n FILE deMon.mkl WILL BE WRITTEN";
                    break;
                case "RESTART_OPT_STEP":
                    header += "\n FILE deMon.wfn WILL BE WRITTEN";
                    break;
                default:
                    header += "\n INPUT ORIENTATION in BOHR\n";
                    String noH = "  No."; // 6u
                    String atomH = "ATOM"; // 8i
                    String xH = "X"; // 26u
                    String yH = "Y"; // 41u
                    String zH = "Z"; // 56u
                    String zAtomH = "Z-ATOM"; // 65
                    String massH = "MASS"; // 74u
                    String optxyzH = "OPTXYZ"; // 80u
                    header += "\n  NO.  ATOM        X              Y              Z        Z-ATOM   MASS   OPTXYZ\n";
                    printer.printCoords(header, atoms);
                    break;
            }
        }

    }






    public void printInitResume() throws IOException {
        String message = "";

        message += " " + Parser.repeat("*", 79) + "\n\n";
        message += " " + title + "\n\n";
        message += " " + Parser.repeat("*", 79) + "\n\n";

        message += "\n *** MOLECULE DATA ***\n";
        message += "\n MOLECULAR CHARGE: " + charge ;
        message += "\n MOLECULAR MULTIPLICITY: " + multiplicity ;
        message += "\n NUMBER OF ATOMS: " + parser.atomos.size() ;
        message += "\n NUMBER OF ELECTRONS: " + n_elec + "\n\n";

        message += " *** SCF OPTIONS ***\n";
//        if(true)
//            message += "\n SCF THEORY: RESTRICTED KOHN-SHAM (RKS)";
        switch(scfTheory) {
            case "NONE":
                message += "\n SCF THEORY: NOT USED";
                break;
            case "GUESSONLY":
                message += "\n SCF THEORY: GUESS ONLY CALCULATION";
                break;
            case "SCFMAX":
                message += "\n SCF THEORY: ENERGY ONLY CALCULATION (NO SCF)";
                break;
            case "OKS_ROKS":
                message += "\n SCF THEORY: RESTRICTED OPEN KOHN-SHAM (ROKS)";
                break;
            case "OKS_UKS":
                message += "\n SCF THEORY: UNRESTRICTED KOHN-SHAM (UKS)";
                break;
            default:
                message += "\n SCF THEORY: RESTRICTED KOHN-SHAM (RKS)";
                break;
        }
        message += "\n DIAGONALIZER: " + matdia ;
        switch(startGuess) {
            case "HUECKEL":
                message += "\n START GUESS: HUECKEL DENSITY";
                break;
            case "FERMI":
                message += "\n START GUESS: FERMI DENSITY";
                break;
            case "PROJCT":
                message += "\n START GUESS: PROJECTION";
                break;
            case "RESTART":
                message += "\n START GUESS: RESTART DENSITY";
                break;
            default:
                message += "\n START GUESS: CORE DENSITY";
                break;
        }
        message += "\n SCF TOLERANCE: " + scfTol ;
        if(scfAdj)
            message += "\n SCF TOLERANCE WILL BE TIGHTEN";
        message += "\n POTENTIAL DENSITY: " + (cdxc? "FITTED DENSITY": "ORBITAL DENSITY");
        message += "\n NUMBER OF SCF ITERATIONS: " + scfMax ;
        message += "\n DIIS ACCELERATION: " + (diis? "ON": "OFF");
        if(l_shift)
            message += "\n AQUI VA ALGO";
        message += "\n\n";

        printer.printMessage(message);
    }

    public void establishAtomsArray() {
        atoms = new Atomo[parser.atomos.size()];
        parser.atomos.toArray(atoms);
    }

    public void printSymmetryAnalisis() throws IOException {
        String header;

        if(trj_type.equals("NEW") || trj_type.equals("RESTART")) {
            header = "\n\n *** GEOMETRY ***\n";
            switch(whichGeom) {
                case "ORIGINAL_O":
                    header += "\n FILE deMon.mol WILL BE WRITTEN";
                    break;
                case "ORIGINAL_PES":
                    header += "\n FILE deMon.mkl WILL BE WRITTEN";
                    break;
                case "RESTART_OPT_STEP":
                    header += "\n FILE deMon.wfn WILL BE WRITTEN";
                    break;
                default:
                    header += "\n INPUT ORIENTATION in BOHR\n";
                    String noH = "  No."; // 6u
                    String atomH = "ATOM"; // 8i
                    String xH = "X"; // 26u
                    String yH = "Y"; // 41u
                    String zH = "Z"; // 56u
                    String zAtomH = "Z-ATOM"; // 65
                    String massH = "MASS"; // 74u
                    String optxyzH = "OPTXYZ"; // 80u
                    header += "\n  NO.  ATOM        X              Y              Z        Z-ATOM   MASS   OPTXYZ\n";
                    printer.printCoords(header, atoms);
                    break;
            }
        }
    }


//    public void generateOrbitalPointers(Atomo[] atoms) {
//        for(int i_atom = 1, i_sto = 0; i_sto < n_sto; i_sto++)
//    }
}
