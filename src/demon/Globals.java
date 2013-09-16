package demon;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * User: daniel
 * Date: 15/04/13
 * Time: 07:54 PM
 */
public class Globals {

    /**
     * globals
     */
    static int max_atom = 350; // maximum number of atoms.
    static int max_shl = 2000; // maximum number of shells.
    static int max_sto = 5000; // Maximum number of contracted STO orbitals.
    static int max_aux = 10000; // Maximum number of auxiliary functions.
    static int max_gto = 10000; // Maximum number of Gaussian functions.
    static int max_aux_set = 2500; // Maximum number of auxiliary function sets.
    static int max_aux_shl = 2500; // Maximum number of auxiliary function shells.
    static int max_con = 21; // Maximum degree of contraction.
    static int max_l_aux = 6; // Maximum L quantum number for auxiliary functions.
    static int max_l_bas = 4; // Maximum L quantum number for orbital basis.
    static int max_mom = 3; // Maximum rank (>0) of electrostatic moments.
    static int max_der = 1; // Maximum rank of analytic energy derivatives.
    static int max_l_ul = max_l_aux+max_der+2*max_l_bas; //  Maximum L quantum number upper limit.
    static int max_opt = 31; // Maximum number of options for input keywords.


    /**
     * mp_start.f
     */
    static int mp_maxi = 6; // maximum number of mpi index.
    static boolean mpp = false; // true if a parallel deMon run is requested.
    static boolean use_mpf = false; // true if the parallel part of a routine is used.
    static boolean used_mat = false; // true if an already distributed matrix is used.
    static String mp_frq_type = "XYZ"; //frequency parallelization type.
    static int me = 0; // number of local process.
    static int master = 0; // number of master process.
    static double mp_scf = 1.0; //scaling factor depending on the number of processes.
    static int nCPU = 1; // total number of processes (CPU).
    static boolean mp_master = true; // true if local process is the master process.
    static boolean mp_slave = false; // true if local process is a slave process.
    static boolean[] mp_ini = new boolean[mp_maxi]; // true if mpi load balancing is initialized.

    /**
     * machine.f
     */
    static double largest_exp = 0.0;
    static double smallest_exp = 0.0;

    /**
     * fileio.f
     */
    static BufferedWriter outFile = openFile("deMon.out");
    static BufferedWriter inFile = openFile("deMon.in");
    static BufferedWriter newFile = openFile("deMon.new");
    static BufferedWriter trjFile = openFile("deMon.trj");
    static BufferedWriter memFile = openFile("deMon.mem");
    static BufferedWriter rstFile = openFile("deMon.rst");
    private static BufferedWriter openFile(String path) {
        BufferedWriter file = null;
        try {
            file = new BufferedWriter(new FileWriter(path));
        } catch (IOException e) {

        }
        return file;
    }

    /**
     * input.f
     */

    /**
     * maxckh.f
     */
    public static void checkParamConsistency(String typeChking) {
        String errorMsg;
        switch(typeChking) {
            case "PARAMETER":
                if(max_atom < 3)
                    errorMsg = "\n *** PARAMETER ERROR: MAXATOM < 3 ***";
                if(max_atom > max_shl)
                    errorMsg = "\n *** PARAMETER ERROR: MAXATOM > MAXSHL ***";
                if(max_shl > max_sto)
                    errorMsg = "\n *** PARAMETER ERROR: MAXSHL > MAXSTO ***";
                if(max_sto > max_aux)
                    errorMsg = "\n *** PARAMETER ERROR: MAXSTO > MAXAUX ***";
                if(max_sto > max_gto)
                    errorMsg = "\n *** PARAMETER ERROR: MAXSTO > MAXGTO ***";
                if(max_con > max_gto)
                    errorMsg = "\n *** PARAMETER ERROR: MAXCON > MAXGTO ***";
                if(max_atom > max_aux_set)
                    errorMsg = "\n *** PARAMETER ERROR: MAXATOM > MAXAUXSET ***";
                if(max_aux_set > max_aux_shl)
                    errorMsg = "\n *** PARAMETER ERROR: MAXAUXSET > MAXAUXSHL ***";
                if(max_gto > max_sto*max_con)
                    errorMsg = "\n *** PARAMETER ERROR: MAXGTO > MAXSTO*MAXCON ***";
                if(max_aux > max_gto*(max_gto+1)/2)
                    errorMsg = "\n *** PARAMETER ERROR: MAXAUX > MAXGTO*(MAXGTO+1)/2 ***";
                if(max_l_bas > max_l_aux)
                    errorMsg = "\n *** PARAMETER ERROR: MAXLBAS > MAXLAUX ***";
                if(max_mom > max_l_aux)
                    errorMsg = "\n *** PARAMETER ERROR: MAXMOM > MAXLAUX ***";
                if(max_l_ul > 16)
                    errorMsg = "\n *** NO MACHINE PRECISION FOR GAMMA FUNCTION ***";
        }
    }

    /**
     * readin.f
     */

    /**
     * default.f
     */
    static int charge = 0; // Molecular charge.
    static boolean sym_geo = false; // True, if symmetry in geometry optimization is enforced.
    static boolean symmetry = false; // If true, use molecular symmetry.
    static boolean alignment = false; //
    static boolean modmap = false; //
    static boolean enantio = false; //
    static boolean uniform = false; //
    static int cublev = 6; //
    static int nblack = 0; //
    static int poptol = 50; //
    static String inptyp = "CARTESIAN"; // default geometry
    static boolean scanpes = false; //
    static boolean relaxing = true; //
    static String scan_end = "UNKNOWN"; //
    static String scan_lab = "UNKNOWN"; //
    static boolean embed = false; //
    static double fix_none = 1e20; //
    static double fix_cups = 1e3;
    static double fix_core = 1e2;
    static double fix_valence = 1.0;
    static String fit_typ = "ALL";
    static double zet_fix = fix_none;
    static boolean[] fix_auxis = new boolean[max_atom]; // If true, corresponding atom possess fixed functions.
    static int multip;

    /**
     * newdrv.f
     */
//    static void newdrv(String action) throws IOException {
//        switch(action) {
//            case "INPUT":
//                if(opt) {
//                    if(zbuild) {
//                        outFile.write("CONSTRuiR Z MATRIX");
//                    } else if(opttyp.equals("CARTESIAN")){
//
//                    }
//                }
//        }
//    }


    /**
     * flags.h
     */
    static boolean SPHORB;

    /**
     * settings
     */
    static double TINIEST = 1e-14;

    /**
     * declaraciones en parameter.h
     */
    static final int MAX_ATOM = 350;
    static final int MAX_AUX_SET = 2500;
    static final int MAX_GTO = 10000;
    static final int MAX_AUX_SHL = 2500;
    static final int MAX_SHELL = 2000;
    static final int MAX_STO = 5000;
    static final int MAX_TB_STO = 10000;
    static final int MAX_L_BAS = 4;
    static final int MAX_L_PTR = StrictMath.max(5, max_l_ul);


    /**
     * declaraciones pointer.h
     */
    static Integer [] lowerLimitShell = new Integer[MAX_ATOM];
    static Integer [] upperLimitShell = new Integer[MAX_ATOM];
    static Integer [] lowerLimitGTO = new Integer[MAX_SHELL];
    static Integer [] upperLimitGTO = new Integer[MAX_SHELL];
    static Integer [] lowerLimitSTO = new Integer[MAX_SHELL];
    static Integer [] upperLimitSTO = new Integer[MAX_SHELL];
    static Integer [] lowerAuxSet = new Integer[MAX_ATOM];
    static Integer [] upperAuxSet = new Integer[MAX_ATOM];
    static Integer [] lowerAuxShell = new Integer[MAX_AUX_SET];
    static Integer [] upperAuxShell = new Integer[MAX_AUX_SET];
    static Integer[][] shellPointer = new Integer[MAX_SHELL][4];
    static Integer[] N_Cartesian_Orbitals = new Integer[MAX_L_PTR];
    static Integer[] ll = new Integer[MAX_ATOM];
    static Integer[] ul = new Integer[MAX_ATOM];
    static Integer[][] STO_PTR = new Integer[MAX_STO][5];
    static Integer [] lowerlimitSPH = new Integer[MAX_AUX_SET];
    static Integer [] upperlimitSPH = new Integer[MAX_AUX_SET];
    static Integer[][] sphericPointer = new Integer[MAX_SHELL][4];
    static int lowerLimitMOS;
    static int upperLimitMOS;
    static int lowerLimitDOS;
    static int upperLimitDOS;


    /**
     * symmetry variables, constants, etc.
     */
    static final int MAX_C_AXIS_ROT_DEGREE = 20;
    static int MAX_C_AXIS = MAX_C_AXIS_ROT_DEGREE + 1;
    static final int MAX_S_AXIS_ROT_DEGREE = 2 * MAX_C_AXIS_ROT_DEGREE + 1;
    static int MAX_S_AXIS = MAX_C_AXIS_ROT_DEGREE * 2;
    static int MAX_MIRROR_PLANES = MAX_C_AXIS_ROT_DEGREE + 1;
    static int HIGEST_C_AXIS_DEGREE;
    static int HIGEST_S_AXIS_DEGREE;
    static Integer [] N_C_AXIS = new Integer[MAX_C_AXIS_ROT_DEGREE]; //NSEC
    static Integer [] N_S_AXIS = new Integer[MAX_S_AXIS_ROT_DEGREE]; //NSES
    static int N_MIRROR_PLANES;
    // verdadero si se encuentra un grupo puntual cubico como T, Th, Td, etc
    static boolean CUBIC;
    static boolean DGROUP;
    static boolean IGROUP;
    static boolean INVERS;
    static boolean LINEAR;
    static boolean MAXIS;
    static boolean OGROUP;
    static boolean SGROUP;
    static boolean SIGMAD;
    static boolean SIGMAH;
    static boolean SIGMAV;
    static boolean TGROUP;
    static int N_GRP_EQ_ATOMS;
    static int pointGroupOrder;
    static Integer [] equivalentAtomsGroup = new Integer[MAX_ATOM];
    static Integer [] groupLowerLimit = new Integer[MAX_ATOM];
    static Integer [] nEquivalentAtoms = new Integer[MAX_ATOM];
    static Integer [] equivalentAtomsGroupOrderedList = new Integer[MAX_ATOM];
    static Integer [] groupUpperLimit = new Integer[MAX_ATOM];
    static Double [][] inputToStandarMatrix = new Double[3][3];
    static Double [][] standardToInputMatrix = new Double[3][3];
    static double molecularMass;
    static double molecularCharge;
    static int nElectrons;
    static int spinMultiplicity;
    static Double [] translationVector = new Double[3];
    static double cartesianCoordinatesTolerance;
    static Double [] ain = new Double[MAX_ATOM];
    static boolean LINDEG;
    /**
     *   1, Cartesian (Z-Matrix).
     *   CCSET = 2, Cartesian standard orientation.
     *   CCSET = 3, Internal (Z-matrix).
     *   CCSET = 4, Cartesian input.
     *   CCSET = 5, Work coordinates.
     *   CCSET = 6, Step coordinates.
     *   CCSET = 7, Cartesian reference coordinates.
     */
    static int COORDINATE_SET;
    static int nVibrations;




    /**
     * declaracion de molecule.h
     */
    static int N_ATOM;
    static int [] nuclearCharges = new int[MAX_ATOM];
    static double [] atomicMass = new double[MAX_ATOM];
    static double [] gaussianContractionCoeficient = new double[MAX_GTO];
    static double [] gaussianExponent = new double[MAX_GTO + MAX_AUX_SET]; //
    static int N_GAUSSIAN = 0; //NGTO
    static int N_SHELLS;
    static int N_ORBITALS;
    static Integer[] MAX_L_FOR_ORB_SHELL = new Integer[MAX_ATOM];
    static int N_SPH;
    static int N_ORB;
    static int N_ORTH;


    /**
     * guessdrv
     */
    static int N_TB_ATOMS;
    static int N_TB_GRPS;
    static Integer[] TB_OF = new Integer[MAX_ATOM];
    static Integer[] TB_GRP_LOWER = new Integer[MAX_ATOM];
    static Integer[] TB_GRP_UPPER = new Integer[MAX_ATOM];
    static Integer[] TB_GRP_LIST = new Integer[MAX_ATOM];
    static Integer[] N_TB_EQU_ATOMS = new Integer[MAX_ATOM];
    static int N_TB_ORBITALS;
    static int N_TB_SHELLS;
    static Integer[] TB_ORBITALS_LOWER = new Integer[MAX_ATOM];
    static Integer[][] TB_ORBITAL_POINTER = new Integer[MAX_ATOM][5];
    static Integer[] TB_ORBITALS_UPPER = new Integer[MAX_ATOM];
    static Integer[][] TB_ATOM_POINTER = new Integer[MAX_ATOM][3];
    static String[] TB_ORBITALS_SYMBOLS = new String[MAX_STO];
    static Double[] TB_AO_CONTRACT_COEFF = new Double[MAX_TB_STO];

    static Double [][] aoccna = new Double[(int) StrictMath.pow(MAX_L_BAS+1, 2)][7];
    static Double [][] aoccnb = new Double[(int) StrictMath.pow(MAX_L_BAS+1, 2)][7];
    static Integer [] elcfga = new Integer[(int) StrictMath.pow(MAX_L_BAS + 1, 2)];
    static Integer [] elcfgb = new Integer[(int) StrictMath.pow(MAX_L_BAS + 1, 2)];


}
