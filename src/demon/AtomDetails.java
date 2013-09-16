package demon;

/**
 * User: daniel
 * Date: 8/04/13
 * Time: 08:23 PM
 */
public class AtomDetails {

    protected String symbol = null;
    protected Double mass;
    protected Double covalent;
    protected Double vander;
    protected String group = null;
    protected Integer atomicNumber;
    protected Double c6;
    public static final int ATMNUMBER = 0x01;
    public static final int MASS = 0x02;
    public static final int COVALENT = 0x04;
    public static final int VANDER = 0x06;
    public static final int GROUP = 0x08;

    public AtomDetails(Integer atomicNumber, String symbol, Double mass, Double covalent, Double vander, String group, Double c6) {
        this.atomicNumber = atomicNumber;
        this.symbol = symbol;
        this.mass = mass;
        this.covalent = covalent;
        this.vander = vander;
        this.group = group;
        this.c6 = c6;
    }

    protected String getSymbol() {
        return symbol;
    }

    protected Double getMass() {
        return mass;
    }

    protected Double getCovalent() {
        return covalent;
    }

    protected Double getVander() {
        return vander;
    }

    protected Integer getAtomicNumber() {
        return atomicNumber;
    }

    protected String getGroup() {
        return group;
    }


    public static AtomDetails[] buildAtomsDefaults () {
        AtomDetails[] datos = new AtomDetails[112];

        datos[0 ] = new AtomDetails(0 , "X", 0.0, 0.0, 0.0, "NONE", 0.0);
        datos[1 ] = new AtomDetails(1 , "H", 1.00794, 0.32, 1.2, "IA", 2.845);
        datos[2 ] = new AtomDetails(2 , "He", 4.002602, 0.93, 1.4, "VIIIA", 1.109);
        datos[3 ] = new AtomDetails(3 , "Li", 6.941, 1.23, 1.82, "IA", 0.787);
        datos[4 ] = new AtomDetails(4 , "Be", 9.012182, 0.9, 2.0, "IIA", 5.278);
        datos[5 ] = new AtomDetails(5 , "B", 10.811, 0.82, 2.0, "IIIA", 121.048);
        datos[6 ] = new AtomDetails(6 , "C", 12.011, 0.77, 1.7, "IVA", 26.36);
        datos[7 ] = new AtomDetails(7 , "N", 14.00674, 0.75, 1.55, "VA", 19.48);
        datos[8 ] = new AtomDetails(8 , "O", 15.9994, 0.73, 1.52, "VIA", 12.415);
        datos[9 ] = new AtomDetails(9 , "F", 18.9984, 0.72, 1.47, "VIIA", 10.518);
        datos[10] = new AtomDetails(10, "Ne", 20.1797, 0.71, 1.54, "VIIIA", 7.092);
        datos[11] = new AtomDetails(11, "Na", 22.989768, 1.54, 2.27, "IA", 3.068);
        datos[12] = new AtomDetails(12, "Mg", 24.305, 1.36, 1.73, "IIA", 12.247);
        datos[13] = new AtomDetails(13, "Al", 26.981539, 1.18, 2.0, "IIIA", 607.85);
        datos[14] = new AtomDetails(14, "Si", 28.0855, 1.11, 2.1, "IVA", 366.281);
        datos[15] = new AtomDetails(15, "P", 30.973762, 1.06, 1.8, "VA", 225.171);
        datos[16] = new AtomDetails(16, "S", 32.066, 1.02, 1.8, "VIA", 171.641);
        datos[17] = new AtomDetails(17, "Cl", 35.4527, 0.99, 1.75, "VIIA", 124.577);
        datos[18] = new AtomDetails(18, "Ar", 39.948, 0.98, 1.88, "VIIIA", 89.929);
        datos[19] = new AtomDetails(19, "K", 39.0983, 2.03, 2.75, "IA", 15.588);
        datos[20] = new AtomDetails(20, "Ca", 40.078, 1.74, 2.0, "IIA", 53.271);
        datos[21] = new AtomDetails(21, "Sc", 44.95591, 1.44, 2.0, "IIIB", 3.529);
        datos[22] = new AtomDetails(22, "Ti", 47.88, 1.32, 2.0, "IVB", 2.528);
        datos[23] = new AtomDetails(23, "V", 50.9415, 1.22, 2.0, "VB", 2.243);
        datos[24] = new AtomDetails(24, "Cr", 51.9961, 1.18, 2.0, "VIB", 1.662);
        datos[25] = new AtomDetails(25, "Mn", 54.93805, 1.17, 2.0, "VIIB", 1.272);
        datos[26] = new AtomDetails(26, "Fe", 55.847, 1.17, 2.0, "VIIIB", 1.151);
        datos[27] = new AtomDetails(27, "Co", 58.9332, 1.16, 2.0, "VIIIB", 1.14);
        datos[28] = new AtomDetails(28, "Ni", 58.6934, 1.15, 1.63, "VIIIB", 1.128);
        datos[29] = new AtomDetails(29, "Cu", 63.546, 1.17, 1.4, "IB", 1.323);
        datos[30] = new AtomDetails(30, "Zn", 65.39, 1.25, null, "IIB", 8.008);
        datos[31] = new AtomDetails(31, "Ga", 69.723, 1.26, 1.87, "IIIA", 427.057);
        datos[32] = new AtomDetails(32, "Ge", 72.61, 1.22, 2.0, "IVA", 338.151);
        datos[33] = new AtomDetails(33, "As", 74.92159, 1.2, 1.85, "VA", 256.927);
        datos[34] = new AtomDetails(34, "Se", 78.96, 1.16, 1.9, "VIA", 233.506);
        datos[35] = new AtomDetails(35, "Br", 79.904, 1.14, 1.85, "VIIA", 196.854);
        datos[36] = new AtomDetails(36, "Kr", 83.8, 1.12, 2.02, "VIIIA", 161.014);
        datos[37] = new AtomDetails(37, "Rb", 85.4678, 2.16, 2.0, "IA", 28.148);
        datos[38] = new AtomDetails(38, "Sr", 87.62, 1.91, 2.0, "IIA", 79.47);
        datos[39] = new AtomDetails(39, "Y", 88.90585, 1.62, 2.0, "IIIB", 14.639);
        datos[40] = new AtomDetails(40, "Zr", 91.224, 1.45, 2.0, "IVB", 9.309);
        datos[41] = new AtomDetails(41, "Nb", 92.90638, 1.34, 2.0, "VB", 8.608);
        datos[42] = new AtomDetails(42, "Mo", 95.94, 1.3, 2.0, "VIB", 6.569);
        datos[43] = new AtomDetails(43, "Tc", 98.0, 1.27, 2.0, "VIIB", 5.059);
        datos[44] = new AtomDetails(44, "Ru", 101.07, 1.25, 2.0, "VIIIB", 5.5);
        datos[45] = new AtomDetails(45, "Rh", 102.9055, 1.25, 2.0, "VIIIB", 4.857);
        datos[46] = new AtomDetails(46, "Pd", 106.42, 1.28, 1.63, "VIIIB", 4.136);
        datos[47] = new AtomDetails(47, "Ag", 107.8682, 1.34, 1.72, "IB", 5.085);
        datos[48] = new AtomDetails(48, "Cd", 112.411, 1.48, 1.58, "IIB", 17.66);
        datos[49] = new AtomDetails(49, "In", 114.82, 1.44, 1.93, "IIIA", 687.064);
        datos[50] = new AtomDetails(50, "Sn", 118.71, 1.41, 2.17, "IVA", 590.699);
        datos[51] = new AtomDetails(51, "Sb", 121.757, 1.4, 2.0, "VA", 485.947);
        datos[52] = new AtomDetails(52, "Te", 127.6, 1.36, 2.06, "VIA", 460.826);
        datos[53] = new AtomDetails(53, "I", 126.90447, 1.33, 1.98, "VIIA", 408.586);
        datos[54] = new AtomDetails(54, "Xe", 131.29, 1.31, 2.16, "VIIIA", 351.585);
        datos[55] = new AtomDetails(55, "Cs", 132.90543, 2.35, 2.0, "IA", 55.478);
        datos[56] = new AtomDetails(56, "Ba", 137.327, 1.98, 2.0, "IIA", 136.217);
        datos[57] = new AtomDetails(57, "La", 138.9055, 1.69, 2.0, "IIIB", 4.71);
        datos[58] = new AtomDetails(58, "Ce", 140.115, 1.65, 2.0, "LANS", 3.815);
        datos[59] = new AtomDetails(59, "Pr", 140.90765, 1.65, 2.0, "LANS", 3.191);
        datos[60] = new AtomDetails(60, "Nd", 144.24, 1.64, 2.0, "LANS", 3.03);
        datos[61] = new AtomDetails(61, "Pm", 145.0, 1.63, 2.0, "LANS", 2.61);
        datos[62] = new AtomDetails(62, "Sm", 150.36, 1.62, 2.0, "LANS", 2.209);
        datos[63] = new AtomDetails(63, "Eu", 151.965, 1.85, 2.0, "LANS", 2.109);
        datos[64] = new AtomDetails(64, "Gd", 157.25, 1.61, 2.0, "LANS", 1.907);
        datos[65] = new AtomDetails(65, "Tb", 158.92534, 1.59, 2.0, "LANS", 1.716);
        datos[66] = new AtomDetails(66, "Dy", 162.5, 1.59, 2.0, "LANS", 1.649);
        datos[67] = new AtomDetails(67, "Ho", 164.93032, 1.58, 2.0, "LANS", 1.595);
        datos[68] = new AtomDetails(68, "Er", 167.26, 1.57, 2.0, "LANS", 1.545);
        datos[69] = new AtomDetails(69, "Tm", 168.93421, 1.56, 2.0, "LANS", 1.285);
        datos[70] = new AtomDetails(70, "Yb", 173.04, 1.74, 2.0, "LANS", 47.195);
        datos[71] = new AtomDetails(71, "Lu", 174.967, 1.56, 2.0, "LANS", 13.842);
        datos[72] = new AtomDetails(72, "Hf", 178.49, 1.44, 2.0, "IVB", 10.036);
        datos[73] = new AtomDetails(73, "Ta", 180.9479, 1.34, 2.0, "VB", 11.93);
        datos[74] = new AtomDetails(74, "W", 183.85, 1.3, 2.0, "VIB", 8.126);
        datos[75] = new AtomDetails(75, "Re", 186.207, 1.28, 2.0, "VIIB", 6.365);
        datos[76] = new AtomDetails(76, "Os", 190.2, 1.26, 2.0, "VIIIB", 4.954);
        datos[77] = new AtomDetails(77, "Ir", 192.22, 1.27, 2.0, "VIIIB", 5.56);
        datos[78] = new AtomDetails(78, "Pt", 195.08, 1.3, 1.72, "VIIIB", 5.066);
        datos[79] = new AtomDetails(79, "Au", 196.96654, 1.34, 1.66, "IB", 7.218);
        datos[80] = new AtomDetails(80, "Hg", 200.59, 1.49, 1.55, "IIB", 21.891);
        datos[81] = new AtomDetails(81, "Tl", 204.3833, 1.48, 1.96, "IIIA", 665.972);
        datos[82] = new AtomDetails(82, "Pb", 207.2, 1.47, 2.02, "IVA", 605.78);
        datos[83] = new AtomDetails(83, "Bi", 208.98037, 1.46, 2.0, "VA", 523.633);
        datos[84] = new AtomDetails(84, "Po", 209.0, 1.46, 2.0, "VIA", 514.357);
        datos[85] = new AtomDetails(85, "At", 210.0, 1.45, 2.0, "VIIA", 473.466);
        datos[86] = new AtomDetails(86, "Rn", 222.0, 1.9, 2.0, "VIIIA", 421.345);
        datos[87] = new AtomDetails(87, "Fr", 223.0, null, null, "IA", 100.451);
        datos[88] = new AtomDetails(88, "Ra", 226.0, null, null, "IIA", 144.928);
        datos[89] = new AtomDetails(89, "Ac", 227.0, null, null, "IIIB", 8.478);
        datos[90] = new AtomDetails(90, "Th", 232.0381, 1.65, 2.0, "ACTS", 5.789);
        datos[91] = new AtomDetails(91, "Pa", 231.03588, null, null, "ACTS", 5.146);
        datos[92] = new AtomDetails(92, "U", 238.0289, 1.42, 1.86, "ACTS", 4.89);
        datos[93] = new AtomDetails(93, "Np", 237.0, 1.34, 2.0, "ACTS", 4.444);
        datos[94] = new AtomDetails(94, "Pu", 244.0, 1.55, 2.0, "ACTS", 3.742);
        datos[95] = new AtomDetails(95, "Am", 243.0, 1.89, 2.0, "ACTS", 3.035);
        datos[96] = new AtomDetails(96, "Cm", 247.0, 2.0, 2.0, "ACTS", 2.554);
        datos[97] = new AtomDetails(97, "Bk", 247.0, 2.0, 2.0, "ACTS", 2.615);
        datos[98] = new AtomDetails(98, "Cf", 251.0, 2.0, 2.0, "ACTS", 2.495);
        datos[99] = new AtomDetails(99, "Es", 252.0, 2.0, 2.0, "ACTS", 2.245);
        datos[100] = new AtomDetails(100, "Fm", 257.0, 2.0, 2.0, "ACTS", 2.193);
        datos[101] = new AtomDetails(101, "Md", 258.0, 2.0, 2.0, "ACTS", 1.966);
        datos[102] = new AtomDetails(102, "No", 259.0, 2.0, 2.0, "ACTS", 1.875);
        datos[103] = new AtomDetails(103, "Lr", 262.0, 2.0, 2.0, "ACTS", 1.833);
        datos[104] = new AtomDetails(104, "Rf", 267.0, 2.0, 2.0, "IVB", 1.833);
        datos[105] = new AtomDetails(105, "Db", 268.0, 2.0, 2.0, "VB", 1.833);
        datos[106] = new AtomDetails(106, "Sg", 271.0, 2.0, 2.0, "VIB", 1.833);
        datos[107] = new AtomDetails(107, "Bh", 270.0, 2.0, 2.0, "VIIB", 1.833);
        datos[108] = new AtomDetails(108, "Hs", 277.0, 2.0, 2.0, "VIIIB", 1.833);
        datos[109] = new AtomDetails(109, "Mt", 276.0, 2.0, 2.0, "VIIIB", 1.833);
        datos[110] = new AtomDetails(110, "Ds", 281.0, 2.0, 2.0, "VIIIB", 1.833);
        datos[111] = new AtomDetails(111, "Rg", 280.0, 2.0, 2.0, "IB", 1.833);

        return datos;
    }


    public static Object getBySymbol (AtomDetails[] source, String symbol, int fieldToSearch) {
        Object data = null;
        for(AtomDetails record : source) {
            if(record.symbol.equals(symbol)) {
                switch (fieldToSearch) {
                    case AtomDetails.ATMNUMBER:
                        data = record.atomicNumber;
                        break;
                    case AtomDetails.MASS:
                        data = record.mass;
                        break;
                    case AtomDetails.COVALENT:
                        data = record.covalent;
                        break;
                    case AtomDetails.VANDER:
                        data = record.vander;
                        break;
                    case AtomDetails.GROUP:
                        data = record.group;
                        break;
                    default:
                        data = record;
                        break;
                }
                break;
            }
        }

        return data;
    }

}
