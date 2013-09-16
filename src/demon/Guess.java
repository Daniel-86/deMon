package demon;

/**
 * User: daniel
 * Date: 14/06/13
 * Time: 02:09 AM
 */
public class Guess {

    Atomo[] atoms;


    public Guess(Atomo[] atoms) {
        this.atoms = atoms;
    }

    public void doTightBinding () {
        tightBindingInit(Globals.N_ATOM);

        int naos = 0;
        Utils.fillVectorWith(Globals.TB_AO_CONTRACT_COEFF, 0.0);
        for(int itbgrp = 0; itbgrp < Globals.N_TB_GRPS; itbgrp++) {
            int iatom = Globals.TB_GRP_LIST[Globals.TB_GRP_LOWER[itbgrp]];
        }

    }

    public void tightBindingInit(int nAtomsTBRegion) {
        Globals.N_TB_ATOMS = nAtomsTBRegion;
        tbList();
        generateTBOrbitalsPointers();
        generateTBOrbitalSymbols();
    }

    private void generateTBOrbitalSymbols() {
        /**
         * construye simbolos orbitales esfericos
         */
        for(int itbo = 0; itbo < Globals.N_TB_ORBITALS; itbo++) {
            int n = Globals.TB_ORBITAL_POINTER[itbo][3];
            int l = Globals.TB_ORBITAL_POINTER[itbo][4];
            int m = Globals.TB_ORBITAL_POINTER[itbo][5];

            // todo llenar el arreglo TB_ORBITALS_SYMBOLS
        }
    }

    public String printNOrbitals() {
        return "NUMBER OF TB ORBITALS:" + Globals.N_TB_ORBITALS;
    }

    /**
     * gentbop.f
     */
    private void generateTBOrbitalsPointers() {
        int lao = 0;

        Globals.N_TB_ORBITALS = 0;
        Globals.N_TB_SHELLS = 0;
        for(int i = 0; i < Globals.N_TB_ATOMS; i++) {
            Globals.TB_ORBITALS_LOWER[i] = Globals.N_TB_ORBITALS + 1;
            String confStr = generateElementConfigString();
            modifyConfigurationString(i, confStr);
            String[] tokens = confStr.split(" ");
            for(int ishl = 0; ishl < tokens.length; ishl++) {
                int nao = Integer.parseInt(tokens[i]);
                switch (tokens[i]) {
                    case "s":
                        lao = 0;
                        break;
                    case "p":
                        lao = 1;
                        break;
                    case "d":
                        lao = 2;
                        break;
                    case "f":
                        lao = 3;
                        break;
                }
                Globals.N_TB_SHELLS++;
                for(int mao = -lao; mao < lao; mao++) {
                    Globals.N_TB_ORBITALS++;
                    Globals.TB_ORBITAL_POINTER[Globals.N_TB_ORBITALS][1] = i;
                    Globals.TB_ORBITAL_POINTER[Globals.N_TB_ORBITALS][1] = Globals.N_TB_SHELLS;
                    Globals.TB_ORBITAL_POINTER[Globals.N_TB_ORBITALS][1] = nao;
                    Globals.TB_ORBITAL_POINTER[Globals.N_TB_ORBITALS][1] = lao;
                    Globals.TB_ORBITAL_POINTER[Globals.N_TB_ORBITALS][1] = mao;
                }
            }
            Globals.TB_ORBITALS_UPPER[i] = Globals.N_TB_ORBITALS;
        }

        /**
         * construye TB atom pointer
         */
        for(int iatom = 0; iatom < Globals.N_TB_ATOMS; iatom++) {
            Globals.TB_ATOM_POINTER[iatom][0] = atoms[iatom].carga.intValue();
            Globals.TB_ATOM_POINTER[iatom][1] = 0;
            for(int itbo = Globals.TB_ORBITALS_LOWER[iatom]; itbo < Globals.TB_ORBITALS_UPPER[iatom]; itbo++)
                Globals.TB_ATOM_POINTER[iatom][1] = StrictMath.max(Globals.TB_ATOM_POINTER[iatom][2], Globals.TB_ORBITAL_POINTER[itbo][2]);
            Globals.TB_ATOM_POINTER[iatom][2] = 0;

            String confStr = generateElementConfigString();
            modifyConfigurationString(iatom, confStr);
            String[] tokens = confStr.split(" ");
            for(int ishl = 0; ishl < tokens.length; ishl++) {
                int nElectron = Integer.parseInt(tokens[0]);
                Globals.TB_ATOM_POINTER[iatom][3] += nElectron;
            }
        }
    }

    private void modifyConfigurationString(int i, String confStr) {

    }

    private String generateElementConfigString() {
        return "ESTO DESPUES";
    }

    /**
     * genera una lista TB para identificar atomos equivalentes en TB
     */
    public void tbList() {
        Globals.N_TB_GRPS = Globals.N_TB_ATOMS;

        for(int i = 0; i < Globals.N_TB_GRPS; i++) {
            Globals.TB_OF[i] = i;
            Globals.TB_GRP_LOWER[i] = i;
            Globals.N_TB_EQU_ATOMS[i] = 1;
            Globals.TB_GRP_LIST[i] = i;
            Globals.TB_GRP_UPPER[i] = i;
        }

        int atomN = 0;
        int tbGroup = 0;

        for(int i = 0; i < Globals.N_TB_ATOMS; i++) {
            if( StrictMath.abs(Globals.ain[atomN] - Globals.ain[i]) < Globals.TINIEST && !isInList(Globals.TB_GRP_LIST, i)) {
                Globals.TB_GRP_UPPER[tbGroup]++;
                Globals.TB_GRP_LIST[Globals.TB_GRP_UPPER[tbGroup]] = i;
                Globals.N_TB_EQU_ATOMS[tbGroup]++;
            }
        }

        if(Globals.TB_GRP_UPPER[Globals.N_TB_GRPS] != Globals.N_TB_ATOMS) {
            for(int i = 0; i < Globals.N_TB_ATOMS; i++) {
                if(!isInList(Globals.TB_GRP_LIST, i)) {
                    atomN = i;
                    tbGroup++;
                    Globals.TB_GRP_LOWER[tbGroup] = Globals.TB_GRP_UPPER[tbGroup-1] + 1;
                    Globals.TB_GRP_LIST[Globals.TB_GRP_LOWER[tbGroup]] = i;
                    Globals.TB_GRP_UPPER[tbGroup] = Globals.TB_GRP_LOWER[tbGroup];
                }
            }
        }

        for(int i = 0; i < tbGroup; i++) {
            for(int j = Globals.TB_GRP_LOWER[i]; j < Globals.TB_GRP_UPPER[i]; i++)
                Globals.TB_OF[Globals.TB_GRP_LIST[j]] = i;
        }
    }

    public boolean isInList(Integer[] list, int val) {
        for (Integer aList : list)
            if (val == aList)
                return true;
        return false;
    }


    public void initAtomicSCF(int nAtom) {
        Globals.molecularCharge = 0;
        Globals.N_ATOM = 1;
        atoms[0].atomicN = atoms[nAtom].atomicN;
        Globals.nElectrons = Globals.TB_ATOM_POINTER[nAtom][3];
        atoms[0].carga = Double.valueOf(Globals.TB_ATOM_POINTER[nAtom][1]);
        Globals.spinMultiplicity = (Globals.nElectrons % 2 ==0)? 2: 1;

        /**
         * define banderas para calculos atomicos esfericos
         */
//        Globals.oks = true;
//        Globals.roks = true;
//        Globals.cdxc = true;
//        if(Globals.lap) {
//            Globals.lap = false;
//            Globals.gga = true;
//            Globals.styp = "PBE96";
//            Globals.ctyp = "PBE";
//        }
//        Globals.diis = true;
//        Globals.scf = true;
//        Globals.sphorb = true;
//        Globals.cfgatom = true;
//        Globals.setaocc = true;
//        Globals.symmetry = true;
//
//        Globals.embed = false;
//        Globals.fixcfg = false;
//        Globals.fixlast = false;
//        Globals.lshift = false;
//        Globals.moex = false;
//        Globals.multipole = false;
//        Globals.restart = false;
//        Globals.rotgrd = false;
//        Globals.smear = false;
//
//        Globals.cdmix = 0.2;
//        Globals.cfgcvc = 50;
//        Globals.fixtol = 0;
//        Globals.grdtol = 1e-6;
//        Globals.scfmax = 100;
//        Globals.scftol = 1e-6;
//        Globals.usesvd = true;
//
//        Globals.cdftyp = "MIXING";
//        Globals.erityp = "NOT DIRECT";
//        Globals.fit_typ = "ALL";
//        Globals.matdia = "JACOBI";
//
//        /**
//         * establece campos externos para cuanto m number splitting
//         */
//        Globals.embed = true;
//        Globals.ccset = 2;
//        Globals.nembed = 1;
//        Globals.qembed[0] = 0.1;
//        atoms[0].coords.set(0, 0d);
//        atoms[0].coords.set(1, 0d);
//        atoms[0].coords.set(2, 0d);
//        Globals.cembed[0] = 0d;
//        Globals.cembed[1] = 0d;
//        Globals.cembed[2] = 100d;

        setThresholds();
        calculateSphericalAtomConfig(nAtom, Globals.aoccna, Globals.aoccnb, Globals.elcfga, Globals.elcfgb);
        loadBasis(nAtom);
        orbital();

    }

    private void orbital() {
        generateCartesianOrbitalPtrs();
        generateCartesianOrbitalSymbols();
        generateSphericOrbitalPtrs();
        generateSphericOrbitalSymbols();

        /**
         * define el numero de orbitales procesados en el calculo SCF
         */
        if(Globals.SPHORB) {
            Globals.N_ORB = Globals.N_SPH;
            Globals.N_ORTH = Globals.N_SPH;
        } else {
            Globals.N_ORB = Globals.N_ORBITALS;
            Globals.N_ORTH = Globals.N_ORBITALS;
        }

        /**
         * establece SCF MO printing pointers
         */
        if(Globals.lowerLimitMOS == 0 && Globals.upperLimitMOS == 0) {
            Globals.lowerLimitMOS = 1;
            Globals.upperLimitMOS = Globals.N_ORB;
        }

        /**
         * establece SCF MO printing pointers
         */
        if(Globals.lowerLimitDOS == 0 && Globals.upperLimitDOS == 0) {
            Globals.lowerLimitDOS = 1;
            Globals.upperLimitDOS = Globals.N_ORB;
        }

        /**
         * imprime la tabla GTO
         */

        /**
         * normaliza los GTO
         */
        normalizeGTO();

        /**
         * normaliza los STO
         */
        normalizeSTO();

        /**
         * imprime la tabla STO
         */
    }

    private void normalizeSTO() {
        //To change body of created methods use File | Settings | File Templates.
    }

    private void normalizeGTO() {
        double norm;
        int ishl, lshl, igto;

        for(ishl = 0; ishl < Globals.N_SHELLS; ishl++) {
            lshl = Globals.sphericPointer[ishl][2];
            norm = StrictMath.pow(2, lshl) * StrictMath.pow(2/StrictMath.PI, 0.75);
            for(igto = Globals.lowerLimitGTO[ishl]; igto < Globals.upperLimitGTO[ishl]; igto++)
                Globals.gaussianContractionCoeficient[igto] =
                        norm
                                * StrictMath.pow(StrictMath.pow(Globals.gaussianExponent[igto], 2*lshl+3), 0.025)
                                * Globals.gaussianContractionCoeficient[igto];
        }
    }

    private void generateSphericOrbitalPtrs() {
        Utils.fillVectorWith(Globals.lowerlimitSPH, 0);
        Utils.fillVectorWith(Globals.upperlimitSPH, 0);
        Utils.fillMatrixWith(Globals.sphericPointer, 0);

        int origin, m, l, n;

        /**
         * inicializa los contadores de orbital esferico
         */
        Globals.N_SPH = 0;

        /**
         * construye el apuntador orbital esferico
         */
        for(int ishl = 0; ishl < Globals.N_SHELLS; ishl++) {
            Globals.lowerlimitSPH[ishl] = Globals.N_SPH + 1;
            origin = Globals.sphericPointer[ishl][0];
            n = Globals.sphericPointer[ishl][1];
            l = Globals.sphericPointer[ishl][2];
            for(m = -l; m < l; m++) {
                Globals.N_SPH++;
                Globals.sphericPointer[Globals.N_SPH][0] = origin;
                Globals.sphericPointer[Globals.N_SPH][1] = ishl;
                Globals.sphericPointer[Globals.N_SPH][2] = n;
                Globals.sphericPointer[Globals.N_SPH][3] = l;
                Globals.sphericPointer[Globals.N_SPH][4] = m;
            }
            Globals.upperlimitSPH[ishl] = Globals.N_SPH;
        }
    }

    private void generateSphericOrbitalSymbols() {
        //To change body of created methods use File | Settings | File Templates.
    }

    private void generateCartesianOrbitalSymbols() {
        //To change body of created methods use File | Settings | File Templates.
    }

    private void generateCartesianOrbitalPtrs() {
        int lsto;
        int ishl;

        Utils.fillVectorWith(Globals.ll, 0);
        Utils.fillVectorWith(Globals.ul, 0);
        Utils.fillVectorWith(Globals.lowerLimitGTO, 0);
        Utils.fillVectorWith(Globals.upperLimitGTO, 0);
        Utils.fillVectorWith(Globals.lowerLimitSTO, 0);
        Utils.fillVectorWith(Globals.upperLimitSTO, 0);
        Utils.fillVectorWith(Globals.lowerLimitShell, 0);
        Utils.fillVectorWith(Globals.upperLimitShell, 0);
        Utils.fillMatrixWith(Globals.STO_PTR, 0);


        /**
         * inicializa los contadores
         */
        int isto = 0;
        Integer[] number = new Integer[Globals.MAX_SHELL];

        /**
         * construye apuntadores STO
         */
        for(ishl = 0; ishl < Globals.N_SHELLS; ishl++) {
            lsto = Globals.shellPointer[ishl][2];
            for(int ax = lsto; ax > 0; ax--)
                for(int ay = lsto-ax; ay > 0; ay--) {
                    int az = lsto-ax-ay;
                    isto++;
                    Globals.STO_PTR[ishl][0] = Globals.shellPointer[ishl][0];
                    Globals.STO_PTR[ishl][1] = ishl;
                    Globals.STO_PTR[ishl][2] = ax;
                    Globals.STO_PTR[ishl][3] = ay;
                    Globals.STO_PTR[ishl][4] = az;
                }
        }

        if(isto != Globals.N_SHELLS) {
            System.out.println("GENOP','ERROR WHILE BUILDING ORBITAL POINTER");
        }

        /**
         * genera los orbitales limite para los atomos
         */
        int iatom = 1;
        number[0] = 0;
        for(isto = 0; isto < Globals.N_ORBITALS; isto++) {
            if(iatom == Globals.STO_PTR[isto][0])
                number[iatom]++;
            else {
                iatom = Globals.STO_PTR[isto][0];
                number[iatom] = 1;
            }
        }
        Globals.ll[0] = 1;
        Globals.ul[0] = number[0];
        for(iatom = 1; iatom < Globals.N_ATOM; iatom++) {
            Globals.ll[iatom] =Globals.ul[iatom+1] + 1;
            Globals.ul[iatom] = Globals.ll[iatom] + number[iatom] - 1;
        }

        /**
         * genera los shell limite para los atomos
         */
        for(iatom = 0; iatom < Globals.N_ATOM; iatom++) {
            number[iatom] = 0;
            ishl = 0;
            for(isto = Globals.ll[iatom]; isto < Globals.ul[iatom]; isto++) {
                if(ishl != Globals.STO_PTR[isto][1]) {
                    number[iatom]++;
                    ishl = Globals.STO_PTR[isto][1];
                }
            }
        }
        Globals.ll[0] = 1;
        Globals.ul[0] = number[0];
        for(iatom = 1; iatom < Globals.N_ATOM; iatom++) {
            Globals.lowerLimitShell[iatom] =Globals.upperLimitShell[iatom+1] + 1;
            Globals.upperLimitShell[iatom] = Globals.lowerLimitShell[iatom] + number[iatom] - 1;
        }

        /**
         * genera orbital limite para los shells
         */
        ishl = 1;
        number[0] = 0;
        for(isto = 0; isto < Globals.N_ORBITALS; isto++) {
            if(ishl == Globals.STO_PTR[isto][1])
                number[ishl]++;
            else {
                ishl = Globals.STO_PTR[isto][1];
                number[ishl] = 1;
            }
        }
        Globals.lowerLimitSTO[0] = 1;
        Globals.upperLimitSTO[0] = number[0];
        for(ishl = 1; ishl < Globals.N_SHELLS; ishl++) {
            Globals.lowerLimitSTO[ishl] =Globals.upperLimitSTO[ishl-1] + 1;
            Globals.upperLimitSTO[ishl] = Globals.lowerLimitSTO[ishl] + number[ishl] - 1;
        }

        /**
         * genera limites gausianos para los orbitales
         */
        Globals.lowerLimitGTO[0] = 1;
        Globals.upperLimitGTO[0] = Globals.shellPointer[0][3];
        for(ishl = 1; ishl < Globals.N_SHELLS; ishl++) {
            Globals.lowerLimitGTO[ishl] =Globals.upperLimitGTO[ishl-1] + 1;
            Globals.upperLimitGTO[ishl] = Globals.lowerLimitGTO[ishl] + Globals.shellPointer[0][3] - 1;
        }

        /**
         * encuentra el L mas grande para cada atomo y para el sistema
         */
        Utils.fillVectorWith(Globals.MAX_L_FOR_ORB_SHELL, 0);
        for(iatom = 0; iatom < Globals.N_ATOM; iatom++) {
            Globals.MAX_L_FOR_ORB_SHELL[iatom] = 0;
            for(ishl = Globals.lowerLimitShell[iatom]; ishl < Globals.upperLimitShell[iatom]; ishl++)
                Globals.MAX_L_FOR_ORB_SHELL[iatom] = StrictMath.max(Globals.MAX_L_FOR_ORB_SHELL[iatom], Globals.shellPointer[ishl][3]);
            Globals.MAX_L_FOR_ORB_SHELL[0] = StrictMath.max(Globals.MAX_L_FOR_ORB_SHELL[0], Globals.MAX_L_FOR_ORB_SHELL[iatom]);
        }
    }

    private void loadBasis(int nAtom) {

        int lshl;

        /**
         * inicializa contadores
         */
        int nMolGto = Globals.N_GAUSSIAN;
        double norm, renorm;

        Globals.N_SHELLS = 0;
        Globals.N_ORBITALS = 0;
        Globals.N_GAUSSIAN = 0;

        /**
         * carga la informacion del conjunto de bases
         */
        for(int ishl = Globals.lowerLimitShell[nAtom]; ishl < Globals.upperLimitShell[nAtom]; ishl++) {
            Globals.N_SHELLS++;
            Globals.shellPointer[Globals.N_SHELLS][1] = 1;
            Globals.shellPointer[Globals.N_SHELLS][2] = Globals.shellPointer[ishl][2];
            Globals.shellPointer[Globals.N_SHELLS][3] = Globals.shellPointer[ishl][3];
            Globals.shellPointer[Globals.N_SHELLS][4] = Globals.shellPointer[ishl][4];
            lshl = Globals.shellPointer[Globals.N_SHELLS][3];
            norm = StrictMath.pow(2, Globals.N_SHELLS) * StrictMath.pow(2/StrictMath.PI, 0.75);
            Globals.N_ORBITALS += Globals.N_Cartesian_Orbitals[Globals.shellPointer[Globals.N_SHELLS][3]];
            for(int igto = 0; igto < Globals.shellPointer[Globals.N_SHELLS][4]; igto++) {
                Globals.gaussianExponent[Globals.N_GAUSSIAN + igto] = Globals.gaussianExponent[Globals.lowerLimitGTO[ishl]-1+igto];
                renorm = norm * StrictMath.pow(StrictMath.pow(Globals.gaussianExponent[Globals.N_GAUSSIAN], 2*lshl+3), 0.025);
                Globals.gaussianContractionCoeficient[Globals.N_GAUSSIAN + igto] = Globals.gaussianContractionCoeficient[Globals.lowerLimitGTO[ishl]-1+igto]/renorm;
            }
            Globals.N_GAUSSIAN += Globals.shellPointer[Globals.N_SHELLS][4];
        }

        /**
         * todo limpiar los 'campos', destruir la informacion de la molecula
         * SHLPTR(NSHL+1:MAXSHL,:) = 0
         GCC(NGTO+1:MAXGTO) = 0.0
         */
    }

    private void calculateSphericalAtomConfig(int nAtom, Double[][] aoccna, Double[][] aoccnb, Integer[] elcfga, Integer[] elcfgb) {

        Integer [] nao =  new Integer[3];
        int lao = 0;
        double full = 0;

        Utils.fillVectorWith(elcfga, 0);
        Utils.fillVectorWith(elcfgb, 0);
        Utils.fillMatrixWith(aoccna, 0d);
        Utils.fillMatrixWith(aoccnb, 0d);
        Utils.fillVectorWith(nao, 0);

        /**
         * obtener la configuracion electronica del elemento
         */
        String confStr = generateElementConfigString();
        modifyConfigurationString(nAtom, confStr);
        String[] tokens = confStr.split(" ");
//        for(int ishl = 0; ishl < tokens.length; ishl++) {
//            Globals.nElectrons = Integer.parseInt(tokens[i]);
//            switch (tokens[i]) {
//                case "s":
//                    lao = 0;
//                    full = 2;
//                    break;
//                case "p":
//                    lao = 1;
//                    full = 6;
//                    break;
//                case "d":
//                    lao = 2;
//                    full = 10;
//                    break;
//                case "f":
//                    lao = 3;
//                    full = 14;
//                    break;
//            }
//            nao[lao]++;
//            for(int mao = -lao; mao < lao; mao++) {
//                Globals.elcfga[lao*lao + lao + mao]++;
//                Globals.elcfgb[lao*lao + lao + mao]++;
//                Globals.aoccna[lao*lao + lao+mao][nao[lao]] = Globals.nElectrons/full;
//                Globals.aoccnb[lao*lao + lao+mao][nao[lao]] = Globals.nElectrons/full;
//            }
//        }
    }

    private void setThresholds() {
//        if(Globals.erityp.equals("DIRECT")) {
//            Globals.erithresh = Globals.eritol/Globals.nElectrons;
//            Globals.rhothresh = Globals.eritol/Globals.nElectrons;
//            Globals.cthresh = StrictMath.min((Globals.eritol/Globals.nElectrons), 1e-12);
//            Globals.pthresh = 1e-5 * Globals.grdtol;
//        } else {
//            Globals.erithresh = StrictMath.min((Globals.eritol/Globals.nElectrons), 1e-10);
//            Globals.rhothresh = StrictMath.min((Globals.eritol/Globals.nElectrons), 1e-10);
//            Globals.cthresh = StrictMath.min((Globals.eritol/Globals.nElectrons), 1e-12);
//            Globals.pthresh = 1e-5 * Globals.grdtol;
//        }
    }

}
