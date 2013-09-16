package demon;

import java.text.DecimalFormat;

/**
 * User: daniel
 * Date: 15/04/13
 * Time: 06:49 PM
 */
public class Symmetry {
    double[][] rotMat = new double[3][3];
    int maxcn = 20; // Maximum rotational degree of a C axis.
    int maxsec = maxcn + 1; // Maximum number of C axes.
    int maxsn = 2*maxcn; // Maximum rotational degree of a S axis.
    int maxses = 2*maxcn + 1; // Maximum number of S axes.
    int maxsig = maxcn + 1; // Maximum number of mirror planes.
    double[][][] sec = new double[maxsec][maxcn-2][3]; // List of C axes (SEC(I,J) => Ith J-fold C axis).
    double[][][] ses = new double[maxses][maxsn-2][3]; // List of S axes (SES(I,J) => Ith J-fold S axis).
    double[][] sig = new double[maxsig][3]; // list of mirror planes
    int ncn = 1; // element counter C.
    int nsn = 1; // element counter S.
    int[] nsec = new int[maxcn-1];
    int[] nses = new int[maxsn-1];
    int nsig = 0;
    double molecularMass = 0;
    double[][][] coord = new double[3][Globals.max_atom][7]; // coordenada, numero atomo, tipo coordenada(cartesiana...)
    double[] shift = new double[3];
    double[][] tenmat = new double[3][3];
    double eps = 1e-10;
    private Double[][] rotationMatrix;
    MineVector [][] cAxisList;
    MineVector [][] sAxisList;
    MineVector [] mirrorPlaneList;

    Double[] tenval;
    Double[][] tenVec;

    public void doMolSym(Atomo[] atomos) {
        for(Atomo a: atomos)
            molecularMass += a.masa;
        for(int i = 0; i < 3; i++) {
            for(Atomo a: atomos) {
                shift[i] += a.masa*a.coords.get(i);
            }
            shift[i] /= molecularMass;
        }
        for(int i = 0; i < 3; i++) {
            for(Atomo a: atomos) {
                a.coords.set(i, shift[i]);
//                a.coords.get(i) = shift[i];
            }
        }
        for(Atomo a: atomos) {
            tenmat[0][0] += a.masa*StrictMath.pow(a.coords.get(1), 2) + StrictMath.pow(a.coords.get(2), 2);
            tenmat[1][0] += a.masa*a.coords.get(1) + a.coords.get(0);
            tenmat[2][0] += a.masa*a.coords.get(2) + a.coords.get(0);
            tenmat[1][1] += a.masa*StrictMath.pow(a.coords.get(0), 2) + StrictMath.pow(a.coords.get(2), 2);
            tenmat[2][1] += a.masa*a.coords.get(2) + a.coords.get(1);
            tenmat[2][2] += a.masa*StrictMath.pow(a.coords.get(0), 2) + StrictMath.pow(a.coords.get(1), 2);
        }
        tenmat[0][0] += eps;
        tenmat[2][2] += eps;
    }

    /**
     *
     * @param matrix matriz a ser simetrizada.
     * @param dimension dimension de la matriz.
     * @param n numero de columnas o filas de la matriz.
     * @param option tipo de simetrizacion.
     */
    public void symmat(double[][] matrix, int dimension, int n, String option) {
        switch(option) {
            case "LOWUP":
                for(int i=0; i<n-1; i++)
                    dcopy(n-1, matrix[i+1][i], 1, matrix[i][i+1], dimension);
        }
    }

    public void dcopy(int n, double dx, int incx, double dy, int incy ) {
        if(n<+0)
            return;
        if(incx!=1 && incy!=1) {
            int ix = 1;
            int iy = 1;
            if(incx<0)
                ix = (-n+1) * incx + 1;
            if(incy<0)
                iy = (-n+1) * incy + 1;
            for(int i=0; i<n; i++) {

            }
        }
    }





    public void fillVectorWith(Object[] vector, Object val) {
        for(int i = 0; i < vector.length; i ++)
                vector[i] = val;
    }


    public void fillMatrixWith(Object[][] matrix, Object val) {
//        double[][] matrix
        int nRows = matrix.length;
        int nCols = matrix[0].length;
        for(int i = 0; i < nRows; i ++)
            for(int j = 0; j < nCols; j++)
                matrix[i][j] = val;
    }


    public void printVector(Object[] vector) {
        for (Object aVector : vector) System.out.print(aVector + " ");
        System.out.println();
    }


    public void printMatrix(Object[][] matrix) {
        int nRows = matrix.length;
        int nCols = matrix[0].length;
        for (Object[] aMatrix : matrix) {
            for (int j = 0; j < nCols; j++)
                System.out.print(aMatrix[j] + " ");
            System.out.println();
        }
    }

    public Object [][] buildIdentityMatrix(int dimension, int type) {
        Object [][] matrix;
        Object val;
        switch(type) {
            case 0:
                matrix = new Integer [dimension][dimension];
                fillMatrixWith(matrix, 0);
                val = 1;
                break;
            case 1:
                matrix = new Double [dimension][dimension];
                fillMatrixWith(matrix, 0.0);
                val = 1.0;
                break;
            default:
                return null;
        }
        for(int i = 0; i < dimension; i ++)
            matrix[i][i] = val;

        return matrix;
    }

    public void init() {
        /**
         * inicializa a 0 la matriz de rotación
         */
        rotationMatrix = new Double[3][3];
        fillMatrixWith(rotationMatrix, 0.0);
//        System.out.println(System.getProperty("line.separator") + "Matriz de rotacion");
//        printMatrix(rotationMatrix);



        /**
         * inicializa la lista de los elementos de simetria
         */
        cAxisList = new MineVector[Globals.MAX_C_AXIS_ROT_DEGREE - 1][Globals.MAX_C_AXIS];
        fillMatrixWith( cAxisList, new MineVector(3).fillWith(0.0) );
//        System.out.println(System.getProperty("line.separator") + "Matriz cAxisList");
//        printMatrix(cAxisList);

        sAxisList = new MineVector[Globals.MAX_S_AXIS_ROT_DEGREE - 1][Globals.MAX_S_AXIS];
        fillMatrixWith( sAxisList, new MineVector(3).fillWith(0.0) );

        Globals.N_MIRROR_PLANES = 0;
        mirrorPlaneList = new MineVector[Globals.MAX_MIRROR_PLANES];
        fillVectorWith(mirrorPlaneList, new MineVector(3).fillWith(0.0));
//        System.out.println(System.getProperty("line.separator") + "MirrorPlaneList");
//        printVector(mirrorPlaneList);
//        for(int i = 2; i < Globals.MAX_C_AXIS_ROT_DEGREE; i++)
//            for(int j = 1; j < Globals.MAX_C_AXIS; j++) {
//                cAxisList[j][i][1] = 0.0;
//                cAxisList[j][i][2] = 0.0;
//                cAxisList[j][i][3] = 0.0;
//            }
//        for(int i = 2; i < Globals.MAX_S_AXIS_ROT_DEGREE; i++)
//            for(int j = 1; j < Globals.MAX_S_AXIS; j++) {
//                sAxisList[j][i][1] = 0.0;
//                sAxisList[j][i][2] = 0.0;
//                sAxisList[j][i][3] = 0.0;
//            }
//        for(int i = 1; i < Globals.MAX_MIRROR_PLANES; i++) {
//            mirrorPlaneList[i][1] = 0.0;
//            mirrorPlaneList[i][2] = 0.0;
//            mirrorPlaneList[i][3] = 0.0;
//        }

        /**
         * inicializa los contadores de los elementos de simetria
         */
        Globals.HIGEST_C_AXIS_DEGREE = 1;
        Globals.HIGEST_S_AXIS_DEGREE = 1;
        fillVectorWith(Globals.N_C_AXIS, 0);
        fillVectorWith(Globals.N_S_AXIS, 0);

        /**
         * inicializa las variables logicas
         */
        Globals.CUBIC = false;
        Globals.DGROUP = false;
        Globals.IGROUP = false;
        Globals.INVERS = false;
        Globals.LINEAR = false;
        Globals.MAXIS = false;
        Globals.OGROUP = false;
        Globals.SGROUP = false;
        Globals.SIGMAD = false;
        Globals.SIGMAV = false;
        Globals.SIGMAH = false;
        Globals.TGROUP = false;

        /**
         * inicializa lista de atomos equivalentes (simetria C1)
         */
//        Globals.N_GRP_EQ_ATOMS = Globals.N_ATOM;
        Globals.N_GRP_EQ_ATOMS = Globals.N_ATOM;
        for(int i = 0; i < Globals.N_GRP_EQ_ATOMS; i++) {
            Globals.equivalentAtomsGroup[i] = i;
            Globals.groupLowerLimit[i] = i;
            Globals.nEquivalentAtoms[i] = 1;
            Globals.equivalentAtomsGroupOrderedList[i] = i;
            Globals.groupUpperLimit[i] = i;
        }
        System.out.println("GRP UPP LIM " + Globals.groupUpperLimit.length);
        Globals.pointGroupOrder = 1;

        /**
         * cargar la orientacion de entrada
         */
//        for(int i = 0; i < Globals.N_ATOM; i++) {
//            Globals.atomicCoordinates[0][i][1] = Globals.atomicCoordinates[0][i][0];
//            Globals.atomicCoordinates[1][i][1] = Globals.atomicCoordinates[1][i][0];
//            Globals.atomicCoordinates[2][i][1] = Globals.atomicCoordinates[2][i][0];
//        }

        /**
         * inicializa la matriz de transformacion para: orientacion de entrada -> orientacion estandar
         */
        Globals.inputToStandarMatrix = (Double[][]) buildIdentityMatrix(3, 1);
//        System.out.println(System.getProperty("line.separator") + "Matriz inputToStandar");
//        printMatrix(Globals.inputToStandarMatrix);

        /**
         * establece la tolerancia para coordenadas atomicas equivalentes
         */
        Globals.cartesianCoordinatesTolerance = 5e-5;

        /**
         * generar los numeros de identificacion atomica. Los valores AIN son codigos de simetria para cada atomo
         * Solo los atomos con el mismo AIN pueden ser equivalentes en simetria.
         */
        for(int i = 0; i < Globals.N_ATOM; i++) {
            Globals.ain[i] = Globals.nuclearCharges[i] + Globals.atomicMass[i];
            System.out.println("antes for en symmetry init " + Globals.lowerLimitShell[i] + "  i:" + i);
            for(int j = Globals.lowerLimitShell[i]; j < Globals.upperLimitShell[i]; j++)
                for(int k = Globals.lowerLimitGTO[j]; k < Globals.upperLimitGTO[j]; k++)
                    Globals.ain[i] += Globals.gaussianContractionCoeficient[k] + Globals.gaussianExponent[k];
            for(int j = Globals.lowerAuxSet[i]; j < Globals.upperAuxSet[i]; j++)
                for(int k = Globals.lowerAuxShell[j]; k < Globals.upperAuxShell[j]; k++)
                    Globals.ain[i] += Globals.gaussianExponent[Globals.N_GAUSSIAN + j];
        }
    }


    public double determineMolecularMass(Atomo[] atoms) {
        double sum = 0.0;

        for(Atomo a: atoms) {
            sum += a.masa;
        }

        return sum;
    }

    public Double[] determineTranslationVector(Atomo[] atoms, double molecularMass) {
        Double [] vector = new Double[3];
        fillVectorWith(vector, 0.0);
        for(Atomo a: atoms)
            for(int i = 0; i < 3; i ++) {
                vector[i] += a.masa * a.coords.get(i) / molecularMass;
            }
        return vector;
    }

    /**
     * es tensor.f
     * determina la simetria molecular
     */
    public void findSymmetry(Atomo[] atoms) {
//        System.out.println("SYMETRY FIND");
        /**
         * suma las masas atomicas
         */
        Globals.molecularMass = determineMolecularMass(atoms);
//        Globals.molecularMass = sumVector(Globals.atomicMass);

        /**
         * encuentra el vector de traslacion del origen
         * (es un vector de 3 elementos, uno por cada coordenada. Cada elemento se define como la sumatoria de los productos
         * de la masa por la coordenada correspondiente al elemento a calcular.
         */
        fillVectorWith(Globals.translationVector, 0.0);
        Globals.translationVector = determineTranslationVector(atoms, Globals.molecularMass);
//        Globals.translationVector = Globals.molecularMass * horizontalVectorXMatrix(Globals.atomicMass, Globals.atomicCoordinates);
        for(Atomo a: atoms){
            a.shiftCenterMass(Globals.translationVector);
            a.printCoords();
        }

        /**
         * construir el tensor molecular de inercia
         */
        Double [][] molecularInertiaTensor = new Double[atoms.length][3];
        fillMatrixWith(molecularInertiaTensor, 0.0);
        for(Atomo a: atoms) {
            molecularInertiaTensor[0][0] += a.masa * (StrictMath.pow(a.coords.get(1), 2) + StrictMath.pow(a.coords.get(2), 2));
            molecularInertiaTensor[1][0] -= a.masa * a.coords.get(1) * a.coords.get(0);
            molecularInertiaTensor[2][0] -= a.masa * a.coords.get(2) * a.coords.get(0);
            molecularInertiaTensor[1][1] += a.masa * (StrictMath.pow(a.coords.get(0), 2) + StrictMath.pow(a.coords.get(2), 2));
            molecularInertiaTensor[2][1] -= a.masa * a.coords.get(2) * a.coords.get(1);
            molecularInertiaTensor[2][2] += a.masa * (StrictMath.pow(a.coords.get(0), 2) + StrictMath.pow(a.coords.get(1), 2));
        }

        /**
         * desplazar los elementos diagonales
         */
        double EPS = 1e-10;
        molecularInertiaTensor[0][0] += EPS;
        molecularInertiaTensor[2][2] -= EPS;

        System.out.println(System.getProperty("line.separator") + "Molecular Inertia Tensor");
        printMatrix(molecularInertiaTensor);

        /**
         * simetrizar el tensor de inercia molecular
         */
        simetrize(molecularInertiaTensor, 3, 3, "lowup");
        System.out.println(System.getProperty("line.separator") + "Molecular Inertia Tensor SIMETRICA");
        printMatrix(molecularInertiaTensor);

        /**
         * imprimir los momentos del tensor de inercia molecular
         * solo si esta activa esa opcion
         */


        /**
         * diagonaliza el tensor de inercia molecular
         */
//        Double[][] tenVec = molecularInertiaTensor;
//        Jacobi jacobi = new Jacobi();
//        jacobi.SolJacobi();
        System.out.println(System.getProperty("line.separator") + "JACOBI");
        tenVec = new Double[][]{{-0.61237202, 0.79056974, 0d}, {0d, 0d, 1d}, {0.79056974, 0.61237202, 0d}};
        printMatrix(tenVec);
        System.out.println(System.getProperty("line.separator") + "TenVal");
        tenval = new Double[]{0.63222680, 1.15408902, 1.78631582};
        printVector(tenval);
//        diagonalizeMatrix(tenVec, 3, 3, "JACOBI", tenval);

        /**
         * se asegura que los ejes principales cumplan la regla de la mano derecha
         */
//        vectorProduct(tenVec[0][0], tenVec[1][0], tenVec[2][0],
//                tenVec[0][1], tenVec[1][1], tenVec[2][1],
//                tenVec[0][2], tenVec[1][2], tenVec[2][2]);
        /**
         * imprimir simetria si es requerido (por defecto no)
         */
//        if(isZeroVector(tenval))
//            System.out.println("tensor de inercia molecular invalido");
//        Globals.LINDEG = true;
    }

    /**
     *
     * @param matrix
     * @param dimension
     * @param nRowOrCol
     * @param option
     */
    public void simetrize(Double [][] matrix, int dimension, int nRowOrCol, String option) {
        switch(option) {
            case "uplow":
                for(int i = 0; i < dimension; i++)
                    for(int j = i+1; j < dimension; j++) {
                        if(i == j)
                            continue;
                        matrix[j][i] = matrix[i][j];
                    }
//                for(int i = 0; i < nRowOrCol-1; i++) {
//                    dcopy(nRowOrCol - i, matrix[i+1][i], 1, matrix[i][i+1], dimension);
//                }
                break;
            case "lowup":
                for(int i = 0; i < dimension; i++)
                    for(int j = i+1; j < dimension; j++) {
                        if(i == j)
                            continue;
                        matrix[i][j] = matrix[j][i];
                    }
//                for(int i = 0; i < nRowOrCol-1; i++) {
//                    dcopy(nRowOrCol - i, matrix[i][i+1], dimension, matrix[i+1][i+1], 1);
//                }
                break;
            case "-uplow":
                for(int i = 0; i < dimension; i++)
                    for(int j = i+1; j < dimension; j++) {
                        if(i == j)
                            continue;
                        matrix[j][i] = -matrix[i][j];
                    }
//                for(int i = 0; i < nRowOrCol-1; i++) {
//                    dcopy(nRowOrCol - i, matrix[i+1][i], 1, matrix[i][i+1], dimension);
//                }
//                for(int i = 1; i < nRowOrCol; i++) {
//                    dscal(-1, -1.0, matrix[0][i], 1);
//                }
                break;
            case "-lowup":
                for(int i = 0; i < dimension; i++)
                    for(int j = i+1; j < dimension; j++) {
                        if(i == j)
                            continue;
                        matrix[i][j] = -matrix[j][i];
                    }
//                for(int i = 0; i < nRowOrCol-1; i++) {
//                    dcopy(nRowOrCol - i, matrix[i][i+1], dimension, matrix[i+1][i], 1);
//                }
//                for(int i = 0; i < nRowOrCol; i++) {
//                    dscal(nRowOrCol-i, -1.0, matrix[i+1][i], 1);
//                }
                break;
            default:
                System.out.println("Error opcion desconocida para SYMMAT");
                break;
        }
    }
//
//    /**
//     * copia un vector x a un vector y, usa loops unrolled para incrementos iguales a 1.
//     */
//    public void dcopy(int n, double[] dx, int incx, double[] dy, int incy) {
//
//        if(n < 0)
//            return;
//
//        int ix = incx < 0? (-n+1)*incx+1: 1;
//        int iy = incx < 0? (-n+1)*incy+1: 1;
//
//        /**
//         * para incrementos distintos o no iguales a 1
//         */
//        if(incx == 1 && incy == 1) {
//            for(int i = 0; i < n; i++, ix = ix+incx, iy = iy+incy)
//                dy[iy] = dx[ix];
//            return;
//        }
//
//        int m = n % 7;
//        if(m != 0) {
//            for( int i = 0; i < m; i++)
//                dy[i] = dx[i];
//            if(n < 7)
//                return;
//        }
//        int mp1 = m + 1;
//        for(int i = 0; i < n; i = i+7) {
//            dy[i] = dx[i];
//            dy[i+1] = dx[i+1];
//            dy[i+2] = dx[i+2];
//            dy[i+3] = dx[i+3];
//            dy[i+4] = dx[i+4];
//            dy[i+5] = dx[i+5];
//            dy[i+6] = dx[i+6];
//        }
//    }
//
//
//    public void dscal(int n, double sa, double []sx, int incx) {
//        if(n < 0 || incx < 0)
//            return;
//
//        int nincx = n*incx;
//        if(incx != 1) {
//            for(int i = 0; i < nincx; i = i + incx)
//                sx[i] = sa*sx[i];
//            return;
//        }
//
//        int m = n%5;
//        if(m!=0) {
//            for(int i = 0; i<m; i++)
//                sx[i] = sa*sx[i];
//            if(n<5)
//                return;
//        }
//        int mp1 = m + 1;
//        for(int i = mp1; i<n; i = i + 5) {
//            sx[i] = sa*sx[i];
//            sx[i+1] = sa*sx[i+1];
//            sx[i+2] = sa*sx[i+2];
//            sx[i+3] = sa*sx[i+3];
//            sx[i+4] = sa*sx[i+4];
//        }
//    }
//
//    public void dscal(int n, double sa, double sx, int incx) {
//        if(n < 0 || incx < 0)
//            return;
//
//        int nincx = n*incx;
//        if(incx != 1) {
//            for(int i = 0; i < nincx; i = i + incx)
//                sx = sa*sx;
//            return;
//        }
//
//        int m = n%5;
//        if(m!=0) {
//            for(int i = 0; i<m; i++)
//                sx = sa*sx;
//            if(n<5)
//                return;
//        }
//    }
//
//
//    public String printInterval() {
//        String result = "";
//
//
//
//        return result;
//    }

    /**
     * SUBROUTINE DIAMAT(MATRIX,DMAT,N,OPTION,EIGVAL)
     * @param matrix
     * @param dimension
     * @param nRowsOrCol
     * @param option
     * @return
     */
    public double[][] diagonalizeMatrix(Double[][] matrix, int dimension, int nRowsOrCol, String option) {
        double [] eigenValues = new double[dimension];
//
//        /**
//         * solo esta implementado con las opciones por defecto, ignora todo el demas procesamiento
//         */
//        int [] iw1 = new int[dimension];
//        double [][] rw = new double[dimension][dimension];
//        double [] rw1 = new double[dimension];
//        int jRot = 0;
//        jacobi(matrix, nRowsOrCol, dimension, iw1);
//        matrix = rw;
//
//        /**
//         * llena de 0 el espacion no usado de la matriz
//         */
//        if(nRowsOrCol < dimension)
//            fillMatrixWith(matrix, 0, nRowsOrCol, 0);
//        fillMatrixWith(matrix, 0, 0, nRowsOrCol);
//
//        return matrix;
        return null;
    }

    public Double[] getDiagonal(Double [][] matrix) {

        Double [] vector;
//        if(matrix[0][0] instanceof Double)
            vector = new Double[matrix.length];
        int nRows = matrix.length;
        int nCols = matrix[0].length;
        for(int i = 0; i < nRows; i ++)
            vector[i] = matrix[i][i];
        return vector;
    }

    public void jacobi(Double [][]h, int n, int ndim, int [] iq) {
        Double [] eigval;
        int iegen = 0;
        Double [][] u = new Double[ndim][ndim];
        Double [] x = new Double[ndim];
        int imi1, ipiv, jpiv;
        double cosine, hdiMin, hdTest, hii, hTemp, p, pyThag, sine, tang, xdif, xMax;

        /**
         * almacena los elementos de la diagonal
         */
        eigval = getDiagonal(h);
        if(iegen == 0)
            u = (Double[][]) buildIdentityMatrix(3, 1);
        int nr =0;

        if(n > 1) {
            /**
             * escanea en busca del elemento diagonal mas grande para cada columna
             * x[i] contiene el elemento mas grande de la columna
             * iq[i] contiene el second subscript que define la posicion del elemento
             */
            for(int i = 1; i < n; i++) {
                x[i] = 0.0;
                imi1 = i - 1;
                for(int j = 0; j < imi1; j++) {
                    if(x[i] < StrictMath.abs(h[i][j])) {
                        x[i] = StrictMath.abs(h[i][j]);
                        iq[i] = j;
                    }
                }
            }
            /**
             * establece el indicador de "apagado"
             */
            hdTest = 1e38;

            /**
             * encuentra el x[i] maximo para definirlo como el elemento pivote
             * y realiza una prueba de "fin del problema"
             */
            xMax = 0.0;
            for(int i = 1; i < n; i++) {
                if(xMax < x[i]) {
                    xMax = x[i];
                    ipiv = iq[i];
                    jpiv = i;
                }
            }
            double eps = 1e-10;
            if(xMax > 0.0) {//ir a 1000
                if(hdTest <= 0.0 && xMax <= hdTest) {// si no salta a 148
                    hdiMin = StrictMath.abs(eigval[1]);
                    for(int i = 1; i < n; i ++)
                        if(hdiMin > StrictMath.abs(eigval[i]))
                            hdiMin = StrictMath.abs(eigval[i]);
                    hdTest = hdiMin * eps;
//                    if(hdTest >= xMax) // salta a 1000
                }
//                148
                nr++;
//                xdif = eigval[ipiv][jpiv];
//                p = xdif/(2.0*h[ipiv][jpiv]);

//                tang = 1.0 / (p++);
            }
        }
    }


    public void molsym(Atomo[] atoms) {
        /**
         * calcula el tensor molecular de inercia
         */
        findSymmetry(atoms);

        /**
         * tolerancia de los eigenvalores degenerados para el tensor molecular de inercia
         */
//        Double[] tenval = new Double[3];
        int factor = (int)StrictMath.log10(sumSelfVector(tenval)/3) - 2;
        double tentol = Globals.cartesianCoordinatesTolerance * StrictMath.pow(10, StrictMath.max(0, factor));

        /**
         * usa la informacion de simetria de los eigenvalores del tensor de inercia molecular
         */
//        Double [][]tenVec = new Double[3][3];
        if(tenval[2] - tenval[0] < tentol)
            Globals.CUBIC = true;
        else if(tenval[2] - tenval[1] < tentol) {
            if(tenval[0] < tentol)
                Globals.LINEAR = true;
            transformCartesians(2, 0, 1, atoms, tenVec);
        } else {
//            System.out.println("TRANSFORM CART");
            transformCartesians(0, 1, 2, atoms, tenVec);
        }
        Double[][] auxMatrix = (Double[][]) buildIdentityMatrix(3, 1);
//        System.out.println("MAX C AXIS ROT DEGREE: " + Globals.MAX_C_AXIS_ROT_DEGREE);
        for(int i = 1; i < Globals.MAX_C_AXIS_ROT_DEGREE; i++) {
//            System.out.println("ITER CAXIS " + i);
            if(rotateVector(i+1, auxMatrix[2], atoms, 0)) addSec(i, auxMatrix[2]);
            if(rotateVector(i+1, auxMatrix[0], atoms, 0)) addSec(i, auxMatrix[0]);
            if(rotateVector(i+1, auxMatrix[1], atoms, 0)) addSec(i, auxMatrix[1]);
        }

        for(int i = 1; i < Globals.MAX_S_AXIS_ROT_DEGREE; i++) {
            if(rotateVector(i+1, auxMatrix[2], atoms, 1)) addSes(i, auxMatrix[2]);
            if(rotateVector(i+1, auxMatrix[0], atoms, 1)) addSes(i, auxMatrix[0]);
            if(rotateVector(i+1, auxMatrix[1], atoms, 1)) addSes(i, auxMatrix[1]);
        }

        if(sigma(auxMatrix[2], atoms)) addSig(auxMatrix[2]);
        if(sigma(auxMatrix[0], atoms)) addSig(auxMatrix[0]);
        if(sigma(auxMatrix[1], atoms)) addSig(auxMatrix[1]);

        if(Globals.HIGEST_C_AXIS_DEGREE >= 2 || Globals.HIGEST_S_AXIS_DEGREE >= 4) {
            Globals.MAXIS = true;
            findGoupsWithAxialSymmetry();
        } else
            lowsym();


        /**
         * encuentra los ejes S restantes
         */
        if(!Globals.LINEAR) {
            for(int i = 0; i < Globals.HIGEST_C_AXIS_DEGREE; i++)
                for(int j = 0; j < Globals.N_C_AXIS[i]; j++) {
                    double ax = (double) cAxisList[j][i].vector[0];
                    double ay = (double) cAxisList[j][i].vector[1];
                    double az = (double) cAxisList[j][i].vector[2];
                    Double[] auxVec = {ax, ay, az};
                    if(rotateVector(i, auxVec, atoms, 1))
                        addSes(i, auxVec);
                    if(rotateVector(2*i, auxVec, atoms, 1))
                        addSes(2*i, auxVec);
                }
        }

        /**
         * remover ejes S2 redundantes
         */
        if(Globals.N_S_AXIS[1] > 0) {
            Globals.N_S_AXIS[0] -= Globals.N_S_AXIS[1];
            Globals.N_S_AXIS[1] = 0;
        }

        System.out.println(System.getProperty("line.separator") + "\nC_AXIS");
        printMatrix(cAxisList);

        System.out.println(System.getProperty("line.separator") + "\nS_AXIS");
        printVector(sAxisList);

        System.out.println(System.getProperty("line.separator") + "\nSIGMA");
        printVector(mirrorPlaneList);
    }

    public void lowsym() {

    }

    public void findGoupsWithAxialSymmetry() {

    }


    public String printCoords(Atomo[] atoms) {
        DecimalFormat coordsFormat = new DecimalFormat("###0.000000");
        DecimalFormat massFormat = new DecimalFormat("####.###");

//        printMessage(header);

        int indxAtom = 0;
        String ret = "";
        for(Atomo a: atoms) {
            int asfd = 0;
            indxAtom++;
            String xVal = coordsFormat.format(a.coords.get(0));
            String yVal = coordsFormat.format(a.coords.get(1));
            String zVal = coordsFormat.format(a.coords.get(2));
            String massVal = massFormat.format(a.masa);
            ret += "\n" + Parser.repeat(" ", 5-String.valueOf(indxAtom).length()) + indxAtom +
                    "  " + a.label +
                    Parser.repeat(" ", 26-xVal.length()-8-a.label.length()) + xVal +
                    Parser.repeat(" ", 41-yVal.length()-26) + yVal +
                    Parser.repeat(" ", 56-zVal.length()-41) + zVal +
                    Parser.repeat(" ", 65-String.valueOf(a.atomicN).length()-56) + a.atomicN +
                    Parser.repeat(" ", 74-massVal.length()-65) + massVal +
                    Parser.repeat(" ", 80-String.valueOf(asfd).length()-74) + asfd
            ;
        }
        ret += "\n\n";
        return ret;
//        bw.write("\n\n");
    }

    public void addSec(int dimension, Double[] vector) {
        normalizeVector(vector);
        for(int i = 0; i < Globals.N_C_AXIS[dimension]; i++) {
            double secX = (double) cAxisList[i][dimension].vector[0];
            double secY = (double) cAxisList[i][dimension].vector[1];
            double secZ = (double) cAxisList[i][dimension].vector[2];
            double scapro = secX*vector[0] + secY*vector[1] + secZ*vector[2];
            if(StrictMath.abs( StrictMath.abs(scapro)-1.0) < 0.5*Globals.cartesianCoordinatesTolerance)
                return;
            if(dimension > Globals.HIGEST_C_AXIS_DEGREE)
                Globals.HIGEST_C_AXIS_DEGREE = dimension;
            if(Globals.N_C_AXIS[dimension] < Globals.MAX_C_AXIS) {
                Globals.N_C_AXIS[0] += 1;
                Globals.N_C_AXIS[dimension] += 1;
                cAxisList[Globals.N_C_AXIS[dimension]][dimension].vector[0] = vector[0];
                cAxisList[Globals.N_C_AXIS[dimension]][dimension].vector[1] = vector[1];
                cAxisList[Globals.N_C_AXIS[dimension]][dimension].vector[2] = vector[2];
                System.out.println("ADSEC X:" + vector[0] + " Y:" + vector[1] + " Z:" + vector[2]);
            } else
                System.out.println("Error addsec: macimo numero de ejes C excedido");
        }
    }


    public void addSes(int dimension, Double[] vector) {
        normalizeVector(vector);
        for(int i = 0; i < Globals.N_S_AXIS[dimension]; i++) {
            double secX = (double) sAxisList[i][dimension].vector[0];
            double secY = (double) sAxisList[i][dimension].vector[1];
            double secZ = (double) sAxisList[i][dimension].vector[2];
            double scapro = secX*vector[0] + secY*vector[1] + secZ*vector[2];
            if(StrictMath.abs( StrictMath.abs(scapro)-1.0) < 0.5*Globals.cartesianCoordinatesTolerance)
                return;
            if(dimension > Globals.HIGEST_S_AXIS_DEGREE)
                Globals.HIGEST_S_AXIS_DEGREE = dimension;
            if(Globals.N_S_AXIS[dimension] < Globals.MAX_C_AXIS) {
                Globals.N_S_AXIS[0] += 1;
                Globals.N_S_AXIS[dimension] += 1;
                cAxisList[Globals.N_S_AXIS[dimension]][dimension].vector[0] = vector[0];
                cAxisList[Globals.N_S_AXIS[dimension]][dimension].vector[1] = vector[1];
                cAxisList[Globals.N_S_AXIS[dimension]][dimension].vector[2] = vector[2];
            } else
                System.out.println("Error addsec: macimo numero de ejes C excedido");
        }
    }

    public void addSig(Double[] vector) {
        normalizeVector(vector);
        for(int i = 0; i < Globals.N_MIRROR_PLANES; i++) {
            double sigX = (double) mirrorPlaneList[i].vector[0];
            double sigY = (double) mirrorPlaneList[i].vector[1];
            double sigZ = (double) mirrorPlaneList[i].vector[2];
            double scapro = sigX*vector[0] + sigY*vector[1] + sigZ*vector[2];
            if(StrictMath.abs( StrictMath.abs(scapro)-1.0) < 0.5*Globals.cartesianCoordinatesTolerance)
                return;
            if(Globals.N_MIRROR_PLANES < Globals.MAX_MIRROR_PLANES) {
                Globals.N_MIRROR_PLANES++;
                mirrorPlaneList[Globals.N_MIRROR_PLANES].vector[0] = vector[0];
                mirrorPlaneList[Globals.N_MIRROR_PLANES].vector[1] = vector[1];
                mirrorPlaneList[Globals.N_MIRROR_PLANES].vector[2] = vector[2];
            } else
                System.out.println("Error addsec: macimo numero de ejes C excedido");
        }
    }

    /**
     *
     * @param baseAngle
     * @param vector
     * @param atoms
     * @param axisType 0 para ejes C, 1 para ejes S
     * @return
     */
    public boolean rotateVector(double baseAngle, Double[] vector, Atomo[] atoms, int axisType) {
//        System.out.println("ENTRA CAXIS " + baseAngle);
//        printVector(vector);
        boolean isOk = false;
        Double phi, ax = 0d, ay = 0d, az = 0d;
        Double[] vec = vector;
        if(!isZeroVector(vector)) {
            phi = 2*StrictMath.PI/baseAngle;
            vec = normalizeVector(vector);
            Double[][] rotationMatrix = buildRotationMatrix(phi, vec, (int) baseAngle);
//            if(baseAngle == 2d) {
//                System.out.println("\nROT matrix caxis");
//                printVector(vec);
//                printMatrix(rotationMatrix);
//            }
            int asdf  = 0;
            for(Atomo a: atoms) {
                ax = rotationMatrix[0][0]*a.coords.get(0) + rotationMatrix[0][1]*a.coords.get(1) + rotationMatrix[0][2]*a.coords.get(2);
                ay = rotationMatrix[1][0]*a.coords.get(0) + rotationMatrix[1][1]*a.coords.get(1) + rotationMatrix[1][2]*a.coords.get(2);
                az = rotationMatrix[2][0]*a.coords.get(0) + rotationMatrix[2][1]*a.coords.get(1) + rotationMatrix[2][2]*a.coords.get(2);
                if(baseAngle == 2d) {
                    System.out.println("CAXIS COOORDS");
                    System.out.println(a.coords.get(0) + " " + a.coords.get(1) + " " + a.coords.get(2));
                    System.out.println("ROTATION MATRIX");
                    printMatrix(rotationMatrix);
                    System.out.println("\taxyz " +ax + " " + ay + " " + az + " IATOM " + ++asdf + "\n");
//                    printVector(vec);
//                    printMatrix(rotationMatrix);
                }
            }
            Double[] reflected;
            if(axisType == 1) {
                Double[] originalVector = new Double[]{ax, ay, az};
                reflected = reflectVector(originalVector, vec);
            }
            if(determineEquivalentAtom(0, ax, ay, az, atoms) != 0)
                return false;
        }
        return true;
    }

    public boolean sigma(Double[] vector, Atomo[] atoms) {
        if(isZeroVector(vector))
            return false;
        int i = 0;
        for(Atomo a: atoms) {
            Double[] original = {a.coords.get(0), a.coords.get(1), a.coords.get(2)};
            Double[] reflected = reflectVector(original, vector);
            if(determineEquivalentAtom(i, reflected[0], reflected[1], reflected[2], atoms) != 0)
                return false;
        }
        return true;
    }

    public Double[] reflectVector(Double[] vector, Double[] planeNormal) {
        Double[] reflected = new Double[vector.length];
        if(isZeroVector(planeNormal))
            fillVectorWith(reflected, 0);
        else {
            normalizeVector(planeNormal);
            Double scapro = vector[0]*planeNormal[0] + vector[1]*planeNormal[1] + vector[2]*planeNormal[2];
            reflected[0] = vector[0] - 2*scapro*planeNormal[0];
            reflected[1] = vector[1] - 2*scapro*planeNormal[1];
            reflected[2] = vector[2] - 2*scapro*planeNormal[2];
        }
        return reflected;
    }

    public Double[] normalizeVector(Double[] vector) {
        if(isZeroVector(vector))
            return null;
        double magnitud = StrictMath.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
        vector[0] /= magnitud;
        vector[1] /= magnitud;
        vector[2] /= magnitud;
        return vector;
    }


    public boolean isZeroVector(Object [] vector) {
        for(Object o: vector)
            if(o == 0 || o == 0d)
                return true;
            return false;
    }


    public Double[][] buildRotationMatrix(double angle, Double[] vector, int index) {
        Double[] vec = vector;
        if(isZeroVector(vector))
            return (Double[][]) buildIdentityMatrix(3, 1);
        vec = normalizeVector(vector);
//        System.out.println("NORMALIZED");
//        printVector(vec);
        double phiCosine = StrictMath.cos(angle);
        double phiSine = StrictMath.sin(angle);
        double cosT = 1.0 - phiCosine;
        if(index == 2) {
            System.out.println("ANGLE " + angle + " cos " + phiCosine + " sin " + phiSine + " soct " + cosT);
            System.out.println("y " + vec[1] + " cosT " + cosT + " cos" + phiCosine);
        }
        return new Double[][]{
                {
                        vec[0] * vec[0] * cosT + phiCosine,
                        vec[0] * vec[1] * cosT - vec[2] * phiSine,
                        vec[0] * vec[2] * cosT + vec[1] * phiSine},
                {
                        vec[1] * vec[0] * cosT + vec[2] * phiSine,
                        vec[1] * vec[1] * cosT + phiCosine,
                        vec[1] * vec[2] * cosT - vec[0] * phiSine},
                {
                        vec[2] * vec[0] * cosT - vec[1] * phiSine,
                        vec[2] * vec[1] * cosT + vec[0] * phiSine,
                        vec[2] * vec[2] * cosT + phiCosine}
        };
    }


    public Double[][] buildRotationMatrix(double angle, Double[] vector) {
        Double[] vec = vector;
        if(isZeroVector(vector))
            return (Double[][]) buildIdentityMatrix(3, 1);
        vec = normalizeVector(vector);
//        System.out.println("NORMALIZED");
//        printVector(vec);
        double phiCosine = StrictMath.cos(angle);
        double phiSine = StrictMath.sin(angle);
        double cosT = 1.0 - phiCosine;
//        if(index == 2)
//            System.out.println("ANGLE " + angle + " cos " + phiCosine + " sin " + phiSine + " soct " + cosT);
        return new Double[][]{
                {
                        vec[0] * vec[0] * cosT + phiCosine,
                        vec[0] * vec[1] * cosT - vec[2] * phiSine,
                        vec[0] * vec[2] * cosT + vec[1] * phiSine},
                {
                        vec[1] * vec[0] * cosT + vec[2] * phiSine,
                        vec[1] * vec[1] * cosT + phiCosine,
                        vec[1] * vec[2] * cosT - vec[0] * phiSine},
                {
                        vec[2] * vec[0] * cosT - vec[1] * phiSine,
                        vec[2] * vec[1] * cosT + vec[0] * phiSine,
                        vec[2] * vec[2] * cosT + phiCosine}
        };
    }


    public int determineEquivalentAtom(int currentAtom, double ax, double ay, double az, Atomo[] atoms) {
        for(int i = 0; i < atoms.length; i++) {
            if( StrictMath.abs(Globals.ain[i])-Globals.ain[currentAtom] > 0.0
                    && StrictMath.abs(ax - atoms[i].coords.get(0)) > 5*Globals.cartesianCoordinatesTolerance
                    && StrictMath.abs(ay - atoms[i].coords.get(1)) > 5*Globals.cartesianCoordinatesTolerance
                    && StrictMath.abs(az - atoms[i].coords.get(2)) > 5*Globals.cartesianCoordinatesTolerance
                    )
                return i;
        }
        return 0;
    }

    public void transformCartesians (int ix, int iy, int iz, Atomo[] atoms, Double[][] tensorVector) {
        Double ax, ay, az;
//        System.out.println("TRACAR COORDS INI");
//        System.out.println(printCoords(atoms));
//        System.out.println("TENVEC");
//        printMatrix(tensorVector);
        for(Atomo a: atoms) {
            ax = a.coords.get(0);
            ay = a.coords.get(1);
            az = a.coords.get(2);
//            System.out.println("XCOORDS " + (tensorVector[0][0]*ax + tensorVector[1][0]*ay + tensorVector[2][0]*az));
            a.coords.set(ix, tensorVector[0][0]*ax + tensorVector[1][0]*ay + tensorVector[2][0]*az);
//            System.out.println("YCOORDS " + (tensorVector[0][1]*ax + tensorVector[1][1]*ay + tensorVector[2][1]*az));
            a.coords.set(iy, tensorVector[0][1]*ax + tensorVector[1][1]*ay + tensorVector[2][1]*az);
//            System.out.println("ZCOORDS " + (tensorVector[0][2]*ax + tensorVector[1][2]*ay + tensorVector[2][2]*az));
            a.coords.set(iz, tensorVector[0][2]*ax + tensorVector[1][2]*ay + tensorVector[2][2]*az);
        }

        Globals.inputToStandarMatrix[ix][0] = tensorVector[0][0];
        Globals.inputToStandarMatrix[ix][1] = tensorVector[1][0];
        Globals.inputToStandarMatrix[ix][2] = tensorVector[2][0];
        Globals.inputToStandarMatrix[iy][0] = tensorVector[0][1];
        Globals.inputToStandarMatrix[iy][1] = tensorVector[1][1];
        Globals.inputToStandarMatrix[iy][2] = tensorVector[2][1];
        Globals.inputToStandarMatrix[iz][0] = tensorVector[0][2];
        Globals.inputToStandarMatrix[iz][1] = tensorVector[1][2];
        Globals.inputToStandarMatrix[iz][2] = tensorVector[2][2];

        System.out.println("TRACAR COORDS");
        System.out.println(printCoords(atoms));
    }

    public double sumSelfVector(Double[] vector) {
        double sum = 0;
        for(int i = 0; i < vector.length; i++)
            sum += vector[i];
        return sum;
    }

    public void determineOrderPointGroup() {
        /**
         * contar todos los elementos de simetria de los grupos de simetria
         * primer E y todos los planos espejo
         */
        Globals.pointGroupOrder = 1 + Globals.N_MIRROR_PLANES;

        /**
         * iterar los ejes C
         */
        for(int i = 1; i < Globals.HIGEST_C_AXIS_DEGREE; i++)
            for(int j = 0; j < Globals.N_C_AXIS[i]; j++)
                for(int k = 0; k < i-1; k++)
                    if(areRedundant(j, i)) Globals.pointGroupOrder += 1;

        /**
         * itera los ejes S
         */
        int step = Globals.SIGMAH? 1: 2;

        for(int i = 1; i < Globals.HIGEST_S_AXIS_DEGREE; i++)
            for(int j = 0; j < Globals.N_S_AXIS[i]; j++)
                for(int k = 0; k < i-1; k++)
                    if(areRedundant(j, i)) Globals.pointGroupOrder += 1;
    }

    public boolean areRedundant(int se1, int se2) {
        for(int i = 0; i < StrictMath.min(se1, se2); i++)
            if(se1 % i == 0 && se2 % i == 0)
                return false;
        return true;
    }

    public void generateSymmetryList(Atomo[] atoms) {
        int step = Globals.SIGMAH? 1: 2;
        int iatom = 0;
        int iequatom;

        do {
            /**
             * itera los planos espejo
             */
            for(int i = 0; i < Globals.N_MIRROR_PLANES; i++){
                Double[] auxVector = {atoms[i].coords.get(0), atoms[i].coords.get(1), atoms[i].coords.get(2)};
                Double[] reflected = reflectVector(auxVector, (Double[]) mirrorPlaneList[i].vector);
                iequatom = determineEquivalentAtom(i, reflected[0], reflected[1], reflected[2], atoms);
                if(iequatom > 0 && !inSymList(iequatom)) {
                    Globals.groupUpperLimit[Globals.N_GRP_EQ_ATOMS]++;
                    Globals.equivalentAtomsGroupOrderedList[Globals.N_GRP_EQ_ATOMS] = iequatom;
                    Globals.nEquivalentAtoms[Globals.N_GRP_EQ_ATOMS]++;
                }
            }

            /**
             * usa todos los ejes C
             */
            for(int j = 0; j < Globals.HIGEST_C_AXIS_DEGREE; j++)
                for(int k = 0; k < Globals.N_C_AXIS[j]; k++) {
                    for(int l = 0; l < j-1; l++) {
                        double phi = 2.0*StrictMath.PI*l/j;
                        Double[] auxVector = {atoms[iatom].coords.get(0), atoms[iatom].coords.get(0), atoms[iatom].coords.get(0)};
                        Double[] rotatedVetor = doRotateVector(auxVector, phi, (Double[]) cAxisList[j][k].vector);
                        iequatom = determineEquivalentAtom(iatom, rotatedVetor[0], rotatedVetor[1], rotatedVetor[2], atoms);
                        if(iequatom > 0 && !inSymList(iequatom)) {
                            Globals.groupUpperLimit[Globals.N_GRP_EQ_ATOMS]++;
                            Globals.equivalentAtomsGroupOrderedList[Globals.N_GRP_EQ_ATOMS] = iequatom;
                            Globals.nEquivalentAtoms[Globals.N_GRP_EQ_ATOMS]++;
                        }
                    }
                }

            /**
             * usa todos los ejes S
             */
            for(int j = 0; j < Globals.HIGEST_S_AXIS_DEGREE; j++)
                for(int k = 0; k < Globals.N_S_AXIS[j]; k++) {
                    for(int l = 0; l < j-1; l++) {
                        double phi = 2.0*StrictMath.PI*l/j;
                        Double[] auxVector = {atoms[iatom].coords.get(0), atoms[iatom].coords.get(0), atoms[iatom].coords.get(0)};
                        Double[] rotatedVetor = doRotateVector(auxVector, phi, (Double[]) sAxisList[j][k].vector);
                        Double[] reflected = reflectVector(rotatedVetor, (Double[]) sAxisList[j][k].vector);
                        iequatom = determineEquivalentAtom(iatom, reflected[0], reflected[1], reflected[2], atoms);
                        if(iequatom > 0 && !inSymList(iequatom)) {
                            Globals.groupUpperLimit[Globals.N_GRP_EQ_ATOMS]++;
                            Globals.equivalentAtomsGroupOrderedList[Globals.N_GRP_EQ_ATOMS] = iequatom;
                            Globals.nEquivalentAtoms[Globals.N_GRP_EQ_ATOMS]++;
                        }
                    }
                }

            System.out.println("GRPUPPLIM " + Globals.groupUpperLimit.length + " atoms.length " + atoms.length + " NGRP " + Globals.N_GRP_EQ_ATOMS + "     " + Globals.groupUpperLimit[Globals.N_GRP_EQ_ATOMS-1]);
            if(Globals.groupUpperLimit[Globals.N_GRP_EQ_ATOMS-1] != atoms.length) {
//              if(2 != atoms.length) {
//            if(Globals.groupUpperLimit[0] != 1) {
                for(int i = 0; i < atoms.length; i++)
                    if(!inSymList(i)) {
                        iatom = i;
                        Globals.N_GRP_EQ_ATOMS++;
                        Globals.groupLowerLimit[Globals.N_GRP_EQ_ATOMS-1] = Globals.groupUpperLimit[Globals.N_GRP_EQ_ATOMS-2] + 1;
                        Globals.groupUpperLimit[Globals.N_GRP_EQ_ATOMS-1] = Globals.groupLowerLimit[Globals.N_GRP_EQ_ATOMS-1];
                        Globals.equivalentAtomsGroupOrderedList[Globals.groupLowerLimit[Globals.N_GRP_EQ_ATOMS-1]] = iatom;
                        break;
                    }
            }
        } while (Globals.groupUpperLimit[Globals.N_GRP_EQ_ATOMS-1] != atoms.length);

        /**
         * genera el vector GRPOF
         */
        for(int i = 0; i < Globals.N_GRP_EQ_ATOMS; i++)
            for(int j = Globals.groupLowerLimit[i]; j < Globals.groupUpperLimit[i]; j++) {
                Globals.equivalentAtomsGroup[Globals.equivalentAtomsGroupOrderedList[Globals.N_GRP_EQ_ATOMS-1]] = iatom;
            }
    }

    public Double[] doRotateVector(Double[] vector, double angle, Double[] plane) {
        if(isZeroVector(plane))
            return new Double[] {0d,0d,0d};
        Double[][] rotationMatrix = buildRotationMatrix(angle, plane);
        return new Double[] {rotationMatrix[0][0]* vector[0] + rotationMatrix[0][1]* vector[1] + rotationMatrix[0][2] * vector[2],
                rotationMatrix[1][0]* vector[0] + rotationMatrix[1][1]* vector[1] + rotationMatrix[1][2] * vector[2],
                rotationMatrix[2][0]* vector[0] + rotationMatrix[2][1]* vector[1] + rotationMatrix[2][2] * vector[2]
        };
    }

    public boolean inSymList(int atomA) {
        for(int i = 0; i < Globals.groupUpperLimit[Globals.N_GRP_EQ_ATOMS-1]; i++)
            if(Globals.equivalentAtomsGroupOrderedList[i] == atomA)
                return true;
        return false;
    }

    /**
     * SYMTRA
     */
    public void generetaSymmetryTransformationMatrix() {
        Globals.standardToInputMatrix = Globals.inputToStandarMatrix;

        System.out.println(System.getProperty("line.separator") + "Standar to input");
        printMatrix(Globals.standardToInputMatrix);
        /**
         * calcular la matriz inversa de inputToStandar, es decir standardToInput
         */
        double determinant = invertMatrix(Globals.standardToInputMatrix);

        System.out.println("DETER " + determinant);
        System.out.println(System.getProperty("line.separator") + "Input to standar");
        printMatrix(Globals.standardToInputMatrix);
        /**
         * checa el determinante resultado de la eliminacion Gauss Jordan
         */
//        if(StrictMath.abs(determinant) < tiniest)
            System.out.println("SYMTRA error al obtener la matriz inversa");
    }


    public Double invertMatrix(Double[][] matrix) {
        double eps = 1e-10;
        double determinant = 1;
        int n  = matrix.length;
        Integer [] l = new Integer[3*n];
        Integer [] m = new Integer[3*n];
        double bigA;

        for(int i = 0; i < n; i++) {
            l[i] = i;
            m[i] = i;
            bigA = matrix[i][i];

            /**
             * encuentra el valor mas grande para la fila y columna
             */
            for(int j = 0; j < n ; j++)
                for(int k = 0; k < n; k++){
                    if(StrictMath.abs(bigA) < StrictMath.abs(matrix[j][k])) {
                        l[i] = k;
                        m[i] = j;
                        bigA = matrix[j][k];
                    }
                }

            /**
             * intercambia las columnas
             */
            int aux = l[i];
            double temp;
            if(aux > i) {
                for(int j = 0; j < n; j++){
                    temp = -matrix[j][i];
                    matrix[j][i] = matrix[j][aux];
                    matrix[j][aux] = temp;
                }
            }

            /**
             * intercambia las filas
             */
            int auxRow = l[i];
            if(auxRow > i) {
                for(int j = 0; j < n; j++){
                    temp = -matrix[i][j];
                    matrix[i][j] = matrix[auxRow][j];
                    matrix[auxRow][j] = temp;
                }
            }

            if(StrictMath.abs(bigA) < eps)
                return 0d;

            /**
             * divide usando el elemento pivote
             */
            for(int j = 0; j < n; j++)
                if(j != i)
                    matrix[i][j] /= -bigA;

            /**
             * reduce filas y columnas excepto para el pivote
             */
            for(int j = 0; j < n; j++)
                for(int k = 0; k < n; k++)
                    if(j != i && k != i)
                        matrix[k][j] = matrix[i][j]*matrix[k][i] + matrix[k][j];

            /**
             * divide usando el elemento pivote
             */
            for(int j = 0; j < n; j++)
                if(j != i)
                    matrix[j][i] /= bigA;
            determinant = StrictMath.max(-1e25, StrictMath.min(1e25, determinant));
            determinant *= bigA;
            matrix[i][i] = 1/bigA;
        }

        /**
         * realize la sustitucion por atras
         */
        int k = n;
        k--;
        double temp;
//        if(i <= 0)
//            return determinant;
//        int i = l[k];
//        if(i > k)
//            for(int j = 0; j < n; j++){
//                temp = matrix[k][j];
//                matrix[k][j] = -matrix[i][j];
//                matrix[i][j] = temp;
//            }
//
//        int j = m[k];
//        if(j <= k)
//        return determinant;
//        for(i = 0; i < n; i++){
//            temp = matrix[i][k];
//            matrix[i][k] = -matrix[i][j];
//            matrix[i][j] = temp;
//        }
        return determinant;
    }

    public void transformDummyAtoms() {
//        if(Globals.N_DUMMY_ATOMS == 0) return;
//
//        /**
//         * trasladar el centro de masa de los atomos dummies
//         */
//        for(Atomo a: Atomo.listDummies) {
//            shiftVector(a.coords, shiftVector);
//            rotateVector(a.coords, matrix);
//        }
    }


    public void symmetryDriver() {
        Globals.COORDINATE_SET = 1;
    }


    /**
     * GETPG
     * obtener la informacion del grupo puntual
     */
    public void findPointGroupInfo() {

    }
}

class Jacobi {

//    Double [][]matriz={{(double) 4, (double) -2, (double) 1},{(double) 1, (double) -5, (double) 32},{(double) 2, (double) 1, (double) 4}};
    Double [][]matriz={{0.9583909562031309, 0d, 0.2526453923123321},{0d, 1.7863158230614995, 0d},{0.2526453923123321, 0d, 0.8279248668583686}};
    double []vector={2,1,3};
    double []vectorR={1,2,3};
    double []x2=vectorR;
    double sumatoria=1;
    int max=50;

//    public Jacobi(Double [][] matriz, Double eigval) {
//
//    }

    public void SolJacobi(){
        int tam = matriz.length;
        for (int y = 0; y < 50; y++) {
            System.out.println("\nvector " + y);
//            for(int t=0;t<max;t++){
                x2=vectorR.clone();
                for (int i = 0; i < tam; i++) {
                    sumatoria=0;
                    for (int s = 0; s < tam; s++) {
                        if(s!=i)sumatoria += matriz[i][s]*x2[s];
                    }
                    vectorR[i]=(vector[i]-sumatoria)/matriz[i][i];
                    System.out.print(" " + vectorR[i]);
                }
//            }

        }
    }
    public static void main(String[] args) {
        Jacobi obj=new Jacobi();
        obj.SolJacobi();
    }




}
