package demon;

/**
 * User: daniel
 * Date: 8/04/13
 * Time: 02:21 PM
 */
public class PhysicConstants {
    static int MAX_L_BAS = 4; //Máximo número de cuanto L (para orbitales base).
    static int MAX_DER = 1; //Máximo rango para las derivadas de energía analítica.
    static int MAX_L_AUX = 6; //Número máximo del cuanto L (para las funciones auxiliares).
    static int MAX_L_UL = MAX_L_AUX + MAX_DER + 2 * MAX_L_BAS; // Número máximo del cuanto L (general).
    static int MAX_L_PTR = StrictMath.max(5, MAX_L_UL); // Número máximo del cuanto L (para la definición del puntero).
    static int MAX_L_BFN = StrictMath.max(MAX_DER, MAX_L_BAS);
    static int MAX_L_ECP = 5;
    static int MAX_L_TRA = StrictMath.max(MAX_L_ECP, MAX_L_BAS) + MAX_L_BFN;
    static double[] covalentRadii = new double[112];
    static double[] vanDerWaalsRadii = new double[112];
    static double[] coefficient6 = new double[112];
    static String[] configuration = new String[112];
    static String[] group = new String[112];
    static String[] atomSymbol = new String[112];
    static double[] mass = new double[112];
    static String[] orbitalSymbol = new String[13];
    static String[] coordinateLabel = new String[3];
    static double LIGHT_SPEED = 2.99792458e8;
    static double ELEMENTARY_CHARGE = 1.6021765314e-19;
    static double ELECTRON_MASS = 9.109382616e-31;
    static double ELECTRIC_FIELD = 8.854187817e-12;
    static double PLANK = 6.626069311e-34;
    static double BOLTZMANN = 1.380650524e-23;
    static double PERMEBILITy = 4.0 * StrictMath.PI * 1e-7;
    static double AVOGADRO = 6.022141510e23;
    static double PARTS_PER_MILLION = 1/1000000.0;

    static int[] NCO = buildNCO(); //Número de orbitales cartesianos.
    static int[] NCO_SET = buildNCOSet(); //Número de orbitales cartesianos para funciones auxiliares.
    static int[][][] COP = new int[MAX_L_PTR][MAX_L_PTR][MAX_L_PTR];
    static int[][][] SOP = new int[MAX_L_PTR][MAX_L_PTR][MAX_L_PTR];

    private static int[] buildNCO() {
        int[] nco = new int[MAX_L_PTR];
        for(int i = 0; i < MAX_L_PTR; i++)
            nco[i] = (i + 1) * (i + 2) / 2;
        return nco;
    }

    private static int[] buildNCOSet() {
        int[] nco_set = new int[MAX_L_PTR];
        nco_set[0] = 1;
        for(int i = 1; i < MAX_L_PTR; i++) {
            nco_set[i] = nco_set[i-1] + NCO[i];
        }
        return nco_set;
    }

    private static int[][][] buildCOPandSOP() {
//        int[][][] cop = new int[MAX_L_PTR][MAX_L_PTR][MAX_L_PTR];
//        int[][][] sop = new int[MAX_L_PTR][MAX_L_PTR][MAX_L_PTR];
        int LA;
        for(int i = 0; i < MAX_L_PTR; i++)
            for(int j = 0; j < MAX_L_PTR; j++)
                for(int k = 0; k < MAX_L_PTR; k ++) {
                    LA = i + j + k;
                    if(LA > MAX_L_PTR)
                        continue;
                    SOP[i][j][k] = 1 + k + (LA-i)*(LA-i+1) / 2;
                    COP[i][j][k] = LA == 0? SOP[i][j][k]: NCO_SET[LA-1] + SOP[i][j][k];
                }
        return SOP;
    }

    private int[][] buildSPH() {
        int[][] sph = new int[MAX_L_PTR][2*MAX_L_PTR+1];
        for(int i = 0; i < MAX_L_PTR; i++)
            for(int j = 1; j < 2*i+1; j++)
                sph[i][j] = i*i + j;
        return sph;
    }

    private int[][] buildSPP() {
        int[][] spp = new int[MAX_L_PTR][MAX_L_PTR];
        for(int i = 0; i < MAX_L_PTR; i++)
            for(int j = 1; j < i; j++)
                spp[i][j] = i * (i+1)/2 + j;
        return spp;
    }

    private int[][] buildPOP() {
        int[][] pop = new int[2*MAX_L_BFN][2*MAX_L_BFN];
        int LA, lb;
        for(int i = 0; i < 2*MAX_L_BFN; i++)
            for(int j = 0; j < 2*MAX_L_BFN; j ++)
                for(int k = 0; k < 2*MAX_L_BFN; k++) {
                    LA = i+j+k;
                    if(LA > 2*MAX_L_BFN)
                        continue;
                    for(int l = 0; l < 2*MAX_L_BFN; l++)
                        for(int m = 0; m < 2*MAX_L_BFN; m++)
                            for(int n = 0; n < 2*MAX_L_BFN; n++) {
                                lb = l+m+n;
                                if(lb > 2*MAX_L_BFN)
                                    continue;
                                pop[SOP[i][j][k]] [COP[l][m][n]] = COP[l][m][n] + SOP[i][j][k] + NCO_SET[lb];
                            }
                }
        return pop;
    }

    static int factorial(int n) {
//        if(n==1)
//            return 1;
//        return n*factorial(n-1);
        return 1;
    }

    static double factorialD(int n) {
//        if(n==1)
//            return 1;
//        return n*factorial(n-1);
        return 1.0;
    }

    static double noverk(int i, int j) {
        return 1.0;
    }

    /**
     * Matrices de transformación de coordenadas cartesianas a esféricas.
     * @return
     */
    private double[][][] cartToSphe() {
        double [][][] cartToSphe = new double[MAX_L_TRA*2][20*MAX_L_TRA+1][MAX_L_PTR*2];
        double [][][] polToSphe = new double[MAX_L_TRA*2][20*MAX_L_TRA+1][MAX_L_PTR*2];
        double [][][] spheToCart = new double[MAX_L_TRA*2][20*MAX_L_TRA+1][MAX_L_PTR*2];
        for(int i = 0; i < MAX_L_TRA; i++)
            for(int j = 1; j < 2*i+1; j++)
                for(int k = 1; k < NCO[i]; k++) {
                    cartToSphe[j][k][i] = 0.0;
                    polToSphe[j][k][i] = 0.0;
                }
        for(int i = 0; i < MAX_L_BAS; i++)
            for(int j = 1; j < 2*i+1; j++)
                for(int k = 1; k < NCO[i]; k++) {
                    cartToSphe[j][k][i] = 0.0;
                    polToSphe[j][k][i] = 0.0;
                }

        for(int l = 0; l < MAX_L_TRA; l++)
            for(int lx = l; lx > 0; lx--) {
                for (int ly = l-lx; ly > 0; ly--) {
                    int lz = l - lx - ly;
                    int ic = SOP[lx][ly][lz];
                    for (int m = -l; m < l; m++) {
                        int ls = l + m + 1;
                        int ma = StrictMath.abs(m);
                        int j = lx + ly - ma;
                        if(j >= 0 && j%2 == 0) {
                            lx /= 2;
                            double s1 = 0.0;
                            for(int i = 0; i < (l-ma)/2; i++) {
                                double s2 =0.0;
                                for(int k = 0; k < j; k++) {
                                    int expo;
                                    double s;
                                    if(m<0 && StrictMath.abs(ma-lx)%2 == 1 ||
                                            m>0 && StrictMath.abs(ma-lx)%2 == 0) {
                                        expo = (ma - lx + 2*k) / 2;
                                        s = StrictMath.pow(-1.0, expo) * StrictMath.sqrt(2.0);
                                    } else if( m==0 && lx%2==0) {
                                        expo = k - lx/2;
                                        s = StrictMath.pow(-1.0, expo);
                                    } else
                                        s = 0.0;
                                    s2 += noverk(i, j) * noverk(ma, lx-2*k) * s;
                                }
                                s1 += noverk(i, j) * noverk(i, j) * StrictMath.pow(-1.0, i) * factorial(2*l*2*i)/factorial(l-ma-2*i) * s2;
                            }
                            cartToSphe[ls][ic][l] = StrictMath.sqrt(factorial(2*lx) * factorial(2*ly) * factorial(2*lz) * factorial(l) * factorial(l-ma) /
                                    factorial(lx) * factorial(ly) * factorial(lz) * factorial(2*l) * factorial(l+ma)) *
                                    s1/StrictMath.pow(2.0, l)*factorial(l);
                        } else
                            cartToSphe[ls][ic][l] = 0.0;
                    }
                }
            }

        for(int l = 0; l < MAX_L_BAS; l++)
            for(int lx1 = l; lx1 > 0; lx1--)
                for(int ly1 = l-lx1; ly1 > 0; ly1--) {
                    int lz1 = l -lx1 - ly1;
                    int ic1 = SOP[lx1][ly1][lz1];
                    double s1 = StrictMath.sqrt(factorial(lx1) * factorial(ly1) * factorial(lz1) / factorial(2 * lx1) * factorial(2 * ly1) * factorial(2 * lz1));
                    for(int lx2 = l; lx2 > 0; lx2--)
                        for(int ly2 = l-lx2; ly2 > 0; ly1--) {
                            int lz2 = l - lx2 - ly2;
                            int ic2 = SOP[lx2][ly2][lz2];
                            int lx = lx1 + lx2;
                            int ly = ly1 + ly2;
                            int lz = lz1 + lz2;
                            if(lx%2 == 0 && ly%2 == 0 && lz%2 == 0) {
                                double s2 = StrictMath.sqrt(factorial(lx2) * factorial(ly2) * factorial(lz2) / factorial(2*lx2) * factorial(2*ly2) * factorial(2*lz2));
                                double s = StrictMath.sqrt(factorial(lx) * factorial(ly) * factorial(lz) * s1 * s2 / factorial(lx/2) * factorial(ly/2) * factorial(lz/2));
                                for(int is = 1; is < 2*l+1; is++)
                                    spheToCart[is][ic1][l] += s*cartToSphe[is][ic2][l];
                            }
                        }
                }

        for(int l1 = 0; l1 < MAX_L_TRA; l1++)
            for(int lx1 = 0; lx1 < l1; lx1++)
                for(int ly1 = 0; ly1 < l1-lx1; ly1++) {
                    int is;
                    double xyz = 0;
                    int lz1 = l1 -lx1 - ly1;
                    int ic1 = COP[lx1][ly1][lz1];
                    for(int l2 = 0; l2 < l1; l2++)  {
                        double s1 = 4 * StrictMath.PI * factorialD(2*l2+1);
                        for(is = 1; is < 2*l2+1; is ++) {
                            xyz = 0.0;
                            for(int lx2 = 0; lx2 < l2; lx2++)
                                for(int ly2 = 0; ly2 < l2-lx2; ly2++) {
                                    int lz2 = l2 - lx2 - ly2;
                                    int ic2 = SOP[lx2][ly2][lz2];
                                    int lx = lx2 + lx1;
                                    int ly = ly2 + ly1;
                                    int lz = lz1 + lz2;
                                    if(lx%2 == 0 && ly%2 == 0 && lz%2 == 0) {
                                        double s2 = factorialD(lx-1) * factorialD(ly-1) * factorialD(lz-1) / factorialD(l1+l2+1);
                                        double s = s1 / (factorialD(2*lx2-1) * factorialD(2*ly2-1) * factorialD(2*lz2-1));
                                        xyz += cartToSphe[is][ic2][l2] * StrictMath.sqrt(s) * s2;
                                    }
                                }
                        }
                        polToSphe[is][ic1][l2] = xyz;
                    }
                }

        return cartToSphe;
    }

    public PhysicConstants() {
        cartToSphe();
    }
}
