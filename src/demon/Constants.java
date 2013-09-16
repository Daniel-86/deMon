package demon;

/**
 * User: daniel
 * Date: 8/04/13
 * Time: 10:48 PM
 */
public class Constants {

    static double degreeToRadian (double degree) {
        return degree * StrictMath.PI / 180;
    }

    static double radianToDegree (double radian) {
        return radian * 180 / StrictMath.PI;
    }

    static long factorial (int number) {
        return 1;
    }

    static double factorialD (double number) {
        return 1.0;
    }

    static final double LIGHT_SPEED = 2.99792458e8;
    static final double CHARGE_ELEM = 1.602176531e-19;
    static final double ELECTRON_MASS = 9.109382616e-31;
    static final double E_FIELD = 8.854187817e-12;
    static final double PLANK = 6.626069311E-34;
    static final double BOLTZ = 1.380650524E-23;
    static final double PERM_VACCUM = 4.0 * StrictMath.PI * 1.0E-7;
    static final double AVOGADRO = 6.022141510E23;
    static final double PPM = 1.0/1000000.0;
    static final double PRESSURE_STD = 100000.0;
    static final double AFINE = 0.5 * PERM_VACCUM * LIGHT_SPEED * CHARGE_ELEM * CHARGE_ELEM / PLANK;
    static final double CONST_GAS = AVOGADRO * BOLTZ;
    static final double CONST_RYDBERG = 0.5 * ELECTRON_MASS * LIGHT_SPEED * AFINE * AFINE / PLANK;
    static final double BOHR_RADIUS = AFINE / (4 * StrictMath.PI * CONST_RYDBERG);
    static final double AMU_TO_KG = 1.66053886E-27;
    static final double AMU_TO_AU = AMU_TO_KG / ELECTRON_MASS;
    static final double ANGST_TO_BOHR = 1e-10 / BOHR_RADIUS;
    static final double FSEC_TO_AU = 4 * StrictMath.PI * CONST_RYDBERG * LIGHT_SPEED *1e-15;
    static final double HRTREE_TO_JOULE = 2 * CONST_RYDBERG * PLANK * LIGHT_SPEED;
    static final double HRTREE_TO_KJMOL = 0.001 * HRTREE_TO_JOULE * AVOGADRO;
    static final double HRTREE_TO_KCALMOL = HRTREE_TO_KJMOL / 4.184;
    static final double HRTREE_TO_EV = HRTREE_TO_JOULE / CHARGE_ELEM;
    static final double HRTREE_TO_HZ = HRTREE_TO_JOULE / PLANK;
    static final double HRTREE_TO_MHZ = 1e6 * HRTREE_TO_HZ;
    static final double HRTREE_TO_WAVENUM = 0.02 * CONST_RYDBERG;
    static final double WAVE_SEC = 100 * LIGHT_SPEED;
    static final double VIBFAC = 5 * StrictMath.sqrt(HRTREE_TO_KJMOL) / (StrictMath.PI * BOHR_RADIUS * LIGHT_SPEED);
    static final double HRTREE_TO_ESU = 1e21 * BOHR_RADIUS * LIGHT_SPEED * CHARGE_ELEM;
}
