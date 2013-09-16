package demon;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

/**
 * User: daniel
 * Date: 4/06/13
 * Time: 12:15 AM
 */
public class DeMonPrinter {

    BufferedWriter bw;
    String fileName = "deMon.out";


    DecimalFormat coordsFormat = new DecimalFormat("###0.000000");
    DecimalFormat massFormat = new DecimalFormat("####.###");


    public DeMonPrinter() throws IOException {
        bw = new BufferedWriter(new FileWriter(fileName));   // Es el archivo de entrada cargado en memoria
    }

    public DeMonPrinter(String fileName) throws IOException {
        this.fileName = fileName;
        bw = new BufferedWriter(new FileWriter(fileName));   // Es el archivo de entrada cargado en memoria
    }

    public void finalize() throws Throwable {
        if(bw != null)
            bw.close();
        bw = null;
        super.finalize();
    }

    public void copyright() throws IOException {
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
//        bw.write(copyright);
        printMessage(copyright);
    }

    public void printMessage(String message) throws IOException {
        bw.write(message);
    }


    public void printCoords(String header, Atomo[] atoms) throws IOException {

        printMessage(header);

        int indxAtom = 0;
        for(Atomo a: atoms) {
            int asfd = 0;
            indxAtom++;
            String xVal = coordsFormat.format(a.coords.get(0));
            String yVal = coordsFormat.format(a.coords.get(1));
            String zVal = coordsFormat.format(a.coords.get(2));
            String massVal = massFormat.format(a.masa);
            bw.write("\n" + Parser.repeat(" ", 5-String.valueOf(indxAtom).length()) + indxAtom +
                    "  " + a.label +
                    Parser.repeat(" ", 26-xVal.length()-8-a.label.length()) + xVal +
                    Parser.repeat(" ", 41-yVal.length()-26) + yVal +
                    Parser.repeat(" ", 56-zVal.length()-41) + zVal +
                    Parser.repeat(" ", 65-String.valueOf(a.atomicN).length()-56) + a.atomicN +
                    Parser.repeat(" ", 74-massVal.length()-65) + massVal +
                    Parser.repeat(" ", 80-String.valueOf(asfd).length()-74) + asfd
            );
        }

        bw.write("\n\n");
    }
}
