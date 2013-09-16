import demon.*;

import java.io.IOException;

public class Main {

	/**
	 * @param args la ubicación del archivo de entrada
	 */
	public static void main(String[] args) throws Throwable {

        DeMon deMon = new DeMon();

        deMon.initMPIExecution();
        deMon.startTimer();
        deMon.timer.doMark(false);
        deMon.printCopyright();
        deMon.printer.printMessage(deMon.timer.printLastStep("FIRST TIMING"));
        deMon.timer.doMark(false);
        deMon.printer.printMessage(deMon.timer.printLastStep("INITIALIZATION"));

        try {
            deMon.parser.parseGeometry();
            deMon.establishAtomsArray();
        } catch (Exception e) {
            System.out.println("Error parser:" + e.getMessage());
            e.printStackTrace();
            return;
        }
        deMon.printInitResume();
        deMon.timer.doMark(false);
        try {
            deMon.printInputGeometry();
        } catch (IOException e) {
            System.out.println("Error printing geometry:" + e.getMessage());
            return;
        }
        deMon.printer.printMessage(deMon.timer.printLastStep("INPUT"));


//        System.out.println("TESTSTEEST");
        deMon.doSymmetryAnalysis();
        deMon.timer .doMark(false);
        deMon.printSymmetryAnalisis();
        deMon.printer.printMessage(deMon.timer.printLastStep("SYMMETRY ANALISIS"));
//        deMon.printer.printMessage(deMon.printSymmetryResults());
//        deMon.timer.printLastStep("SYMMETRY ANALISIS");



        deMon.parser.finalize();
        deMon.printer.finalize();
	}

}
