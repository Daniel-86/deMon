package demon;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Date;

/**
 * User: daniel
 * Date: 1/04/13
 * Time: 11:42 PM
 */
public class Timer {
    String name;
    int status;
    static int NEW = 0, STARTED = 1, CONTINUED = 2, SUSPENDED = 3, STOPPED = 4, RUNNING = 5;
    static DecimalFormat timeElapsedFormat = new DecimalFormat("####.###");
    ArrayList<Long> steps;
    double precision;
    static double SECONDS = 1e-3;

    public Timer(String name) {
        this.name = name;
        this.status = Timer.NEW;
        this.steps = new ArrayList<>();
        this.precision = SECONDS;
    }

    /**
     * Cambia a estatus RUNNING y almacena el tiempo de inicio (si no ha sido iniciado).
     */
    public double start() {
//        System.out.println("STARTING TIMER");
        if(steps.size() < 1)
            steps.add(System.currentTimeMillis());
        status = RUNNING;
        return steps.get(0) * precision;
    }

    /**
     * Registra el tiempo actual y devuelve el tiempo transcurrido desde el tiempo inicial o la marca anterior.
     * @param justStep indica si devolver el tiempo transcurrido desde la marca anterior o el total.
     * @return el tiempo transcurrido.
     */
    public double doMark(boolean justStep) {
        double elapsedTime = -1;
        if(status == RUNNING) {
            Long mark = System.currentTimeMillis();
            steps.add(mark);
            int nSteps = steps.size() - 1;
            elapsedTime = mark - ((justStep && nSteps > 1)? steps.get(nSteps): steps.get(0));
            elapsedTime *= precision;
        }
        return elapsedTime;
    }

    /**
     * Registra tiempo y cambia el estatus a STOPPED.
     * @return el tiempo total del timer.
     */
    public double doStop() {
        double totalTime = doMark(false);
        status = STOPPED;
        return totalTime;
    }

    public void clear() {
        steps.clear();
        status = NEW;
    }

    /**
     * El tiempo transcurrido entre 2 marcas.
     * @param start marca inicial.
     * @param end marca final.
     * @return el tiempo transcurrido o -1 y hubo algún error.
     */
    public double getInterval(int start, int end) {
        if(start < 0 || end <= start || steps.size() < end)
            return -1;
        return (steps.get(end) - steps.get(start)) * precision;
    }


    public String print(String message) {
        System.out.println("PRINTING TIMER" + status);
        if(status == NEW)
            return "El timer está vacío";
        String today = new Date().toString();
//        String elapsed = timeElapsedFormat.format(getInterval(0, steps.size() - 1));
//        return "" +
//                " " + name + Parser.repeat(" ", 57 - today.length() - 10 - name.length()) + "DATE: " + today +
//                "  TOTAL TIME:" +  Parser.repeat(" ", 12 - elapsed.length()) + elapsed +
//                "\n";
        String elapsedTotal = timeElapsedFormat.format(getInterval(0, steps.size()-1));
        String step = "DATE: ";
        String total = "  TOTAL TIME:";
        return "" +
                " " + message + Parser.repeat(" ", 55-message.length()-2 - step.length()-today.length()) +
                step + today +
                total + Parser.repeat(" ", 80 - 55 - total.length() - elapsedTotal.length()) + elapsedTotal +
                "\n\n";
    }

    /**
     * Imprime un resumen del intervalo. La vida del timer se divide en marcas, lo intervalos estan delimitados por las marcas.
     * @param start marca inicial.
     * @param end marca final.
     * @param message titulo para el resumen.
     * @return cadena que contiene el resume.
     */
    public String printInterval(int start, int end, String message) { // 2 35  57
        if(status == NEW || steps.size() < 1 || steps.size() <= end)
            return "El timer está vacío";
        String elapsedStep = timeElapsedFormat.format(getInterval(start, end));
        String elapsedTotal = timeElapsedFormat.format(getInterval(0, end));
        String step = "STEP TIME:";
        String total = "  TOTAL TIME:";
        return "" +
                " " + message + Parser.repeat(" ", 35-message.length()-2) +
                step + Parser.repeat(" ", 55 - 35 - step.length() - elapsedStep.length()) + elapsedStep +
                total + Parser.repeat(" ", 80 - 55 - total.length() - elapsedTotal.length()) + elapsedTotal +
                "\n\n";
    }

    /**
     * Imprime un resumen del intervalo. La vida del timer se divide en marcas, lo intervalos estan delimitados por las marcas.
     * @param start marca inicial.
     * @param message titulo para el resumen.
     * @return cadena que contiene el resume.
     */
    public String printInterval(int start, String message) {
        return printInterval(start, start+1, message);
    }

    public String printLastStep(String message) {
        return printInterval(steps.size()-2, message);
    }
}
