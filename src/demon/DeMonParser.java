package demon;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * User: daniel
 * Date: 4/06/13
 * Time: 01:50 AM
 */
public class DeMonParser {

    BufferedReader br;
    String fileName = "deMon.inp";

    public ArrayList<Atomo> atomos;
    public ArrayList<Atomo> atomosProd;
    public ArrayList<Atomo> atomosReact;

    static HashMap<String, String[][]> keywords = palabrasClave();
    String bloque = "NONE";
    String ultimoBloqueGeom;

    ArrayList<String> opciones;

    String geome;
    String title;
    int charge = 0;
    int multiplicity;
    int n_elec;


    public DeMonParser() throws FileNotFoundException {
        br = new BufferedReader(new FileReader(fileName));   // Es el archivo de entrada cargado en memoria
        bloque = "NONE";
        opciones = new ArrayList<>();
        geome = "";
    }

    public DeMonParser(String fileName) throws FileNotFoundException {
        this.fileName = fileName;
        br = new BufferedReader(new FileReader(fileName));   // Es el archivo de entrada cargado en memoria
        bloque = "NONE";
        opciones = new ArrayList<>();
        geome = "";
    }

    public void finalize() throws Throwable {
        if(br != null)
            br.close();
        br = null;
        super.finalize();
    }




    public void parseGeometry() throws Exception {
        String linea;
        String nombreBloque;   //

        atomos = new ArrayList<>();
        atomosProd = new ArrayList<>();
        atomosReact = new ArrayList<>();
        while ( (linea = br.readLine()) != null) {

            if(linea.length() == 0 || esComentario(linea))
                continue;

            linea = linea.replaceAll("\\s+", " ").trim();
            ArrayList<String> tokens = new ArrayList<>( Arrays.asList(linea.split(" ")) );


			 /* ningun bloque activo */
            if(bloque.equals("NONE")) {
                nombreBloque = toKeyword(tokens.get(0));
                tokens.remove(0);
                if(!keywords.containsKey(nombreBloque)){
                    throw new Exception("Keyword " + nombreBloque + " no reconocida");
                }
                bloque = nombreBloque;
                procesaOpciones(bloque, tokens);
            }

			 /* termina bloque */
            else if( esFinBloque(toKeyword(tokens.get(0))) ) {
                if( bloque.equals("GEOME") || bloque.equals("PRODU") || bloque.equals("REACT") )
                    ultimoBloqueGeom = bloque;
                bloque = toKeyword(tokens.get(0));
                tokens.remove(0);
                procesaOpciones(bloque, tokens);
            }

			 /* cuerpo del bloque */
            else {
                switch(bloque) {
                    case "GEOME":
                        if(opciones.contains("CARTE")) {
                            Atomo atomo = new Atomo(Atomo.CARTESIAN, tokens);
                            atomos.add(atomo);
                            geome = "CARTE";
                        }
                        else if(opciones.contains("ZMATR") || opciones.contains("Z-MAT")) {
                            Atomo atomo = new Atomo(Atomo.ZMATRIX, tokens);
                            atomos.add(atomo);
                            geome = "ZMATR";
                        }
                        break;
                    case "CONST":
                        ArrayList<Atomo> listaAtomos = new ArrayList<>();
                        switch(ultimoBloqueGeom) {
                            case "GEOME":
                                listaAtomos = atomos;
                                break;
                            case "PRODU":
                            case "REACT":
                                listaAtomos.addAll(atomosReact);
                                listaAtomos.addAll(atomosProd);
                                break;
                        }
                        if(geome.equals("CARTE")) {
                            ArrayList<Atomo> atomosConst = Atomo.getAllByLabel(listaAtomos, tokens.get(0).toUpperCase());// System.out.println(atomosConst);
                            String coords = tokens.get(1);
                            int consts = 0;
                            for(int i = 0; i < coords.length(); i++) {
                                switch(coords.toUpperCase().charAt(i)) {
                                    case 'X':
                                        consts |= Atomo.CONST_X;
                                        break;
                                    case 'Y':
                                        consts |= Atomo.CONST_Y;
                                        break;
                                    case 'Z':
                                        consts |= Atomo.CONST_Z;
                                        break;
                                }
                            }
                            for(Atomo a: atomosConst) {
                                a.setConstantes(consts);
                            }
                        } else if(geome.equals("ZMATR")) {
                            String labelConst = tokens.get(0);
                            Double valorConst = Double.parseDouble(tokens.get(1));
                            for(Atomo a: listaAtomos) {
                                if(a.tipoCoords != Atomo.ZMATRIX)
                                    continue;
                                for(ZMatrixCoord zC: a.coordsZ) {
                                    if(zC.porAsignar != null && zC.porAsignar.equals(labelConst)) {
                                        zC.valNum = valorConst;
                                        zC.constLabel = labelConst;
                                    }
                                }
                            }
                        }
                        break;
                    case "VARIA":
                        ArrayList<Atomo> listaAtomosVar = new ArrayList<>();
                        switch(ultimoBloqueGeom) {
                            case "GEOME":
                                listaAtomosVar = atomos;
                                break;
                            case "PRODU":
                                listaAtomosVar = atomosProd;
                                break;
                            case "REACT":
                                listaAtomosVar = atomosReact;
                                break;
                        }
                        if(geome.equals("ZMATR")) {
                            String labelVar = tokens.get(0);
                            Double valorVar = Double.parseDouble(tokens.get(1));
                            for(Atomo a: listaAtomosVar) {
                                if(a.tipoCoords != Atomo.ZMATRIX)
                                    continue;
                                for(ZMatrixCoord zC: a.coordsZ) {
                                    if(zC.porAsignar != null && zC.porAsignar.equals(labelVar)) {
                                        zC.valNum = valorVar;
                                        zC.varLabel = labelVar;
                                    }
                                }
                            }
                        }
                        break;
                    case "PRODU":
                        if(opciones.contains("CARTE")) {
                            Atomo atomo = new Atomo(Atomo.CARTESIAN, tokens);
                            atomosProd.add(atomo);
                            geome = "CARTE";
                        }
                        else if(opciones.contains("ZMATR") || opciones.contains("Z-MAT")) {
                            Atomo atomo = new Atomo(Atomo.ZMATRIX, tokens);
                            atomosProd.add(atomo);
                            geome = "ZMATR";
                        }
                        break;
                    case "REACT":
                        if(opciones.contains("CARTE")) {
                            Atomo atomo = new Atomo(Atomo.CARTESIAN, tokens);
                            atomosReact.add(atomo);
                            geome = "CARTE";
                        }
                        else if(opciones.contains("ZMATR") || opciones.contains("Z-MAT")) {
                            Atomo atomo = new Atomo(Atomo.ZMATRIX, tokens);
                            atomosReact.add(atomo);
                            geome = "ZMATR";
                        }
                        break;
                    default:
                        throw new Exception("Bloque desconocido: " + bloque);
                }
            }
        }
        n_elec = 0;
        for(Atomo a: atomos) {
            n_elec += a.atomicN;
        }
        n_elec -= charge;
        Globals.multip = n_elec%2 == 0? 1: 2;

        Globals.N_ATOM = atomos.size();
    }

    private boolean esComentario(String linea) {
        return linea.startsWith("#");
    }





    private static HashMap<String, String[][]> palabrasClave() {
        HashMap<String, String[][]> palabras = new HashMap<>();

        palabras.put("TITLE", new String[0][0]);
        String [][] aux = { {"CARTE", "ZMATR", "Z-MAT", "MIXED"}, {"ANGST", "BOHR"} };
        palabras.put("GEOME", aux);
        palabras.put("PRODU", aux);
        palabras.put("REACT", aux);
        palabras.put("CONST", new String[0][0]);
        palabras.put("VARIA", new String[0][0]);
        palabras.put("MULTI", new String[0][0]);

        return palabras;
    }




    private void procesaOpciones(String bloque, ArrayList<String> tokens) throws Exception {
        if(!keywords.containsKey(bloque))
            throw new Exception("Bloque " + bloque + " no reconocido");
        switch(bloque) {
            case "TITLE":
                title = tokens.get(0);
                System.out.println("TITULO " + tokens);
                break;
            case "MULTI":
                multiplicity = Integer.parseInt(tokens.get(0));
                System.out.println("MULTIPLICITY " + tokens);
                break;
            case "GEOME":
//                System.out.println("GEOME");
                for(String token: tokens) {
                    String [][] opcionesBloque = keywords.get(bloque);
                    for(String[] opcionGrupo: opcionesBloque) {
                        if( Arrays.asList(opcionGrupo).contains(toKeyword(token)) )
                            opciones.add(toKeyword(token));
                    }
                }
                break;
            case "PRODU":
                for(String token: tokens) {
                    String [][] opcionesBloque = keywords.get(bloque);
                    for(String[] opcionGrupo: opcionesBloque) {
                        if( Arrays.asList(opcionGrupo).contains(toKeyword(token)) )
                            opciones.add(toKeyword(token));
                    }
                }
                break;
            case "REACT":
                for(String token: tokens) {
                    String [][] opcionesBloque = keywords.get(bloque);
                    for(String[] opcionGrupo: opcionesBloque) {
                        if( Arrays.asList(opcionGrupo).contains(toKeyword(token)) )
                            opciones.add(toKeyword(token));
                    }
                }
                break;
        }
    }





    private boolean esFinBloque(String token) {
        return esComentario(token) || token.equals("END") || keywords.containsKey(toKeyword(token));
    }




    /**
     * Formatea @param palabra como una keyword, una Keyword es maximo 5 letras y en mayúsculas.
     * @param palabra una supuesta palabra clave.
     * @return máximo las primeras 5 letras de la palabra, en mayúsculas.
     */
    private String toKeyword(String palabra) {
        return (palabra.length() < 5)? palabra.toUpperCase(): palabra.toUpperCase().substring(0, 5);
    }

}
