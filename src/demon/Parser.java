package demon;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Parsea el archivo de entrada y crea las listas de átomos.
 */
public class Parser {

	static HashMap<String, String[][]> keywords = palabrasClave();
	String archivo;
	String bloque;
	String ultimoBloqueGeom;
	ArrayList<String> opciones;
	public ArrayList<Atomo> atomos;
	public ArrayList<Atomo> atomosProd;
	public ArrayList<Atomo> atomosReact;
	String geome;
    BufferedWriter bw;
    String title;
    int charge = 0;
    int multiplicity;
    int n_elec;
    String matdia = "DSYEV";
    String startGuess = "HUECKEL";
    String scfTheory = "DEFAULT";
    double scfTol = 1e-5;
    boolean scfAdj = false;
    boolean cdxc = true;
    String fileToWrite = "MOLDEN";
    int scfMax = 100;
    boolean diis = true;
    boolean port = true;
    boolean l_shift = false;
    DecimalFormat coordsFormat = new DecimalFormat("###0.000000");
    DecimalFormat massFormat = new DecimalFormat("####.###");
    String trj_type = "NEW";
    String whichGeom = "NORMAL";



    /**
     * Parsea el archivo .inp, crea lista para los átomos; los principales, los reactantes y los productos.
     * @param ruta ubicación del archivo de entrada. Si no se especifica alguno, toma por defecto ./deMon.inp
     */
	public Parser(String ruta) {
		archivo = (ruta != null)? ruta:"deMon.inp" ;
		bloque = "NONE";
		opciones = new ArrayList<>();
		geome = "";
    }


    public void inputResume() throws IOException {
        bw.write(" " + Parser.repeat("*", 79) + "\n\n");
        bw.write(" " + title + "\n\n");
        bw.write(" " + Parser.repeat("*", 79) + "\n\n");

        bw.write("\n *** MOLECULE DATA ***\n");
        bw.write("\n MOLECULAR CHARGE: " + charge );
        bw.write("\n MOLECULAR MULTIPLICITY: " + multiplicity );
        bw.write("\n NUMBER OF ATOMS: " + atomos.size() );
        bw.write("\n NUMBER OF ELECTRONS: " + n_elec + "\n\n");

        bw.write(" *** SCF OPTIONS ***\n");
//        if(true)
//            bw.write("\n SCF THEORY: RESTRICTED KOHN-SHAM (RKS)");
        switch(scfTheory) {
            case "NONE":
                bw.write("\n SCF THEORY: NOT USED");
                break;
            case "GUESSONLY":
                bw.write("\n SCF THEORY: GUESS ONLY CALCULATION");
                break;
            case "SCFMAX":
                bw.write("\n SCF THEORY: ENERGY ONLY CALCULATION (NO SCF)");
                break;
            case "OKS_ROKS":
                bw.write("\n SCF THEORY: RESTRICTED OPEN KOHN-SHAM (ROKS)");
                break;
            case "OKS_UKS":
                bw.write("\n SCF THEORY: UNRESTRICTED KOHN-SHAM (UKS)");
                break;
            default:
                bw.write("\n SCF THEORY: RESTRICTED KOHN-SHAM (RKS)");
                break;
        }
        bw.write("\n DIAGONALIZER: " + matdia );
        switch(startGuess) {
            case "HUECKEL":
                bw.write("\n START GUESS: HUECKEL DENSITY");
                break;
            case "FERMI":
                bw.write("\n START GUESS: FERMI DENSITY");
                break;
            case "PROJCT":
                bw.write("\n START GUESS: PROJECTION");
                break;
            case "RESTART":
                bw.write("\n START GUESS: RESTART DENSITY");
                break;
            default:
                bw.write("\n START GUESS: CORE DENSITY");
                break;
        }
        bw.write("\n SCF TOLERANCE: " + scfTol );
        if(scfAdj)
            bw.write("\n SCF TOLERANCE WILL BE TIGHTEN");
        bw.write("\n POTENTIAL DENSITY: " + (cdxc? "FITTED DENSITY": "ORBITAL DENSITY"));
        bw.write("\n NUMBER OF SCF ITERATIONS: " + scfMax );
        bw.write("\n DIIS ACCELERATION: " + (diis? "ON": "OFF"));
        if(l_shift)
            bw.write("\n AQUI VA ALGO");


        if(trj_type.equals("NEW") || trj_type.equals("RESTART")) {
            bw.write("\n\n *** GEOMETRY ***\n");
            switch(whichGeom) {
                case "ORIGINAL_O":
                    bw.write("\n FILE deMon.mol WILL BE WRITTEN");
                    break;
                case "ORIGINAL_PES":
                    bw.write("\n FILE deMon.mkl WILL BE WRITTEN");
                    break;
                case "RESTART_OPT_STEP":
                    bw.write("\n FILE deMon.wfn WILL BE WRITTEN");
                    break;
                default:
                    int asfd = 0;
                    bw.write("\n INPUT ORIENTATION in BOHR\n");
                    String noH = "  No."; // 6u
                    String atomH = "ATOM"; // 8i
                    String xH = "X"; // 26u
                    String yH = "Y"; // 41u
                    String zH = "Z"; // 56u
                    String zAtomH = "Z-ATOM"; // 65
                    String massH = "MASS"; // 74u
                    String optxyzH = "OPTXYZ"; // 80u
                    bw.write("\n  NO.  ATOM        X              Y              Z        Z-ATOM   MASS   OPTXYZ\n");
                    int indxAtom = 0;
                    for(Atomo a: atomos) {
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
                    break;
            }
        }
        bw.write("\n\n");
    }

    public String symmetryResume() {

//        String resume = "";
//        String resumeTitle;
//
//        if(deMon.options.contains("OPT", "ZBUILD") && deMon.coordsType != deMon.COORDS_TYPE_CARTESIAN)
//            resumeTitle =  "THE GENERATED Z-MATRIX ORIENTATION IS USED";
//        else if(deMon.options.contains("OPT") && deMon.coordsType != deMon.COORDS_TYPE_CARTESIAN)
//            resumeTitle =  "THE Z-MATRIX ORIENTATION IS USED";
//        else
//            resumeTitle =  "THE INPUT ORIENTATION IS USED";
//
//        resume += deMon.printer.printSectionTitle(resumeTitle);
//        resume += deMon.printer.printPointGroup();
//        resume += deMon.printer.printGroupOrder();
//        resume += deMon.printer.printCoords();
        return "";
    }


    /**
     * Este es el método que parsea el archivo. El archivo está dividido en bloques
     * @throws Exception
     */
	public void parseGeom() throws Exception {
		
		 BufferedReader br = new BufferedReader(new FileReader(archivo));   // Es el archivo de entrada cargado en memoria
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
		 br.close();
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



    public void copyright() throws IOException {
        bw = new BufferedWriter(new FileWriter("deMon.out"));   // Es el archivo de entrada cargado en memoria
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
        bw.write(copyright);
//        bw.close();
    }

    /**
     * función para repetir una cadena n veces (java no tiene nada parecido).
     * @param s cadena
     * @param n número de veces a repetir.
     * @return la cadena concatenada n veces con sigo misma.
     */
    public static String repeat(String s, int n) {//System.out.println("ASDF: "+s+", "+n);
        StringBuilder sb = new StringBuilder(s.length() * n);
        for (int i = 0; i < n; i++)
            sb.append(s);
        return sb.toString();
    }



    public void closeFiles() throws IOException {
        bw.close();
    }



    public BufferedWriter getOut() {
        return bw;
    }


}
