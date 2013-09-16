package demon;

import java.util.ArrayList;

public class Atomo {
	
	final static int CONST_X = 0x01;
	final static int CONST_Y = 0x02;
	final static int CONST_Z = 0x04;
	final static int CARTESIAN = 0x01;
	final static int ZMATRIX = 0x02;
    public static AtomDetails[] atomsDefaults = AtomDetails.buildAtomsDefaults();
	
	String label;
	int tipoCoords;
	ArrayList<Double> coords;
	ArrayList<ZMatrixCoord> coordsZ;
	Double carga;
	Double masa;
	int constante;
    int atomicN;
    Double ain;
	
	
	public Atomo(int tipoCoords, ArrayList<String> datos) throws Exception {
		this.tipoCoords = tipoCoords;
		switch(tipoCoords) {
			case CARTESIAN:
				if(datos.size() < 4)
					throw new Exception("Error GEOMETRIA: faltan datos");
				label = datos.get(0);
				coords = new ArrayList<>(3);
				for(int i = 1; i < 4; i++ ) {
					coords.add( Double.parseDouble(datos.get(i)) );
				}
				carga = (datos.size() >= 5)? Double.parseDouble(datos.get(4)): defaultCarga();
				masa = (datos.size() >= 6)? Double.parseDouble(datos.get(5)): defaultMasa();
                atomicN = (int) AtomDetails.getBySymbol(atomsDefaults, label.replaceAll("\\d+", ""), AtomDetails.ATMNUMBER);
				break;
			case ZMATRIX:
				label = datos.get(0);
				coordsZ = new ArrayList<>();
				for(int i = 1; i < datos.size(); i+=2) {
					coordsZ.add(new ZMatrixCoord(datos.get(i), datos.get(i+1)));
				}
                carga = defaultCarga();
                masa = defaultMasa();
                atomicN = (int) AtomDetails.getBySymbol(atomsDefaults, label.replaceAll("\\d+", ""), AtomDetails.ATMNUMBER);
				break;
		}
		constante = 0;
        ain = carga + masa;
	}
	
	
	
	private Double defaultCarga() {
		String simbolo = label.replaceAll("\\d+", "");
		//return cargas.get(simbolo);
        return (Double) AtomDetails.getBySymbol(atomsDefaults, simbolo, AtomDetails.COVALENT);
	}



	private Double defaultMasa() {
		String simbolo = label.replaceAll("\\d+", "");
		//return masas.get(simbolo);
        return (Double)AtomDetails.getBySymbol(atomsDefaults, simbolo, AtomDetails.MASS);
	}
	
	
	
	public String toString() {
        return (tipoCoords == Atomo.CARTESIAN)? label + " " + coords.toString() + " " + carga + " " + masa + " CONSTANTES: " + constante: label + " " + coordsZ.toString();
	}

	
	
	public static ArrayList<Atomo> getAllByLabel(ArrayList<Atomo> atomos, String label) {
		ArrayList<Atomo> atoms = new ArrayList<>();
		boolean porLabel = label.matches("[A-Z]+\\d+");// System.out.println("busca por label: " + porLabel + "\t" + label);
		for(Atomo atomo: atomos) {
			if(porLabel && atomo.label.equals(label)) {
				atoms.add(atomo);
			}
			if(!porLabel && atomo.label.startsWith(label))
				atoms.add(atomo);
		}
		return atoms;
	}
	
	
	
	public void setConstantes(int constant) {
		constante = constant;
	}


    public void shiftCenterMass(Double[] translationVector) {
        for(int i = 0; i < 3; i++)
            this.coords.set(i, this.coords.get(i) - translationVector[i]);
    }

    public void printCoords() {
        String coordsString= "";
        for(Double c: this.coords) {
            coordsString += c + " ";
        }
        System.out.println(coordsString);
    }
}
