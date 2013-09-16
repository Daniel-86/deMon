package demon;


public class ZMatrixCoord {
	
	String label;
	Double valNum;
	String porAsignar;
	String constLabel;
	String varLabel;
	
	
	public ZMatrixCoord(String label, String val) {
		this.label = label;
		try {
			valNum = Double.parseDouble(val);
		} catch (Exception e) {
			porAsignar = val;
		}
	}
	
	
	
	public String toString() {
		return label + " " + valNum + " ";
	}

}
