package demon;

/**
 * User: daniel
 * Date: 20/05/13
 * Time: 09:04 PM
 */
public class MineVector {

    Object [] vector;

    public MineVector(int dimension) {
        vector = new Object[dimension];
    }

    public String toString() {
        String vectorString = "[ ";
        for(Object val: vector) {
            vectorString += val + " ";
        }
        vectorString += "]";
        return vectorString;
    }

    public MineVector fillWith(Object val) {
        for(int i = 0; i < vector.length; i ++) {
            vector[i] = val;
        }
        return this;
    }

}
