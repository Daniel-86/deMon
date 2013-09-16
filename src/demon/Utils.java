package demon;

/**
 * User: daniel
 * Date: 16/06/13
 * Time: 02:29 PM
 */
public class Utils {

    public static void fillVectorWith(Object[] vector, Object val) {
        for(int i = 0; i < vector.length; i ++)
            vector[i] = val;
    }


    public static void fillMatrixWith(Object[][] matrix, Object val) {
//        double[][] matrix
        int nRows = matrix.length;
        int nCols = matrix[0].length;
        for(int i = 0; i < nRows; i ++)
            for(int j = 0; j < nCols; j++)
                matrix[i][j] = val;
    }


    public static void printVector(Object[] vector) {
        for (Object aVector : vector) System.out.print(aVector + " ");
        System.out.println();
    }


    public static void printMatrix(Object[][] matrix) {
        int nRows = matrix.length;
        int nCols = matrix[0].length;
        for (Object[] aMatrix : matrix) {
            for (int j = 0; j < nCols; j++)
                System.out.print(aMatrix[j] + " ");
            System.out.println();
        }
    }

    public static Object [][] buildIdentityMatrix(int dimension, int type) {
        Object [][] matrix;
        Object val;
        switch(type) {
            case 0:
                matrix = new Integer [dimension][dimension];
                fillMatrixWith(matrix, 0);
                val = 1;
                break;
            case 1:
                matrix = new Double [dimension][dimension];
                fillMatrixWith(matrix, 0.0);
                val = 1.0;
                break;
            default:
                return null;
        }
        for(int i = 0; i < dimension; i ++)
            matrix[i][i] = val;

        return matrix;
    }
}
