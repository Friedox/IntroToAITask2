import java.util.Arrays;
import java.util.Map;
import java.util.Scanner;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

class Matrix {
    private int nCol;
    private int nRow;
    protected double[][] matrix;

    public Matrix(int nRow, int nCol) {
        this.nCol = nCol;
        this.nRow = nRow;
        this.matrix = new double[nRow][nCol];
    }

    public Matrix(double[][] matrix) {
        this.nCol = matrix[0].length;
        this.nRow = matrix.length;
        this.matrix = matrix;
    }

    public Matrix(double[] vector) {
        this.nCol = 1;
        this.nRow = vector.length;
        this.matrix = new double[nRow][nCol];
        for(int i = 0; i < nRow; i++) {
            this.matrix[i][0] = vector[i];
        }
    }

    public int getnCol() {
        return nCol;
    }

    public int getnRow() {
        return nRow;
    }

    public void multiplyByScalar(double scalar) {
        for(int i = 0; i < nRow; i++) {
            for(int j = 0; j < nCol; j++) {
                this.matrix[i][j] *= scalar;
            }
        }
    }

    public Matrix add(Matrix other) {
        if(this.nCol != other.getnCol() || this.nRow != other.getnRow()) {
            throw new IllegalArgumentException("Matrices must have the same dimensions");
        }
        for(int i = 0; i < nRow; i++) {
            for(int j = 0; j < nCol; j++) {
                this.matrix[i][j] += other.matrix[i][j];
            }
        }

        return new Matrix(this.matrix);
    }


    public Matrix multiply(Matrix other) {
        if (this.nCol != other.getnRow()) {
            throw new IllegalArgumentException("Matrices must have the same dimensions");
        }

        Matrix product = new Matrix(this.nRow, other.getnCol());

        for(int i = 0; i < this.nRow; i++) {
            for(int j = 0; j < other.getnCol(); j++) {
                for(int k = 0; k < this.nCol; k++) {
                    product.matrix[i][j] += this.matrix[i][k]*other.matrix[k][j];
                }
            }
        }
        return product;
    }

    public Matrix inverse() {
        if (this.nCol != this.nRow) {
            throw new IllegalArgumentException("Matrix must be square");
        }
        double determinant = this.determinant();
        if (determinant == 0) {
            throw new IllegalArgumentException("Matrix must be invertible");
        }


        Matrix inverse = new Matrix(this.nRow, this.nCol);
        for(int i = 0; i < this.nRow; i++) {
            for(int j = 0; j < this.nCol; j++) {
                inverse.matrix[i][j] = Math.pow(-1, i+j)*this.subMatrix(i, j).determinant()/determinant;
            }
        }

        return inverse.transpose();

    }

    public double determinant() {
        if (this.nCol != this.nRow) {
            throw new IllegalArgumentException("Matrix must be square");
        }
        double determinant = 0;
        if (this.nCol == 1) {
            determinant = this.matrix[0][0];
        } else if (this.nCol == 2) {
            determinant = this.matrix[0][0]*this.matrix[1][1] - this.matrix[0][1]*this.matrix[1][0];
        } else {
            for (int i = 0; i < this.nCol; i++) {
                determinant += Math.pow(-1, i)*this.matrix[0][i]*this.subMatrix(0, i).determinant();
            }
        }
        return determinant;
    }

    public Matrix subMatrix(int i, int j) {
        if (i > nRow || j > nCol) {
            throw new IllegalArgumentException("Index out of bounds");
        }
        Matrix subMatrix = new Matrix(this.nRow - 1, this.nCol - 1);
        int row = 0;
        int col = 0;
        for (int k = 0; k < this.nRow; k++) {
            col = 0;
            if (k != i) {
                for (int l = 0; l < this.nCol; l++) {
                    if (l != j) {
                        subMatrix.matrix[row][col] = this.matrix[k][l];
                        col++;
                    }
                }
                row++;
            }
        }
        return subMatrix;
    }

    public String toString() {
        String s = "";
        for (int i = 0; i < this.nRow; i++) {
            s += "[";
            for (int j = 0; j < this.nCol; j++) {
                s += this.matrix[i][j];
                if (j != this.nCol - 1) {
                    s += ", ";
                }
            }
            s += "]\n";
        }
        return s;
    }

    public Matrix subtract(Matrix cN) {
        if(this.nCol != cN.getnCol() || this.nRow != cN.getnRow()) {
            throw new IllegalArgumentException("Matrices must have the same dimensions");
        }
        Matrix difference = new Matrix(this.nRow, this.nCol);
        for(int i = 0; i < nRow; i++) {
            for(int j = 0; j < nCol; j++) {
                difference.matrix[i][j] = this.matrix[i][j] - cN.matrix[i][j];
            }
        }
        return difference;
    }

    public int rank() {
        int rank = 0;
        Matrix reducedRowEchelonForm = this.reducedRowEchelonForm();
        for(int i = 0; i < this.nRow; i++) {
            if(reducedRowEchelonForm.matrix[i][i] != 0) {
                rank++;
                break;
            }
        }
        return rank;
    }

    public Matrix transpose() {
        double[][] transpose = new double[this.nCol][this.nRow];
        for(int i = 0; i < this.nRow; i++) {
            for(int j = 0; j < this.nCol; j++) {
                transpose[j][i] = this.matrix[i][j];
            }
        }
        return new Matrix(transpose);
    }
    public Matrix reducedRowEchelonForm() {

        Matrix reducedRowEchelonForm = new Matrix(this.matrix);
        int pivot = 0;
        for(int i = 0; i < this.nRow; i++) {
            if(pivot >= this.nCol) {
                break;
            }
            if(reducedRowEchelonForm.matrix[i][pivot] == 0) {
                for(int j = i + 1; j < this.nRow; j++) {
                    if(reducedRowEchelonForm.matrix[j][pivot] != 0) {
                        reducedRowEchelonForm.swapRows(i, j);
                        break;
                    }
                }
            }
            if(reducedRowEchelonForm.matrix[i][pivot] != 0) {
                reducedRowEchelonForm.multiplyRow(i, 1/reducedRowEchelonForm.matrix[i][pivot]);
                for(int j = 0; j < this.nRow; j++) {
                    if(i != j && reducedRowEchelonForm.matrix[j][pivot] != 0) {
                        reducedRowEchelonForm.addRow(i, j, -reducedRowEchelonForm.matrix[j][pivot]);
                    }
                }
                pivot++;
            }
        }
        return reducedRowEchelonForm;
    }

    public void swapRows(int i, int j) {
        double[] temp = this.matrix[i];
        this.matrix[i] = this.matrix[j];
        this.matrix[j] = temp;
    }

    public void multiplyRow(int i, double scalar) {
        for(int j = 0; j < this.nCol; j++) {
            this.matrix[i][j] *= scalar;
        }
    }

    public static Matrix identityMatrix(int n) {
        Matrix identityMatrix = new Matrix(n, n);
        for(int i = 0; i < n; i++) {
            identityMatrix.matrix[i][i] = 1;
        }
        return identityMatrix;
    }

    public void addRow(int i, int j, double scalar) {
        for(int k = 0; k < this.nCol; k++) {
            this.matrix[j][k] += scalar*this.matrix[i][k];
        }
    }

    public double[] getColumn(int i) {
        double[] column = new double[this.nRow];
        for(int j = 0; j < this.nRow; j++) {
            column[j] = this.matrix[j][i];
        }
        return column;
    }

    public Matrix getRow(int i) {
        Matrix row = new Matrix(1, this.nCol);
        for(int j = 0; j < this.nCol; j++) {
            row.matrix[0][j] = this.matrix[i][j];
        }
        return row;
    }

    public void setRow(double[] row, int i) {
        this.matrix[i] = row;
    }

    public void setCol(double[] col, int i) {
        for(int j = 0; j < this.nRow; j++) {
            this.matrix[j][i] = col[j];
        }
    }

    public static double distance(Matrix a, Matrix b){
        if (a.getnCol() != 1 || b.getnCol() != 1) {
            System.out.println("NOW WORKABLE");
        }
        Matrix c = a.subtract(b);

        double result = 0;

        for (int i = 0; i < a.getnCol(); i++) {
            result += c.matrix[1][i] * c.matrix[1][i];
        }

        return sqrt(result);
    }
}

public class Main {

    public static void InteriorPoint(int numRows, int numCols, Matrix C, Matrix A, double[] b, double accuracy, double alpha, Matrix X_zero) {
        System.out.println("STARTED");

        int counter = 0;
        while (true) {
            if (counter >= 10000) {
                System.out.println("The problem does not have solution!");
                break;
            }

            Matrix D = Matrix.identityMatrix(numCols);


            for (int i = 0; i < numCols; i++) {
                D.multiplyRow(i, X_zero.matrix[i][0]);
            }

            Matrix A_hat = A.multiply(D);

            Matrix C_hat = D.multiply(C);

            Matrix P = Matrix.identityMatrix(numCols).subtract(A_hat.transpose()
                    .multiply(A_hat.multiply(A_hat.transpose()).inverse())
                    .multiply(A_hat));

            Matrix C_p = P.multiply(C_hat);

            Matrix I_vector = new Matrix(numCols, 1);

            for (int i = 0; i < numCols; i++) {
                I_vector.matrix[i][0] = 1;
            }

            double minC_p = Double.MAX_VALUE;

            for (int i = 0; i < numCols; i++) {
                if (C_p.matrix[i][0] < minC_p) {
                    minC_p = C_p.matrix[i][0];
                }
            }

            double coeff = alpha/abs(minC_p);

            C_p.multiplyByScalar(coeff);

            Matrix X_hat = I_vector.add(C_p);

            Matrix X_new = new Matrix(numCols, 1);

            for (int i = 0; i < X_new.getnRow(); i++) {
                X_new.matrix[i][0] = X_hat.matrix[i][0]*D.matrix[i][i];
            }

            if (Matrix.distance(X_new, X_zero) < accuracy) {
                System.out.printf("Results with alpha = %f\n", alpha);
                System.out.printf("Iterations: %d\n", counter);
                System.out.println(X_new);
                break;
            }

            System.out.printf("Iteration: %d\n", counter);
            System.out.println(X_new);
            X_zero = X_new;

            counter++;
        }
    }


    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        double alpha1 = 0.5;

        double alpha2 = 0.9;

        System.out.print("Enter the number of constraints: ");
        int numRows = scanner.nextInt();
        System.out.print("Enter the number of variables: ");
        int numCols = scanner.nextInt();

        double[] cTemp = new double[numCols];
        System.out.println("Enter the coefficients of the objective function vector C:");
        for (int j = 0; j < numCols; j++) {
            cTemp[j] = scanner.nextDouble();
        }

        Matrix C = new Matrix(cTemp);

        double[][] aTemp = new double[numRows][numCols];
        System.out.println("Enter the coefficients of the constraint matrix A:");
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                aTemp[i][j] = scanner.nextDouble();
            }
        }

        Matrix A = new Matrix(aTemp);

        double[] b = new double[numRows];
        System.out.println("Enter the right-hand side vector b:");
        for (int i = 0; i < numRows; i++) {
            b[i] = scanner.nextDouble();
        }

        System.out.print("Enter the approximation accuracy: ");
        double accuracy = Double.parseDouble(scanner.next());

        double[] X_zero_temp = new double[numCols];
        System.out.println("Enter the coefficients of the X zero:");
        for (int j = 0; j < numCols; j++) {
            X_zero_temp[j] = scanner.nextDouble();
        }

        Matrix X_zero = new Matrix(X_zero_temp);

        boolean allZeros = true;
        for (double element : cTemp) {
            if (element != 0.0) {
                allZeros = false;
                break;
            }
        }

        boolean foundNegative = false;
        for (double element : b) {
            if (element < 0.0) {
                foundNegative = true;
                break;
            }
        }

        if (allZeros || foundNegative) {
            System.out.println("The method is not applicable!");
        } else {
            InteriorPoint(numRows, numCols, C, A, b, accuracy, 0.5, X_zero);

            InteriorPoint(numRows, numCols, C, A, b, accuracy, 0.9, X_zero);
        }
    }
}