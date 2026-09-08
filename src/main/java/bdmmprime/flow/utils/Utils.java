package bdmmflow.utils;

import com.flag4j.Matrix;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.linear.*;
import org.jblas.DoubleMatrix;

import java.util.Arrays;
import java.util.Random;

import org.hipparchus.linear.SchurTransformer;

public class Utils {

    /**
     * Returns a random square matrix of the given dimension.
     */
    public static RealMatrix getRandomMatrix(int dimension) {
        return getRandomMatrix(dimension, new Random().nextInt());
    }

    /**
     * Returns a random square matrix of the given dimension.
     */
    public static RealMatrix getRandomMatrix(int dimension, int seed) {
        RealMatrix randomMatrix = new BlockRealMatrix(dimension, dimension);
        Random random = new Random(seed);

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                randomMatrix.setEntry(i, j, random.nextGaussian());
            }
        }

        return randomMatrix;
    }

    public static RealMatrix computeRegularMinimizer(RealMatrix A) {
        // Compute SVD
        SingularValueDecomposition svd = new SingularValueDecomposition(A);
        RealMatrix V = svd.getV();
        double[] sigma = svd.getSingularValues();
        int n = sigma.length;

        // Compute product of singular values directly (avoid streams)
        double logProd = 0.0;
        for (double s : sigma) logProd += Math.log(s);
        double geomMean = Math.exp((Math.log(1.0) + logProd) / n); // d=1.0 inline

        // Compute V * Σ⁻¹ efficiently (scale columns of V by 1/σᵢ)
        double[][] vData = V.getData();
        for (int j = 0; j < n; j++) {
            double scale = geomMean / sigma[j];  // combine scalarMultiply(g_d) + Σ⁻¹
            for (int i = 0; i < n; i++) {
                vData[i][j] *= scale;
            }
        }
        RealMatrix Xstar = MatrixUtils.createRealMatrix(vData).multiply(V.transpose());

        return Xstar;
    }

    /**
     * Returns a RealMatrix matrix filled with the values in the array in column-major order.
     */
    public static RealMatrix toMatrix(double[] array, int n) {
        return Utils.toMatrix(array, n, 0);
    }

    /**
     * Returns a RealMatrix matrix filled with the values in the array in column-major order.
     */
    public static RealMatrix toMatrix(double[] array, int n, int offset) {
        RealMatrix matrix = new BlockRealMatrix(n, n);
        double[] column = new double[n];
        for (int j = 0; j < n; j++) {
            System.arraycopy(array, j * n + offset, column, 0, n);
            matrix.setColumn(j, column);
        }
        return matrix;
    }

    public static RealMatrix toMatrix(DoubleMatrix source) {
        RealMatrix destination = new BlockRealMatrix(source.rows, source.columns);
        for (int i = 0; i < source.rows; i++) {
            for (int j = 0; j < source.columns; j++) {
                destination.setEntry(i, j, source.get(i, j));
            }
        }
        return destination;
    }

    public static RealMatrix toMatrix(org.hipparchus.linear.RealMatrix source) {
        return new BlockRealMatrix(source.getData());
    }

    public static org.hipparchus.linear.RealMatrix toHipparchusMatrix(RealMatrix source) {
        return new org.hipparchus.linear.BlockRealMatrix(source.getData());
    }

    public static DoubleMatrix toMatrix(RealMatrix source) {
        DoubleMatrix destination = new DoubleMatrix(source.getRowDimension(), source.getColumnDimension());
        for (int i = 0; i < source.getRowDimension(); i++) {
            for (int j = 0; j < source.getColumnDimension(); j++) {
                destination.put(i, j, source.getEntry(i, j));
            }
        }
        return destination;
    }

    public static Matrix toFlag4JMatrix(RealMatrix source) {
        Matrix destination = new Matrix(source.getRowDimension(), source.getColumnDimension());
        for (int i = 0; i < source.getRowDimension(); i++) {
            for (int j = 0; j < source.getColumnDimension(); j++) {
                destination.set(source.getEntry(i, j), i, j);
            }
        }
        return destination;
    }

    public static RealMatrix toMatrix(Matrix source) {
        RealMatrix destination = new BlockRealMatrix(source.numRows(), source.numCols());
        for (int i = 0; i < source.numRows(); i++) {
            for (int j = 0; j < source.numCols(); j++) {
                destination.setEntry(i, j, source.get(i, j));
            }
        }
        return destination;
    }

    /**
     * Fills the given array with the values in the given matrix in column-major order.
     */
    public static void fillArray(RealMatrix matrix, double[] array) {
        Utils.fillArray(matrix, array, 0);
    }

    /**
     * Fills the given array with the values in the given matrix in column-major order.
     */
    public static void fillArray(RealMatrix matrix, double[] array, int offset) {
        int rows = matrix.getRowDimension();
        int cols = matrix.getColumnDimension();
        for (int j = 0; j < cols; j++) {
            System.arraycopy(matrix.getColumn(j), 0, array, j * rows + offset, rows);
        }
    }

    /**
     * Converts the given array to a RealVector vector.
     */
    public static RealVector toVector(double[] array) {
        return new ArrayRealVector(array);
    }

    /**
     * Scales the given array in-place such that the maximum value is 1.
     *
     * @param array the array to scale.
     * @return the log of the scaling factor applied.
     */
    public static double rescale(double[] array) {
        return Utils.rescale(array, 0);
    }

    /**
     * Scales the given array in-place such that the maximum value is 1. Assumes that
     * the array has already been scaled by the given previousLogFactor.
     *
     * @param array             the array to scale.
     * @param previousLogFactor the log of the previous scaling factor.
     * @return the log of the overall scaling factor applied.
     */
    public static double rescale(double[] array, double previousLogFactor) {
        double max = Math.abs(Arrays.stream(array).max().orElse(1.0));

        if (max == 0.0) return previousLogFactor;

        for (int i = 0; i < array.length; i++) {
            array[i] /= max;
        }

        return Math.log(max) + previousLogFactor;
    }

    /**
     * Scales the given vector in-place such that the maximum value is 1. Assumes that
     * the array has already been scaled by the given previousLogFactor.
     *
     * @param array             the array to scale.
     * @param previousLogFactor the log of the previous scaling factor.
     * @return the log of the overall scaling factor applied.
     */
    public static double rescale(RealVector array, double previousLogFactor) {
        double max = Math.abs(array.getMaxValue());
        if (max == 0.0) return previousLogFactor;

        array.mapMultiplyToSelf(1.0 / max);
        return Math.log(max) + previousLogFactor;
    }

    public static RealMatrix expm(RealMatrix A) {
        SchurTransformer schur = new SchurTransformer(Utils.toHipparchusMatrix(A));
        RealMatrix Q = Utils.toMatrix(schur.getP());
        RealMatrix T = Utils.toMatrix(schur.getT());

        RealMatrix expT = expmUpperTriangular(T);

        return Q.multiply(expT).multiply(Q.transpose());
    }

    public static RealMatrix expmUpperTriangular(RealMatrix T) {
        int n = T.getRowDimension();
        RealMatrix F = MatrixUtils.createRealMatrix(n, n);

        // Diagonal: exp(lambda)
        for (int i = 0; i < n; i++) {
            F.setEntry(i, i, Math.exp(T.getEntry(i, i)));
        }

        // Parlett recurrence
        for (int p = 1; p < n; p++) {
            for (int i = 0; i < n - p; i++) {
                int j = i + p;

                double sum = 0.0;
                for (int k = i + 1; k < j; k++) {
                    sum += T.getEntry(i, k) * F.getEntry(k, j)
                            - F.getEntry(i, k) * T.getEntry(k, j);
                }

                double tii = T.getEntry(i, i);
                double tjj = T.getEntry(j, j);

                if (Math.abs(tii - tjj) < 1e-12) {
                    // repeated eigenvalue
                    F.setEntry(i, j,
                            T.getEntry(i, j) * F.getEntry(i, i) + sum);
                } else {
                    F.setEntry(i, j,
                            (T.getEntry(i, j) * (F.getEntry(j, j) - F.getEntry(i, i)) + sum)
                                    / (tjj - tii));
                }
            }
        }

        return F;
    }

    public static double getHermitianSpread(RealMatrix matrix) {
        try {
            RealMatrix hermitian = matrix.add(matrix.transpose()).scalarMultiply(0.5);
            EigenDecomposition eigenDecomposition = new EigenDecomposition(hermitian);
            double minEV = Arrays.stream(eigenDecomposition.getRealEigenvalues()).min().orElseThrow();
            double maxEV = Arrays.stream(eigenDecomposition.getRealEigenvalues()).max().orElseThrow();
            return maxEV - minEV;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

}
