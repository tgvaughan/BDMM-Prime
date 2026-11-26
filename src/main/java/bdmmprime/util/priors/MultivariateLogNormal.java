/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.util.priors;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealVector;
import org.apache.commons.math4.legacy.linear.Array2DRowRealMatrix;
import org.apache.commons.math4.legacy.linear.EigenDecomposition;
import org.apache.commons.math4.legacy.linear.LUDecomposition;
import org.apache.commons.math4.legacy.linear.RealMatrix;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class MultivariateLogNormal extends Distribution {

    public Input<List<RealVector<?>>> xsInput = new Input<>("x",
            "RealVector to which distribution is applied",
            new ArrayList<>());

    public Input<RealVector<?>> MInput = new Input<>("M",
            "M parameter of the multivariate lognormal distribution.",
            Input.Validate.REQUIRED);

    public Input<RealVector<?>> SInput = new Input<>("S",
            "S parameter of the multivariate lognormal distribution.",
            Input.Validate.REQUIRED);

    List<RealVector<?>> xs;
    RealVector<?> M, S;
    int n, m;

    double[] y, mu;
    RealMatrix sigma;

    /** absolute tolerance of the symmetry check on the covariance matrix */
    static final double SYMMETRY_TOLERANCE = 1e-9;

    public MultivariateLogNormal() { }

    @Override
    public void initAndValidate() {
        xs = xsInput.get();
        M = MInput.get();
        S = SInput.get();

        n = xs.size();

        if (n==0)
            throw new IllegalArgumentException("Need at least one x input.");

        m = xs.getFirst().size();

        for (int i=1; i<n; i++) {
            if (xs.get(i).size() != m)
                throw new IllegalArgumentException("Each x must have the same dimension.");
        }

        if (M.size() != n)
            throw new IllegalArgumentException("M must have n elements, " +
                    "where n is the number of x inputs.");

        if (S.size() != n*n)
            throw new IllegalArgumentException("S must have n*n elements, " +
                    "where n is the number of x inputs.");

        y = new double[n];
        mu = new double[n];
        sigma = new Array2DRowRealMatrix(n, n);
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        for (int j=0; j<n; j++) {
            mu[j] = M.get(j);
            for (int k=0; k<n; k++) {
                sigma.setEntry(j, k, S.get(j*n + k));
            }
        }

        // Check that sigma is symmetric (required for covariance matrices):
        for (int j=0; j<n; j++)
            for (int k=j+1; k<n; k++)
                if (!(Math.abs(sigma.getEntry(j, k) - sigma.getEntry(k, j)) <= SYMMETRY_TOLERANCE)) // also rejects NaN
                    throw new IllegalArgumentException("Covariance matrix is not symmetric.");

        // Check that sigma is positive definite (required for covariance matrices):
        double[] evals = new EigenDecomposition(sigma).getRealEigenvalues();
        for (double eval : evals)
            if (eval < 0.0)
                return Double.NEGATIVE_INFINITY;

        // only an exactly zero pivot counts as singular, not the library's default of 1e-11
        LUDecomposition lu = new LUDecomposition(sigma, Double.MIN_VALUE);
        if (!lu.getSolver().isNonSingular())
            throw new IllegalArgumentException("Matrix is singular.");
        RealMatrix sigmaInv = lu.getSolver().getInverse();
        double sigmaDet = lu.getDeterminant();

        for (int i=0; i<m; i++) {
            double ySum = 0.0;
            for (int j=0; j<n; j++) {
                y[j] = Math.log(xs.get(j).get(i)) - M.get(j);
                ySum += y[j];
            }

            // y^T sigmaInv y
            double[] yTSigmaInv = sigmaInv.preMultiply(y);
            double quadForm = 0.0;
            for (int j=0; j<n; j++)
                quadForm += yTSigmaInv[j] * y[j];

            logP += -quadForm
                    - 0.5*n*Math.log(Math.PI)
                    - 0.5*sigmaDet
                    - ySum; // results from change of variables
        }

        return logP;
    }


    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Sampling from MultiVariateLogNormal is not implemented.");
    }

    /**
     * Main method for debugging only
     *
     * @param args unused
     */
    public static void main(String[] args) {

        RealVectorParam<Real> x1Param = new RealVectorParam<>(new double[] {1.0}, Real.INSTANCE);
        RealVectorParam<Real> x2Param = new RealVectorParam<>(new double[] {1.0}, Real.INSTANCE);

        MultivariateLogNormal mvl = new MultivariateLogNormal();
        mvl.initByName(
                "x", x1Param, "x", x2Param,
                "M", new RealVectorParam<>(new double[]{0.1, 0.2}, Real.INSTANCE),
                "S", new RealVectorParam<>(new double[]{1, 0.1, 0.1, 1}, Real.INSTANCE)
        );

        System.out.println(mvl.calculateLogP());
    }
}
