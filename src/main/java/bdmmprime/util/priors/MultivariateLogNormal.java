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
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.Property;

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

    DoubleMatrix2D y, mu, sigma;

    Algebra al = Algebra.DEFAULT;

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

        y = new DenseDoubleMatrix2D(n,1);
        mu = new DenseDoubleMatrix2D(n, 1);
        sigma = new DenseDoubleMatrix2D(n, n);
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        for (int j=0; j<n; j++) {
            mu.set(j, 0, M.get(j));
            for (int k=0; k<n; k++) {
                sigma.setQuick(j,k, S.get(j*n + k));
            }
        }

        // Check that sigma is symmetric (required for covariance matrices):
        if (!Property.DEFAULT.isSymmetric(sigma))
            throw new IllegalArgumentException("Covariance matrix is not symmetric.");

        // Check that sigma is positive definite (required for covariance matrices):
        EigenvalueDecomposition ed = new EigenvalueDecomposition(sigma);
        DoubleMatrix1D evals = ed.getRealEigenvalues();
        for (int i=0; i<evals.size(); i++)
            if (evals.getQuick(i) < 0.0)
                return Double.NEGATIVE_INFINITY;

        DoubleMatrix2D sigmaInv = Algebra.DEFAULT.inverse(sigma);

        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++)
                y.set(j, 0, Math.log(xs.get(j).get(i)) - M.get(j));

            logP += -al.mult(al.mult(al.transpose(y), sigmaInv), y).getQuick(0,0)
                    - 0.5*n*Math.log(Math.PI)
                    - 0.5* al.det(sigma)
                    - y.zSum(); // results from change of variables
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
