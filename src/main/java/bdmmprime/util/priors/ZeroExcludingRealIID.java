package bdmmprime.util.priors;

import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.IID;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.type.RealVector;
import beast.base.spec.type.Scalar;

import java.util.ArrayList;
import java.util.List;

public class ZeroExcludingRealIID extends IID<RealVector<NonNegativeReal>, Scalar<Real, Double>, Double> {

    List<Integer> indices;

    public ZeroExcludingRealIID() { super(); }

    public ZeroExcludingRealIID(RealVector<? extends NonNegativeReal> param, ScalarDistribution<?, Double> dist) {
        initByName("param", param, "dist", dist);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // Set up index map
        indices = new ArrayList<>();
        for (int toIdx=0; toIdx<param.size(); toIdx++) {
            if (param.get(toIdx) != 0.0)
                indices.add(toIdx);
        }
    }

    @Override
    protected double calcLogP(Double... value) {
        refresh(); // this make sure distribution parameters are updated if they are sampled during MCMC

        if (value == null)
            throw new IllegalArgumentException("IID requires param, but it is null ! ");
        if (value.length != dimension())
            throw new IllegalArgumentException("Values dimension != parameter dimension !");
        double logP = 0.0;
        for (int idx : indices) {
            logP += dist.density(value[idx]);
        }
        return logP;
    }

}
