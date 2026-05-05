package bdmmprime.util.priors;

import beast.base.core.Input;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.IID;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.type.RealVector;
import beast.base.spec.type.Scalar;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SmartZeroExcludingRealIID extends IID<RealVector<NonNegativeReal>, Scalar<Real, Double>, Double> {

    public Input<RealVector<? extends Real>> classesToExcludeInput = new Input<>("classesToExclude",
            "Elements with these value will be excluded from the vector.");

    List<Integer> indices;

    public SmartZeroExcludingRealIID() { super(); }

    public SmartZeroExcludingRealIID(RealVector<? extends NonNegativeReal> param, ScalarDistribution<?, Double> dist) {
        initByName("param", param, "distr", dist);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // Making the classesToExclude values already "seen" causes them not
        // to be added to the index list:
        Set<Double> seenValues = new HashSet<>();
        if (classesToExcludeInput.get() != null)
            seenValues.addAll(classesToExcludeInput.get().getElements());

        // Set up index map
        indices = new ArrayList<>();
        for (int toIdx=0; toIdx<param.size(); toIdx++) {
            double argVal = param.get(toIdx);
            if (argVal != 0.0 && !seenValues.contains(argVal)) {
                indices.add(toIdx);
                seenValues.add(argVal);
            }
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
