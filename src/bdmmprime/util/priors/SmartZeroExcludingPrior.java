package bdmmprime.util.priors;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.RealParameter;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@Description("Like ZeroExcludingPrior, but only applies the prior to unique " +
        "nonzero elements of the parameter")
public class SmartZeroExcludingPrior extends Prior {

    public Input<Function> classesToExcludeInput = new Input<>("classesToExclude",
            "Elements having this value will be excluded from the prior calculation");

    List<Integer> indices;

    @Override
    public void initAndValidate() {

        indices = new ArrayList<>();

        // Making the classesToExclude values already "seen" causes them not
        // to be added to the index list:
        Set<Double> seenValues = new HashSet<>();
        if (classesToExcludeInput.get() != null) {
            for (double value : classesToExcludeInput.get().getDoubleValues())
                seenValues.add(value);
        }

        for (int i=0; i<m_x.get().getDimension(); i++) {
            double thisValue = m_x.get().getArrayValue(i);
            if (thisValue != 0.0 && !seenValues.contains(thisValue)) {
                indices.add(i);
                seenValues.add(thisValue);
            }
        }

        super.initAndValidate();
    }

    public double calculateLogP() {
        Function x = m_x.get();
        if (x instanceof RealParameter) {
            // test that parameter is inside its bounds
            double l = ((RealParameter) x).getLower();
            double h = ((RealParameter) x).getUpper();
            for (int i : indices) {
                double value = x.getArrayValue(i);
                if (value < l || value > h) {
                    logP = Double.NEGATIVE_INFINITY;
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }

        logP = 0.0;
        for (int i : indices)
            logP += dist.logDensity(x.getArrayValue(i));

        if (logP == Double.POSITIVE_INFINITY) {
            logP = Double.NEGATIVE_INFINITY;
        }

        return logP;
    }
}
