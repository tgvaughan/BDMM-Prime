package bdmmprime.util;

import beast.core.Description;
import beast.core.Function;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.math.distributions.Prior;

@Description("Like Prior, but excludes any zero elements of the input parameter " +
        "from the calculation.")
public class ZeroExcludingPrior extends Prior {

    public double calculateLogP() {
        Function x = m_x.get();
        if (x instanceof RealParameter || x instanceof IntegerParameter) {
            // test that parameter is inside its bounds
            double l = 0.0;
            double h = 0.0;
            if (x instanceof RealParameter) {
                l = ((RealParameter) x).getLower();
                h = ((RealParameter) x).getUpper();
            } else {
                l = ((IntegerParameter) x).getLower();
                h = ((IntegerParameter) x).getUpper();
            }
            for (int i = 0; i < x.getDimension(); i++) {
                double value = x.getArrayValue(i);
                if (value < l || value > h) {
                    logP = Double.NEGATIVE_INFINITY;
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }

        logP = 0.0;
        for (double elementVal : x.getDoubleValues()) {
            if (elementVal == 0.0)
                continue;

            logP += dist.logDensity(elementVal);
        }

        if (logP == Double.POSITIVE_INFINITY) {
            logP = Double.NEGATIVE_INFINITY;
        }

        return logP;
    }
}
