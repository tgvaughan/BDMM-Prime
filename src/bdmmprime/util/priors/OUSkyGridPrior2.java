package bdmmprime.util.priors;

import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.State;

import java.util.*;

public class OUSkyGridPrior2 extends Distribution {

    public Input<Function> xInput = new Input<>("x",
            "Parameter to place prior on.", Input.Validate.REQUIRED);

    public Input<Function> MInput = new Input<>("M",
            "M parameter for log-normal distribution.", Input.Validate.REQUIRED);

    public Input<Function> SInput = new Input<>("S",
            "S parameter for log-normal distribution.", Input.Validate.REQUIRED);

    public Input<Boolean> meanInRealSpaceInput = new Input<>("meanInRealSpace",
            "Same as meanInRealSpace input for log normal distribution.",
            false);

    public Input<Function> thetaInput = new Input<>("theta",
            "Relaxation parameter for O-U process.", Input.Validate.REQUIRED);

    public Input<List<Double>> classesToExcludeInput = new Input<>("classToExclude",
            "Elements having this value will be excluded from the prior calculation",
            new ArrayList<>());

    public OUSkyGridPrior2() { }

    Function x;
    int n;
    List<Integer> indices;

    final double logOneOnSqrt2Pi = -0.5*Math.log(2*Math.PI);

    @Override
    public void initAndValidate() {
        x = xInput.get();
        n = x.getDimension();

        indices = new ArrayList<>();

        // Making the classesToExclude values already "seen" causes them not
        // to be added to the index list:
        Set<Double> seenValues = new HashSet<>(classesToExcludeInput.get());

        for (int i=0; i<n; i++) {
            double thisValue = x.getArrayValue(i);
            if (thisValue != 0.0 && !seenValues.contains(thisValue)) {
                indices.add(i);
                seenValues.add(thisValue);
            }
        }
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        // Parameters for O-U process:
        double M = MInput.get().getArrayValue();
        double S = SInput.get().getArrayValue();
        double theta = thetaInput.get().getArrayValue();

        if (meanInRealSpaceInput.get())
            M = Math.log(M);

        double expNegTheta = Math.exp(-theta);
        double S2 = S*S;

        double var = S2*(1.0 - expNegTheta*expNegTheta);
        double logGausNorm = logOneOnSqrt2Pi - 0.5*Math.log(var);

        // Keep track of previous value:
        double prevEl = Math.log(x.getArrayValue(indices.get(0)));

        // Log normal distribution for initial value
        logP +=  logOneOnSqrt2Pi - Math.log(S) - 0.5*(prevEl-M)*(prevEl-M)/S2 - prevEl;

        for (int i : indices) {
            double el = Math.log(x.getArrayValue(i));
            double mean = prevEl*expNegTheta + M*(1.0 - expNegTheta);
            double delta = el - mean;

            logP += logGausNorm - 0.5*delta*delta/var - el;

            prevEl = el;
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
        throw new UnsupportedOperationException("Sampling from OUSkygridPrior not supported.");
    }
}
