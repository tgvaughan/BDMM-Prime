package bdmmprime.util.priors;

import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.State;

import java.util.List;
import java.util.Random;

public class SkyGridPrior extends Distribution {

    public Input<Function> xInput = new Input<>("x",
            "Parameter to place prior on.", Input.Validate.REQUIRED);

    public Input<Function> MInput = new Input<>("M",
            "M parameter for log normal distribution of first element.",
            Input.Validate.REQUIRED);

    public Input<Function> SInput = new Input<>("S",
            "S parameter for log normal distribution of first element.",
            Input.Validate.REQUIRED);

    public Input<Function> sigmaInput = new Input<>("sigma",
            "Standard deviation of increment priors.", Input.Validate.REQUIRED);

    public SkyGridPrior() { }

    Function x;
    int n;

    final double logOneOnSqrt2Pi = -0.5*Math.log(2*Math.PI);

    @Override
    public void initAndValidate() {
        x = xInput.get();
        n = x.getDimension();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        double sigma = sigmaInput.get().getArrayValue();
        double M = MInput.get().getArrayValue();
        double S = SInput.get().getArrayValue();

        double prevEl = Math.log(x.getArrayValue(0));

        // Log normal distribution for initial element:
        logP += logOneOnSqrt2Pi - Math.log(S) - prevEl - 0.5*(prevEl - M)*(prevEl - M)/S/S;

        double logGausNorm = logOneOnSqrt2Pi - Math.log(sigma);
        double sigma2 = sigma*sigma;

        for (int i=1; i<n; i++) {
            double el = Math.log(x.getArrayValue(i));
            double delta = el - prevEl;

            logP += logGausNorm - el - 0.5*delta*delta/sigma2;

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
        throw new UnsupportedOperationException("Sampling from SkygridPrior not supported.");
    }
}
