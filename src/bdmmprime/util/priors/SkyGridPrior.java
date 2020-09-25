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

    public Input<Function> sigmaInput = new Input<>("sigma",
            "Standard deviation of increment priors.", Input.Validate.REQUIRED);

    public SkyGridPrior() { }

    Function x, sigma;
    int n;

    final double logOneOnSqrt2Pi = -0.5*Math.log(2*Math.PI);

    @Override
    public void initAndValidate() {
        x = xInput.get();
        sigma = sigmaInput.get();
        n = x.getDimension();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        double logGausNorm = logOneOnSqrt2Pi - Math.log(sigma.getArrayValue());
        double sigma2 = sigma.getArrayValue()*sigma.getArrayValue();

        double prevEl = Math.log(x.getArrayValue(0));

        for (int i=1; i<n; i++) {
            double el = Math.log(x.getArrayValue(i));
            double delta = el - prevEl;

            logP += logGausNorm - 0.5*delta*delta/sigma2;
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
