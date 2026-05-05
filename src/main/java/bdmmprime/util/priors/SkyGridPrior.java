package bdmmprime.util.priors;

import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;

import java.util.List;
import java.util.Random;

public class SkyGridPrior extends Distribution {

    public Input<RealVector<?>> xInput = new Input<>("x",
            "Parameter to place prior on.", Input.Validate.REQUIRED);

    public Input<RealScalar<?>> MInput = new Input<>("M",
            "M parameter for log normal distribution of first element.",
            Input.Validate.REQUIRED);

    public Input<RealScalar<?>> SInput = new Input<>("S",
            "S parameter for log normal distribution of first element.",
            Input.Validate.REQUIRED);

    public Input<RealScalar<?>> sigmaInput = new Input<>("sigma",
            "Standard deviation of increment priors.", Input.Validate.REQUIRED);

    public SkyGridPrior() { }

    RealVector<?> x;
    int n;

    final double logOneOnSqrt2Pi = -0.5*Math.log(2*Math.PI);

    @Override
    public void initAndValidate() {
        x = xInput.get();
        n = x.size();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        double sigma = sigmaInput.get().get();
        double M = MInput.get().get();
        double S = SInput.get().get();

        double prevEl = Math.log(x.get(0));

        // Log normal distribution for initial element:
        logP += logOneOnSqrt2Pi - Math.log(S) - prevEl - 0.5*(prevEl - M)*(prevEl - M)/S/S;

        double logGausNorm = logOneOnSqrt2Pi - Math.log(sigma);
        double sigma2 = sigma*sigma;

        for (int i=1; i<n; i++) {
            double el = Math.log(x.get(i));
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
        throw new UnsupportedOperationException("Sampling from SkyGridPrior not supported.");
    }
}
