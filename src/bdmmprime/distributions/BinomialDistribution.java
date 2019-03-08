package bdmmprime.distributions;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

@Description("Binomial distribution, used as prior  Pr(k; n; p)=\\binom{n}{k} p^k (1-p)^{n-k}" +
        "If the input x is a multidimensional parameter, each of the dimensions is considered as a " +
        "separate independent component.")
public class BinomialDistribution extends ParametricDistribution {
    final public Input<RealParameter> pInput = new Input<>("p", "probability p parameter, defaults to 0.5");
    final public Input<RealParameter> trialsInput = new Input<>("trials", "number of trials parameter, defaults to 1");

    static org.apache.commons.math.distribution.BinomialDistribution dist = new BinomialDistributionImpl(1, 0.5);


    // Must provide empty constructor for construction by XML. Note that this constructor DOES NOT call initAndValidate();
    public BinomialDistribution() {
    }

    public BinomialDistribution(RealParameter p, RealParameter trials) {

        try {
            initByName("p", p);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("Failed to initByName p parameter when constructing BinomialDistribution instance.");
        }

        try {
            initByName("trials", trials);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("Failed to initByName trials parameter when constructing BinomialDistribution instance.");
        }
    }

    @Override
    public void initAndValidate() {
        refresh();
    }

    /**
     * make sure internal state is up to date *
     */
    @SuppressWarnings("deprecation")
    void refresh() {
        double p;
        if (pInput.get() == null) {
            p = 0.5;
        } else {
            p = pInput.get().getValue();
            if (p < 0) {
                p = 0;
            } else if (p > 1){
                p = 1;
            }
        }
        dist.setProbabilityOfSuccess(p);

        int trials;
        if (trialsInput.get() == null) {
            trials = 1;
        } else {
            trials = trialsInput.get().getValue().intValue();
            if (trials < 1) {
                trials = 1;
            }
        }
        dist.setNumberOfTrials(trials);
    }

    @Override
    public org.apache.commons.math.distribution.Distribution getDistribution() {
        refresh();
        return dist;
    }

}