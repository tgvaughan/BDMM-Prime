package bdmmprime.util.priors;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.IID;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.type.RealVector;
import beast.base.spec.type.Scalar;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class ZeroExcludingRealIID extends IID<RealVector<NonNegativeReal>, Scalar<Real, Double>, Double> {

    public Input<RealVector<? extends NonNegativeReal>> argInput = new Input<>("arg",
            "RealVector from which to exclude zeros.",
            Input.Validate.REQUIRED);

    RealVector<? extends Real> arg;

    List<Integer> indices;

    PositiveReal domain;

    @Override
    public void initAndValidate() {

        arg = argInput.get();

        // Set up index map
        indices = new ArrayList<>();
        for (int toIdx=0; toIdx<arg.size(); toIdx++) {
            if (arg.get(toIdx) != 0.0)
                indices.add(toIdx);
        }

        // Set up domain
        if (arg.getDomain() instanceof PositiveReal castedArgDomain)
            domain = castedArgDomain;
        else
            domain = PositiveReal.INSTANCE;
    }

    public int size() {
        return indices.size();
    }

    public double get(int i) {
        return arg.get(indices.get(i));
    }

    public List<Double> getElements() {
        List<Double> vals = new ArrayList<>();
        for (int i=0; i<size(); i++)
            vals.add(get(i));
        return vals;
    }

    @Override
    public void init(PrintStream out) {

    }

    @Override
    public void log(long sample, PrintStream out) {

    }

    @Override
    public void close(PrintStream out) {

    }
}
