package bdmmprime.util.priors;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealVector;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SmartZeroExcludingRealVector extends CalculationNode implements RealVector<PositiveReal>, Loggable {

    public Input<RealVector<? extends NonNegativeReal>> argInput = new Input<>("arg",
            "RealVector from which to exclude zeros.",
            Input.Validate.REQUIRED);

    public Input<RealVector<? extends Real>> classesToExcludeInput = new Input<>("classesToExclude",
            "Elements with these value will be excluded from the vector.");


    RealVector<? extends Real> arg;

    List<Integer> indices;

    PositiveReal domain;

    @Override
    public void initAndValidate() {

        arg = argInput.get();

        // Making the classesToExclude values already "seen" causes them not
        // to be added to the index list:
        Set<Double> seenValues = new HashSet<>();
        if (classesToExcludeInput.get() != null) {
            seenValues.addAll(classesToExcludeInput.get().getElements());
        }

        // Set up index map
        indices = new ArrayList<>();
        for (int toIdx=0; toIdx<arg.size(); toIdx++) {
            double argVal = arg.get(toIdx);
            if (argVal != 0.0 && !seenValues.contains(argVal)) {
                indices.add(toIdx);
                seenValues.add(argVal);
            }
        }

        // Set up domain
        if (arg.getDomain() instanceof PositiveReal castedArgDomain)
            domain = castedArgDomain;
        else
            domain = PositiveReal.INSTANCE;
    }

    @Override
    public int size() {
        return indices.size();
    }

    @Override
    public double get(int i) {
        return arg.get(indices.get(i));
    }

    @Override
    public List<Double> getElements() {
        List<Double> vals = new ArrayList<>();
        for (int i=0; i<size(); i++)
            vals.add(get(i));
        return vals;
    }

    @Override
    public PositiveReal getDomain() {
        return domain;
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
