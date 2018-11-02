package bdmm.distributions;

import beast.core.BEASTObject;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class Skyline extends CalculationNode {

    public Input<RealParameter> changeTimesInput = new Input<>("changeTimes",
            "Parameter containing change times for skyline function.");

    public Input<Boolean> timesAreRelativeInput = new Input<>("timesAreRelative",
            "True if times are relative to tree height. (Default false.)",
            false);

    public Input<Boolean> timesReversedInput = new Input<>("timesReversed",
            "True if times are reversed (ages instead of times after origin).",
            false);


    public Input<RealParameter> rateValuesInput = new Input<>("rateValues",
            "Parameter specifying rate values through time.",
            Input.Validate.REQUIRED);

    boolean reverse, relative;


    @Override
    public void initAndValidate() {
        reverse = timesReversedInput.get();
        relative = timesAreRelativeInput.get();
    }

    public boolean timesAreRelative() {
        return timesAreRelativeInput.get();
    }

    public boolean timesReversed() {
        return timesReversedInput.get();
    }

    List<Double> changeTimes = new ArrayList<>();

    public List<Double> getSanitizedChangeTimes(double maxTime) {
		changeTimes.clear();

        if (!reverse && changeTimesInput.get().getValue(0) != 0.0) {
            throw new RuntimeException("First time in interval times parameter should always be zero.");
        }

        int dim = changeTimesInput.get().getDimension();

        double end;
        for (int i = (reverse?0:1); i < dim; i++) {
            end = reverse ? (maxTime - changeTimesInput.get().getValue(dim - i - 1)) : changeTimesInput.get().getValue(i);
            if (relative) end *= maxTime;
            if (end < maxTime) changeTimes.add(end);
        }

        end = maxTime;

        changeTimes.add(end);

        Collections.sort(changeTimes);

        return changeTimes;
	}
}
