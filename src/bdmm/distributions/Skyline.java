package bdmm.distributions;

import beast.core.BEASTObject;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

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


    @Override
    public void initAndValidate() {

    }
}
