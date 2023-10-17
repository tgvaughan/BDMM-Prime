package bdmmprime.util.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Scale elements of a change time vector, maintaining order " +
        "of elements.")
public class ChangeTimeScaler extends Operator {

    public Input<RealParameter> parameterInput = new Input<>("parameter",
            "Change time parameter to scale.",
            Input.Validate.REQUIRED);

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Maximum scale factor used to scale an element.", 0.8);

    RealParameter param;
    double alphaMin, alphaMax;

    @Override
    public void initAndValidate() {
        param = parameterInput.get();
        alphaMin = Math.min(scaleFactorInput.get(), 1.0 / scaleFactorInput.get());
        alphaMax = 1.0/alphaMin;
    }

    @Override
    public double proposal() {

        double f = alphaMin + Randomizer.nextDouble()*(alphaMax - alphaMin);

        int idx = Randomizer.nextInt(param.getDimension());

        double lower = idx>0
                ? param.getValue(idx-1)
                : param.getLower();
        double upper = idx<param.getDimension() - 1
                ? param.getValue(idx+1)
                : param.getUpper();

        double x = param.getValue(idx)*f;

        if (param.getValue(idx)>0.0
                && x >= lower
                && x <= upper) {

            param.setValue(idx, x);
            return -Math.log(f);
        } else
            return Double.NEGATIVE_INFINITY;
    }

}
