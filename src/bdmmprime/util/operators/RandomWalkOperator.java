package bdmmprime.util.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

public class RandomWalkOperator extends Operator {

    public Input<Double> windowSizeInput = new Input<>("windowSize",
            "Size of window used to select new value.", 1.0);

    public Input<RealParameter> parameterInput = new Input<>("parameter",
            "Parameter on which to operate", Input.Validate.REQUIRED);

    RealParameter param;
    double w;

    @Override
    public void initAndValidate() {
        param = parameterInput.get();
        w = windowSizeInput.get();
    }

    @Override
    public double proposal() {
        int pIdx = Randomizer.nextInt(param.getDimension());
        param.setValue(pIdx, param.getValue(pIdx) + (Randomizer.nextDouble()-0.5)*w);

        return 0;
    }
}
