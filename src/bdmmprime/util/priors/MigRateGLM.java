package bdmmprime.util.priors;

import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;

public class MigRateGLM extends CalculationNode implements Function {

    public Input<Function> flightMatrixInput = new Input<>("flightMatrix",
            "Matrix of numbers of flights between different locations",
            Input.Validate.REQUIRED);

    public Input<Function> scaleParamInput = new Input<>("scaleParam",
            "Scale parameter.", Input.Validate.REQUIRED);

    Function flightMatrix, scaleParam;

    @Override
    public void initAndValidate() {
        flightMatrix = flightMatrixInput.get();
        scaleParam = scaleParamInput.get();
    }

    @Override
    public int getDimension() {
        return flightMatrix.getDimension();
    }

    @Override
    public double getArrayValue(int i) {
        return scaleParam.getArrayValue()*flightMatrix.getArrayValue(i);
    }
}
