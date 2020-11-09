package bdmmprime.util.priors;

import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;

public class MigRateGLM extends CalculationNode implements Function {

    public Input<Function> covariateMatrixInput = new Input<>("covariateMatrix",
            "Matrix of covariate, e.g. numbers of flights between different locations",
            Input.Validate.REQUIRED);

    public Input<Function> scaleParamInput = new Input<>("scaleParam",
            "Scale parameter.", Input.Validate.REQUIRED);

    public Input<Function> errorInput = new Input<>("errorTerm",
            "Error term.", Input.Validate.OPTIONAL);

    Function covariateMatrix, scaleParam, errorTerm;

    @Override
    public void initAndValidate() {
        covariateMatrix = covariateMatrixInput.get();
        scaleParam = scaleParamInput.get();
        errorTerm = errorInput.get();
    }

    @Override
    public int getDimension() {
        return covariateMatrix.getDimension();
    }

    @Override
    public double getArrayValue(int i) {
        double lograte = scaleParam.getArrayValue()*covariateMatrix.getArrayValue(i);

        if (errorTerm!=null)
            lograte += errorTerm.getArrayValue(i);

        return Math.exp(lograte);
    }
}
