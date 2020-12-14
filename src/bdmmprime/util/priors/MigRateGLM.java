package bdmmprime.util.priors;

import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;

public class MigRateGLM extends CalculationNode implements Function {

//    public Input<Function> covariateMatrixInput = new Input<>("covariateMatrix",
//            "Matrix of covariate, e.g. numbers of flights between different locations",
//            Input.Validate.REQUIRED);

    public Input<Function> covariateListInput = new Input<>("covariateList",
            "List of matrix of covariates, e.g. numbers of flights between different locations",
            Input.Validate.REQUIRED);

    public Input<Function> scaleParamInput = new Input<>("scaleParam",
            "Scale parameter.", Input.Validate.REQUIRED);

    public Input<Function> errorInput = new Input<>("errorTerm",
            "Error term.", Input.Validate.OPTIONAL);

    int covariateSize;

    Function covariateList, scaleParam, errorTerm;
//    Function covariateMatrix, scaleParam, errorTerm;

    @Override
    public void initAndValidate() {
//        covariateMatrix = covariateMatrixInput.get();
        covariateList = covariateListInput.get();
        scaleParam = scaleParamInput.get();
        errorTerm = errorInput.get();
        covariateSize = covariateList.getDimension()/(scaleParam.getDimension());
    }

    @Override
    public int getDimension() {
        return covariateSize;
//        return covariateList.getDimension()/(scaleParam.getDimension());
    }

    @Override
    public double getArrayValue(int i) {
        double lograte = 0;
//        double lograte = scaleParam.getArrayValue()*covariateMatrix.getArrayValue(i);
        for (int j = 0; j < scaleParam.getDimension(); j++){
            lograte += scaleParam.getArrayValue(j)*covariateList.getArrayValue(j * covariateSize + i);
        }

//        if (errorTerm!=null)
//            lograte += errorTerm.getArrayValue(i);

        return Math.exp(lograte);
    }
}
