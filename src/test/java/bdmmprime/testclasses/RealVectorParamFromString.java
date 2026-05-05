package bdmmprime.testclasses;

import beast.base.spec.domain.Domain;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;

public class RealVectorParamFromString<T extends Real> extends RealVectorParam<T> {

    public RealVectorParamFromString(String initString, T domain) {
        String[] spl = initString.split(" ");
        double[] valArray = new double[spl.length];
        for (int i=0; i<valArray.length; i++)
            valArray[i] = Double.parseDouble(spl[i]);

        super(valArray, domain);
    }
}
