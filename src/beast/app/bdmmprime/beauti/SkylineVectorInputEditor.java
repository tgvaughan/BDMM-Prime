package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.SkylineVectorParameter;
import beast.app.beauti.BeautiDoc;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;

public class SkylineVectorInputEditor extends SkylineInputEditor {

    SkylineVectorParameter skylineVector;

    public SkylineVectorInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return SkylineVectorParameter.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption isExpandOption, boolean addButtons) {

        skylineVector = (SkylineVectorParameter) input.get();

        super.init(input, beastObject, itemNr, isExpandOption, addButtons);
    }

    @Override
    SkylineValuesTableModel getValuesTableModel() {
        return new SkylineVectorValuesTableModel(skylineVector.typeSetInput.get(), true, 1);
    }

    @Override
    void ensureParamsConsistent() {

        skylineVector.typeSetInput.get().initAndValidate();

        int nTypes = skylineVector.typeSetInput.get().getNTypes();
        int nIntervals = skylineVector.getChangeCount() + 1;

        RealParameter valuesParam = skylineVector.skylineValuesInput.get();
        int valuesPerInterval = valuesParam.getDimension() / nIntervals;

        if (valuesPerInterval == 1 || valuesPerInterval == nTypes) {
            skylineVector.initAndValidate();
            return;
        }

        StringBuilder valueBuilder = new StringBuilder();

        for (int interval=0; interval<nIntervals; interval++) {
            for (int typeIdx = 0; typeIdx < nTypes; typeIdx++) {
                valueBuilder.append(" ");

                if (typeIdx < valuesPerInterval)
                    valueBuilder.append(valuesParam.getValue(interval*valuesPerInterval + typeIdx));
                else
                    valueBuilder.append(valuesParam.getValue(interval*valuesPerInterval + (valuesPerInterval-1)));
            }
        }

        valuesParam.valuesInput.setValue(valueBuilder.toString(), valuesParam);
        valuesParam.initAndValidate();

        skylineVector.initAndValidate();
    }

    @Override
    String getChangeTimesParameterID() {
        int idx = skylineVector.getID().indexOf("SV");
        String prefix = skylineVector.getID().substring(0, idx);
        String suffix = skylineVector.getID().substring(idx+2);

        return prefix + "ChangeTimes" + suffix;
    }
}
