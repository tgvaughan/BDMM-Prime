package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.SkylineMatrixParameter;
import bdmmprime.parameterization.SkylineVectorParameter;
import beast.app.beauti.BeautiDoc;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import javax.swing.*;
import javax.swing.border.EtchedBorder;

public class SkylineMatrixInputEditor extends SkylineInputEditor {

    SkylineMatrixParameter skylineMatrix;

    public SkylineMatrixInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return SkylineMatrixParameter.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption isExpandOption, boolean addButtons) {

        skylineMatrix = (SkylineMatrixParameter) input.get();

        super.init(input, beastObject, itemNr, isExpandOption, addButtons);
    }

    @Override
    SkylineValuesTableModel getValuesTableModel() {
        return new SkylineMatrixValuesTableModel(skylineMatrix.typeSetInput.get(), true, 0);
    }

    @Override
    void ensureParamsConsistent() {

        skylineMatrix.typeSetInput.get().initAndValidate();

        int nTypes = skylineMatrix.typeSetInput.get().getNTypes();
        int nIntervals = skylineMatrix.getChangeCount() + 1;

        RealParameter valuesParam = skylineMatrix.skylineValuesInput.get();
        int valuesPerInterval = valuesParam.getDimension() / nIntervals;

        if (valuesPerInterval == 1 || valuesPerInterval == nTypes*(nTypes-1)) {
            skylineMatrix.initAndValidate();
            return;
        }

        StringBuilder valueBuilder = new StringBuilder();

        for (int interval=0; interval<nIntervals; interval++) {
            for (int type1Idx = 0; type1Idx < nTypes; type1Idx++) {
                for (int type2Idx = 0; type2Idx < nTypes; type2Idx++) {
                    if (type1Idx == type2Idx)
                        continue;

                    valueBuilder.append(" 0.0");
                }
            }
        }

        valuesParam.valuesInput.setValue(valueBuilder.toString(), valuesParam);
        valuesParam.initAndValidate();

        skylineMatrix.initAndValidate();
    }

    @Override
    void loadFromModel() {
        if (skylineMatrix.getNTypes()==1) {
            mainInputBox.setVisible(false);
            Box box = Box.createHorizontalBox();
            box.add(new JLabel("Insufficient types in model."));
            box.add(makeHorizontalFiller());
            add(box);
            return;
        }

        super.loadFromModel();
    }

    @Override
    String getChangeTimesParameterID() {
        int idx = skylineMatrix.getID().indexOf("SM");
        String prefix = skylineMatrix.getID().substring(0, idx);
        String suffix = skylineMatrix.getID().substring(idx+2);

        return prefix + "ChangeTimes" + suffix;
    }
}
